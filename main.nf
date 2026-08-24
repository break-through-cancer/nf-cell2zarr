#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONVERT }      from './modules/local/convert'
include { ADD_KEYS }     from './modules/local/add_keys'
include { ADD_LAYER }    from './modules/local/add_layer'
include { CONSOLIDATE }  from './modules/local/consolidate'
include { BUILD_STRATA } from './modules/local/build_strata'

workflow NF_CELL2ZARR {

    if (!params.input) {
        error "params.input is required: the path to an .h5ad file"
    }

    def h5ad = file(params.input, checkIfExists: true)

    if (!h5ad.name.endsWith('.h5ad')) {
        error "nf-cell2zarr accepts .h5ad input only; got '${h5ad.name}'. " +
              "Loom and Seurat files in a Single-Cell Dataset are not supported."
    }

    def cheap_keys = (params.add_keys   ?: '').tokenize(',')*.trim().findAll { it }
    def layer_keys = (params.add_layers ?: '').tokenize(',')*.trim().findAll { it }

    ch_input = Channel.of( tuple([id: h5ad.simpleName], h5ad) )

    CONVERT(ch_input)

    // One task for all cheap keys: they are fast to write, so staging the
    // base store dominates and splitting them would pay that cost per key.
    ch_key_groups = Channel.empty()
    if (cheap_keys) {
        ADD_KEYS(
            CONVERT.out.zarr.map { meta, zarr -> tuple(meta, zarr, h5ad, cheap_keys) }
        )
        ch_key_groups = ADD_KEYS.out.groups
    }

    // Genuine fan-out: layers go through the two-phase path and are expensive
    // enough that a parallel task pays for its staging.
    ch_layer_groups = Channel.empty()
    if (layer_keys) {
        ADD_LAYER(
            CONVERT.out.zarr
                .combine(Channel.fromList(layer_keys))
                .map { meta, zarr, layer -> tuple(meta, zarr, h5ad, layer) }
        )
        ch_layer_groups = ADD_LAYER.out.groups
    }

    ch_groups = ch_key_groups
        .mix(ch_layer_groups)
        .map { meta, groups -> groups }
        .flatten()
        .collect()
        .ifEmpty( [] )

    CONSOLIDATE(CONVERT.out.zarr, ch_groups)

    // build-strata needs obs column names nobody can know without having seen
    // the dataset, so it runs only when the user supplied them.
    if (params.atomic_axes) {
        BUILD_STRATA(CONSOLIDATE.out.zarr)
    }
}

workflow {
    NF_CELL2ZARR()
}
