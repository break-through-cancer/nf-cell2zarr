process CONVERT {
    tag   "${meta.id}"
    label 'process_convert'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}.zarr"), emit: zarr
    path "convert.log",                       emit: log

    script:
    def two_phase = params.two_phase           ? '--two-phase'                                        : ''
    def n_top     = params.n_top_genes         ? "--n-top-genes ${params.n_top_genes}"                : ''
    def norm      = params.normalize           ? '--normalize'                                        : ''
    def enc       = params.encoding_config     ? "--encoding-config ${params.encoding_config}"        : ''
    def obsm_cs   = params.obsm_cell_chunk_size ? "--obsm-cell-chunk-size ${params.obsm_cell_chunk_size}" : ''
    """
    mkdir -p scratch

    /app/.venv/bin/cell2zarr convert ${h5ad} ${meta.id}.zarr \\
        ${two_phase} \\
        --cell-chunk-size ${params.cell_chunk_size} \\
        --var-chunk-size ${params.var_chunk_size} \\
        --shard-size ${params.shard_size} \\
        --dtype ${params.dtype} \\
        --sparse-format ${params.sparse_format} \\
        ${n_top} ${norm} ${enc} ${obsm_cs} \\
        --temp-dir scratch \\
        --log-file convert.log

    rm -rf scratch
    """

    stub:
    """
    mkdir -p ${meta.id}.zarr
    cat > ${meta.id}.zarr/zarr.json <<'JSON'
    {"zarr_format": 3, "node_type": "group", "attributes": {"encoding-type": "anndata"}}
    JSON
    touch convert.log
    """
}
