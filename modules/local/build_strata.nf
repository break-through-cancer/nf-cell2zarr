process BUILD_STRATA {
    tag   "${meta.id}"
    label 'process_strata'

    stageInMode 'copy'

    input:
    tuple val(meta), path(zarr_store, stageAs: 'store.zarr')

    output:
    tuple val(meta), path("${meta.id}.zarr"), emit: zarr
    path "build_strata.log",                  emit: log

    script:
    def axes = (params.atomic_axes ?: '')
        .tokenize(',')*.trim()
        .findAll { it }
        .collect { "--atomic-axes ${it}" }
        .join(' ')

    def coarse = (params.coarse ?: '')
        .tokenize(';')*.trim()
        .findAll { it }
        .collect { "--coarse ${it}" }
        .join(' ')

    def max_card = params.max_atomic_cardinality
        ? "--max-atomic-cardinality ${params.max_atomic_cardinality}"
        : ''
    """
    # build-strata has no --log-file, so its output is captured with tee. A pipeline
    # returns the LAST command's status, and Nextflow's default shell is `bash -ue`
    # without pipefail — without this, a failing build-strata would be masked by tee
    # exiting 0 and the task would be reported as successful.
    set -o pipefail

    /app/.venv/bin/cell2zarr build-strata store.zarr \\
        ${axes} ${coarse} ${max_card} \\
        2>&1 | tee build_strata.log

    mv store.zarr ${meta.id}.zarr
    """

    stub:
    """
    mkdir -p ${meta.id}.zarr/uns/strata/atomic
    echo '{"zarr_format": 3, "node_type": "group", "attributes": {}}' > ${meta.id}.zarr/zarr.json
    touch build_strata.log
    """
}
