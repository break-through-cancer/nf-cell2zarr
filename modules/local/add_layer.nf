process ADD_LAYER {
    tag   "${meta.id}:${layer}"
    label 'process_add_large'

    stageInMode 'copy'

    input:
    tuple val(meta), path(zarr_store, stageAs: 'store.zarr'), path(h5ad), val(layer)

    output:
    tuple val(meta), path("groups/*"), emit: groups
    path "add_layer_${layer}.log",     emit: log

    script:
    """
    mkdir -p scratch

    /app/.venv/bin/cell2zarr add ${h5ad} store.zarr \\
        --key layers/${layer} \\
        --no-consolidate \\
        --dtype ${params.dtype} \\
        --temp-dir scratch \\
        --log-file add_layer_${layer}.log

    rm -rf scratch

    mkdir -p groups
    cp -R store.zarr/layers groups/layers
    """

    stub:
    """
    mkdir -p groups/layers/${layer}
    echo '{"zarr_format": 3, "node_type": "group", "attributes": {}}' > groups/layers/zarr.json
    echo '{"zarr_format": 3, "node_type": "group", "attributes": {}}' > groups/layers/${layer}/zarr.json
    touch add_layer_${layer}.log
    """
}
