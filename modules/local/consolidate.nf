process CONSOLIDATE {
    tag   "${meta.id}"
    label 'process_add_small'

    stageInMode 'copy'

    input:
    tuple val(meta), path(zarr_store, stageAs: 'store.zarr')
    path groups, stageAs: 'groups/*'

    output:
    tuple val(meta), path("${meta.id}.zarr"), emit: zarr
    path "consolidate.log",                   emit: log

    script:
    """
    # Merge each emitted group into the store. Several ADD_LAYER tasks each
    # emit a directory named 'layers', so copy directory *contents* rather
    # than replacing the directory, or all but one layer would be lost.
    for group in groups/*; do
        [ -e "\$group" ] || continue
        name="\$(basename "\$group")"
        mkdir -p "store.zarr/\$name"
        cp -R "\$group/." "store.zarr/\$name/"
    done

    /app/.venv/bin/cell2zarr consolidate store.zarr --log-file consolidate.log

    mv store.zarr ${meta.id}.zarr
    """

    stub:
    """
    mkdir -p ${meta.id}.zarr
    echo '{"zarr_format": 3, "node_type": "group", "consolidated_metadata": {"metadata": {}}}' > ${meta.id}.zarr/zarr.json
    touch consolidate.log
    """
}
