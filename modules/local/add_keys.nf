process ADD_KEYS {
    tag   "${meta.id}"
    label 'process_add_small'

    stageInMode 'copy'

    input:
    tuple val(meta), path(zarr_store, stageAs: 'store.zarr'), path(h5ad), val(keys)

    output:
    tuple val(meta), path("groups/*"), emit: groups
    path "add_keys.log",               emit: log

    script:
    def key_list = keys.join(' ')
    """
    for key in ${key_list}; do
        /app/.venv/bin/cell2zarr add ${h5ad} store.zarr \\
            --key "\$key" \\
            --no-consolidate \\
            --log-file add_keys.log
    done

    mkdir -p groups
    for key in ${key_list}; do
        top="\${key%%/*}"
        if [ ! -e "groups/\$top" ]; then
            cp -R "store.zarr/\$top" "groups/\$top"
        fi
    done
    """

    stub:
    def key_list = keys.join(' ')
    """
    mkdir -p groups
    for key in ${key_list}; do
        top="\${key%%/*}"
        mkdir -p "groups/\$top"
        echo '{"zarr_format": 3, "node_type": "group", "attributes": {}}' > "groups/\$top/zarr.json"
    done
    touch add_keys.log
    """
}
