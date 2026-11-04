process PREPARE {
    tag   "${meta.id}"
    label 'process_prepare'

    container params.prepare_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${meta.id}.prepared.h5ad"), emit: h5ad
    path "provenance.json",                            emit: provenance
    path "prepare.log",                                emit: log

    script:
    def norm = params.normalize ? '--normalize' : ''
    // Regexes are single-quoted so the shell leaves their backslashes alone;
    // Groovy interpolation does not touch backslashes inside a param value.
    """
    prepare_h5ad.py \\
        --input ${h5ad} \\
        --output ${meta.id}.prepared.h5ad \\
        --provenance provenance.json \\
        --n-hvg ${params.n_hvg} \\
        --n-pcs ${params.n_pcs} \\
        --n-neighbors ${params.n_neighbors} \\
        --leiden-resolution ${params.leiden_resolution} \\
        --library-column '${params.library_column}' \\
        --site-regex '${params.site_regex}' \\
        --patient-regex '${params.patient_regex}' \\
        --patient-prefix '${params.patient_prefix}' \\
        --tcr-column '${params.tcr_column}' \\
        --atomic-axes '${params.atomic_axes ?: ''}' \\
        ${norm} \\
        --threads ${task.cpus} \\
        --pipeline-revision '${workflow.commitId ?: workflow.revision ?: 'local'}' \\
        --container '${task.container}' \\
        2> >(tee prepare.log >&2)
    """

    stub:
    """
    touch ${meta.id}.prepared.h5ad
    echo '{"stub": true}' > provenance.json
    touch prepare.log
    """
}
