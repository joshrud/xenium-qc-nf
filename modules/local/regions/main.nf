process REGIONS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle'), path(h5ad), path(annotation), path(alignment)
    path markers

    output:
    tuple val(meta), path("regions"), emit: dir

    script:
    def m   = markers ? "--markers ${markers}" : ''
    def ann = annotation ? "--annotation ${annotation} --annotation-units ${params.annotation_units}" : ''
    def aln = alignment ? "--he-alignment ${alignment}" : ''
    """
    xqc regions --bundle bundle --outdir regions --filtered-h5ad ${h5ad} ${m} ${ann} ${aln} \\
        --leiden-resolution ${params.leiden_resolution}
    """
}
