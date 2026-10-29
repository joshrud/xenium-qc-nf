process MORPHOLOGY_QC {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle'), path(he_image), path(alignment)
    path markers

    output:
    tuple val(meta), path("morphology"), emit: dir

    script:
    def m = markers ? "--markers ${markers}" : ''
    def he = (he_image && alignment) ? "--he-image ${he_image} --he-alignment ${alignment}" : ''
    """
    export POLARS_MAX_THREADS=${task.cpus}
    xqc morphology --bundle bundle --outdir morphology ${m} ${he}
    """
}
