process ALIGNMENT_QC {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle'), path(keypoints), path(alignment), path(reference)

    output:
    tuple val(meta), path("alignment"), emit: dir

    script:
    def kp  = keypoints ? "--keypoints ${keypoints}" : ''
    def aln = alignment ? "--he-alignment ${alignment}" : ''
    def ref = reference ? "--reference-alignment ${reference}" : ''
    """
    export POLARS_MAX_THREADS=${task.cpus}
    xqc alignment --bundle bundle --outdir alignment ${kp} ${aln} ${ref}
    """
}
