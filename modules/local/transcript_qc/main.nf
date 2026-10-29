process TRANSCRIPT_QC {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle')

    output:
    tuple val(meta), path("transcripts"), emit: dir

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus}
    xqc transcripts --bundle bundle --outdir transcripts --bin-um ${params.bin_um} --seam-window ${params.seam_window}
    """
}
