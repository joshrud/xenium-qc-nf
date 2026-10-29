process RUN_METRICS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle')

    output:
    tuple val(meta), path("run_metrics"), emit: dir

    script:
    """
    xqc run-metrics --bundle bundle --outdir run_metrics
    """
}
