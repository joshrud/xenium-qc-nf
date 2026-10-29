process CELL_QC {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle')
    path markers

    output:
    tuple val(meta), path("cells"), emit: dir
    tuple val(meta), path("cells/cells_qc.parquet"), emit: cells_qc
    tuple val(meta), path("cells/filtered.h5ad"), emit: h5ad

    script:
    def m = markers ? "--markers ${markers}" : ''
    """
    export POLARS_MAX_THREADS=${task.cpus}
    xqc cells --bundle bundle --outdir cells ${m} --min-counts ${params.min_counts} \\
        --min-cells-per-gene ${params.min_cells_per_gene} --mad-k ${params.mad_k} --target-sum ${params.target_sum}
    """
}
