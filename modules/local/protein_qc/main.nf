process PROTEIN_QC {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle'), path(cells_qc)
    path markers

    output:
    tuple val(meta), path("protein"), emit: dir

    script:
    def m = markers ? "--markers ${markers}" : ''
    """
    xqc protein --bundle bundle --outdir protein --cells-qc ${cells_qc} ${m}
    """
}
