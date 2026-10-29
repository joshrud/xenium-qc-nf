process QC_REPORT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(step_dirs), path(multiqc_html)

    output:
    tuple val(meta), path("report/${meta.id}_xenium_qc_report.html"), emit: html
    tuple val(meta), path("report"), emit: dir

    script:
    def mqc = multiqc_html ? "--multiqc-report ../multiqc/multiqc_report.html" : ''
    def p = params.findAll { k, v -> !(v instanceof Map) && !k.contains('-') }
                  .collectEntries { k, v -> [k, v == null ? null : v.toString()] }
    def pjson = groovy.json.JsonOutput.toJson(p)
    """
    cat <<'PARAMS_EOF' > params.json
${pjson}
PARAMS_EOF
    xqc report --step-dirs ${step_dirs} --outdir report --sample ${meta.id} ${mqc} --params-json params.json
    """
}
