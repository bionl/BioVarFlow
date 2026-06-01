// modules/db_qc_export.nf
// ───────────────────────────────────────────────────────────────────────────
// DB QC Export — per-sample QC metrics for downstream ingestion decisions.
//
// Reuses Sarek's existing QC outputs (samtools stats OR markdup metrics
// for alignment; mosdepth for coverage). Only runs bcftools stats on the
// consensus VCF.
//
// This module is advisory only — it never blocks or filters samples.
// ───────────────────────────────────────────────────────────────────────────
nextflow.enable.dsl=2

params.outdir    = params.outdir    ?: "results"
params.scriptdir = params.scriptdir ?: "${workflow.projectDir}/scripts"

/********************  PROCESSES  ********************/

process DB_QC_VCFSTATS {
    tag { "${meta.sample}" }

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta),
            path("${meta.sample}.bcftools_stats.txt"),
            path("${meta.sample}.callers.txt")

    script:
    def sample = meta.sample
    """
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        bcftools index -t ${vcf} 2>/dev/null || true
    fi

    bcftools stats ${vcf} > ${sample}.bcftools_stats.txt

    if bcftools view -h ${vcf} | grep -q 'ID=CALLERS'; then
        bcftools query -f '%INFO/CALLERS\\n' ${vcf} > ${sample}.callers.txt
    else
        touch ${sample}.callers.txt
    fi
    """
}

process DB_QC_EXPORT_JSON {
    tag { "${meta.sample}" }
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        tuple val(meta),
            path(alignment_file),
            val(alignment_type),
            path(mosdepth_summary),
            path(mosdepth_global_dist),
            path(bcftools_stats),
            path(callers)
        each path(script)

    output:
        tuple val(meta), path("${meta.sample}.qc.json"), path("${meta.sample}.qc.tsv")

    script:
    def sample = meta.sample
    def assay  = meta.assay ?: "NA"
    def aln_arg = alignment_type == 'samtools_stats'
        ? "--samtools-stats ${alignment_file}"
        : "--markdup-metrics ${alignment_file}"
    """
    python ${script} \\
        ${aln_arg} \\
        --mosdepth-summary ${mosdepth_summary} \\
        --mosdepth-global-dist ${mosdepth_global_dist} \\
        --bcftools-stats ${bcftools_stats} \\
        --callers ${callers} \\
        --sample ${sample} \\
        --assay ${assay} \\
        --outdir .
    """
}

process DB_QC_AGGREGATE {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        path(tsvs)

    output:
        path("qc_summary.tsv")

    script:
    """
    head -1 \$(ls *.qc.tsv | head -1) > qc_summary.tsv
    tail -q -n +2 *.qc.tsv >> qc_summary.tsv
    """
}

/********************  WORKFLOW  ********************/

workflow DB_QC_EXPORT {
    take:
        vcf_ch                // tuple(meta, vcf) — consensus VCF
        alignment_qc_ch       // tuple(meta, file, type) — samtools_stats or markdup_metrics
        mosdepth_summary_ch   // tuple(meta, summary) — from Sarek
        mosdepth_dist_ch      // tuple(meta, dist) — from Sarek

    main:
        qc_script_ch = Channel
            .fromPath("${params.scriptdir}/compute_qc_metrics.py")
            .first()

        // Only variant stats needs to run — alignment & coverage come from Sarek
        DB_QC_VCFSTATS(vcf_ch)

        // Join all outputs by meta key → single tuple per sample
        // alignment_qc_ch carries a type tag through the join
        qc_input_ch = alignment_qc_ch                                 // (meta, file, type)
            .join(mosdepth_summary_ch)                                 // + (summary)
            .join(mosdepth_dist_ch)                                    // + (dist)
            .join(DB_QC_VCFSTATS.out)                                  // + (bcfstats, callers)
        // result: (meta, alignment_file, alignment_type, summary, dist, bcfstats, callers)

        DB_QC_EXPORT_JSON(qc_input_ch, qc_script_ch)

        // Aggregate all per-sample TSVs into one summary
        all_tsvs = DB_QC_EXPORT_JSON.out
            .map { meta, json, tsv -> tsv }
            .collect()

        DB_QC_AGGREGATE(all_tsvs)

    emit:
        qc_json    = DB_QC_EXPORT_JSON.out   // tuple(meta, json, tsv) per sample
        qc_summary = DB_QC_AGGREGATE.out     // single aggregated TSV
}
