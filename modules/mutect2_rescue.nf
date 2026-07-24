// modules/mutect2_rescue.nf
//
// Post-hoc rescue of contamination-flagged variants after FilterMutectCalls.
// Reads mutect2.filtered.vcf.gz, applies 5-criterion rescue logic, and emits
// mutect2.filtered.rescued.vcf.gz (bgzipped + tabix indexed).

nextflow.enable.dsl=2

params.outdir    = params.outdir    ?: "results"
params.scriptdir = params.scriptdir ?: "${workflow.projectDir}/scripts"

process MUTECT2_RESCUE {
    tag { sample }
    publishDir "${params.outdir}/${sample}/somatic", mode: 'copy'

    input:
        tuple val(sample), path(filtered_vcf)
        path rescue_script

    output:
        tuple val(sample), path("${sample}.mutect2.filtered.rescued.vcf.gz"),
                           path("${sample}.mutect2.filtered.rescued.vcf.gz.tbi")

    script:
    """
    python3 ${rescue_script} ${filtered_vcf} rescued.vcf

    bgzip -c rescued.vcf > ${sample}.mutect2.filtered.rescued.vcf.gz
    tabix -p vcf ${sample}.mutect2.filtered.rescued.vcf.gz
    """
}
