// modules/vep.nf
nextflow.enable.dsl=2

// -------- Parameters (used by processes) --------
params.bed         = params.bed         ?: "${workflow.projectDir}/data/ACMG_SF_MANE_exons_50bp.bed"
// Somatic reporting panel — restricts the somatic (Mutect2) call set to the
// genes we report on, the same way params.bed scopes the germline ACMG SF set.
params.somatic_bed = params.somatic_bed ?: "${workflow.projectDir}/data/Somatic_125genes_MANE_50bp.bed"
params.outdir      = params.outdir      ?: "results"
params.scriptdir   = params.scriptdir   ?: "${workflow.projectDir}/scripts"
params.template_dir= params.template_dir?: "${workflow.projectDir}/scripts/template-files"

params.run_vep     = params.run_vep     ?: true
params.min_dp   = params.min_dp   ?: 10
params.min_qual = params.min_qual ?: 10
// VEP resource params expected from main/config:
// params.vep_fasta, params.revel_vcf, params.alpha_missense_vcf, params.clinvar_vcf

// Reference for `bcftools norm -f` on BOTH the germline and somatic paths.
//
// Resolved lazily, not at module load: params.fasta is set by
// external/sarek/main.nf when it is included, and the igenomes map is the
// reliable fallback since config-parse time always populates it.
//
// Deliberately NOT params.ref_fasta -- main.nf defaults that to params.vep_fasta,
// which is the Ensembl VEP reference whose contigs are named 1/2/.../MT. Passing
// it here fails on every record, since the BAMs are chr-prefixed GATK.
def resolveNormFasta() {
    def f = params.norm_fasta ?: params.fasta
    if (!f && params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        f = params.genomes[params.genome].fasta
    }
    if (!f) {
        error "❌ No reference for bcftools norm.\n" +
              "   Tried --norm_fasta, params.fasta, params.genomes[${params.genome}].fasta.\n" +
              "   genome=${params.genome}  genomes_loaded=${params.genomes ? params.genomes.size() : 0}\n" +
              "   Set --norm_fasta to the GATK assembly the BAMs were aligned against.\n" +
              "   Do NOT point this at the Ensembl VEP fasta -- its contigs are 1/2/.../MT."
    }
    return f
}


/********************  PROCESSES (unchanged logic, publish to per-sample dirs)  ********************/

process BedFilterVCF {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  input:
    tuple val(meta), path(vcf)
    path bed
  output:
    tuple val(meta), path("${meta.sample}.bed_filtered.vcf.gz")
  script:
    def sample = meta.sample
  """
  tabix -p vcf $vcf || bcftools index -t $vcf
  bcftools view -f PASS -R $bed $vcf -Oz -o ${sample}.bed_filtered.vcf.gz
  tabix -p vcf ${sample}.bed_filtered.vcf.gz
  """
}

process NormalizeVCF {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  input:
    tuple val(meta), path(vcf)
    path fasta
    path fai
  output:
    tuple val(meta), path("${meta.sample}.normalized.vcf.gz")
  script:
    def sample = meta.sample
  """
  # -m -any  splits multiallelic sites into one record per ALT allele.
  # -f       left-aligns indels and trims shared flanking bases. Without it,
  #          alleles stay non-parsimonious (e.g. TAA>GAA instead of T>G), which
  #          breaks position-keyed annotation lookups -- SpliceAI silently
  #          returns nothing for the affected records, and ClinVar misses indels.
  #
  # The reference MUST be the one the BAMs were aligned to (chr-prefixed GATK
  # assembly), never the Ensembl VEP fasta -- see resolveNormFasta().
  #
  # No -c override: bcftools exits on a REF mismatch by default, which is what
  # we want -- a mismatch means the wrong reference, and continuing would
  # silently corrupt allele representations.
  bcftools norm -m -any -f $fasta $vcf -Oz -o ${sample}.normalized.vcf.gz
  tabix -p vcf ${sample}.normalized.vcf.gz
  """
}

process FilterVCF {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  input:
    tuple val(meta), path(vcf)
  output:
    tuple val(meta), path("${meta.sample}.filtered.vcf.gz")
  script:
    def sample = meta.sample
  """
  bcftools view -i 'FORMAT/DP >= ${params.min_dp} && QUAL >= ${params.min_qual}' $vcf -Oz -o ${sample}.filtered.vcf.gz
  tabix -p vcf ${sample}.filtered.vcf.gz
  """
}


process BedFilterBAM {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  input:
    tuple val(meta), path(vcf), path(bam)
    path bed
  output:
    tuple val(meta), path("${meta.sample}.bed_filtered.bam"), path("${meta.sample}.bed_filtered.bam.bai")
  script:
    def sample = meta.sample
  """
  samtools view -L $bed -b -@ 16 $bam -o tmp.bam
  samtools sort -o ${sample}.bed_filtered.bam tmp.bam
  samtools index ${sample}.bed_filtered.bam
  """
}

process CoverageSummary {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy', pattern: "*_coverage_summary.sorted.txt"
  input:
    tuple val(meta), path(bam)
    path bed
  output:
    tuple val(meta), path("${meta.sample}_coverage_summary.sorted.txt")
  script:
    def sample = meta.sample
  """
  bedtools coverage -a $bed -b $bam -d > ${sample}_coverage_per_base.txt
  awk '{
    key=\$1":"\$2"-"\$3; total[key]++
    if(\$5>=20) c20[key]++; if(\$5>=30) c30[key]++; if(\$5>=50) c50[key]++; if(\$5>=100) c100[key]++
  } END {
    for (k in total)
      printf "%s\\t>=20x:%.2f%%\\t>=30x:%.2f%%\\t>=50x:%.2f%%\\t>=100x:%.2f%%\\n", k,(c20[k]/total[k])*100,(c30[k]/total[k])*100,(c50[k]/total[k])*100,(c100[k]/total[k])*100
  }' ${sample}_coverage_per_base.txt > ${sample}_coverage_summary.txt
  sort -t: -k1,1 -k2,2n ${sample}_coverage_summary.txt > ${sample}_coverage_summary.sorted.txt
  """
}

process R1R2Ratio {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  input:
    tuple val(meta), path(bam), path(bai)
    path bed
  output:
    tuple val(meta), path("${meta.sample}_r1r2_per_exon.tsv")
  script:
    def sample = meta.sample
  """
  while read chrom start end ref_name; do
    region="\${chrom}:\${start}-\${end}"
    counts=\$(samtools view -F 0x904 $bam "\$region" | \
      awk '{flag=\$2; if(and(flag,64)) r1++; if(and(flag,128)) r2++} END {if(r1+r2>0) printf("%d\\t%d\\t%.3f\\n", r1, r2, r1/(r1+r2)); else print "0\\t0\\tNA"}')
    echo -e "\${chrom}\\t\${start}\\t\${end}\\t\${ref_name}\\t\${counts}"
  done < $bed > ${sample}_r1r2_per_exon.tsv
  """
}

process ForwardReverseRatio {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  input:
    tuple val(meta), path(bam), path(bai)
    path bed
  output:
    tuple val(meta), path("${meta.sample}_frstrand_per_exon.tsv")
  script:
    def sample = meta.sample
  """
  while read chrom start end ref_name; do
    region="\${chrom}:\${start}-\${end}"
    counts=\$(samtools view -F 0x904 $bam "\$region" | \
      awk '{flag=\$2; if(and(flag,16)) rev++; else fwd++} END {if(fwd+rev>0){frac=rev/(fwd+rev); bal=(fwd/(fwd+rev)<rev/(fwd+rev)?fwd/(fwd+rev):rev/(fwd+rev)); printf("%d\\t%d\\t%.3f\\t%.3f\\n",fwd,rev,frac,bal)} else print "0\\t0\\tNA\\tNA"}')
    echo -e "\${chrom}\\t\${start}\\t\${end}\\t\${ref_name}\\t\${counts}"
  done < $bed > ${sample}_frstrand_per_exon.tsv
  """
}

process SamtoolsFlagstat {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
  input:
    tuple val(meta), path(bam)
  output:
    tuple val(meta), path("${meta.sample}_flagstat.txt")
  script:
    def sample = meta.sample
  """
  samtools flagstat $bam > ${sample}_flagstat.txt
  """
}

process SamtoolsStats {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
  input:
    tuple val(meta), path(bam)
  output:
    tuple val(meta), path("${meta.sample}_stats.txt")
  script:
    def sample = meta.sample
  """
  samtools stats $bam > ${sample}_stats.txt
  """
}
process MosdepthRun {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay 
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'

  input:
    tuple val(meta), path(bam), path(bai)
    path bed

  output:
    tuple val(meta),
      path("${meta.sample}.mosdepth.summary.txt"),
      path("${meta.sample}.thresholds.bed.gz"),
      path("${meta.sample}.quantized.bed.gz")
      //path("${sample}_coverage_summary.overall.txt")

  script:
    def sample = meta.sample
  """
  echo "[\$(date -Is)] Starting mosdepth for ${sample}" >&2
  set -euo pipefail
  cp $bam ./${sample}.bam
  cp $bai ./${sample}.bam.bai
  export MOSDEPTH_Q0=LT20
  export MOSDEPTH_Q1=GE20_LT30
  export MOSDEPTH_Q2=GE30

  mosdepth --no-per-base --by $bed --thresholds 10,20,30,50,100 --quantize 0:20:30: --fast-mode $sample $bam
  echo "[\$(date -Is)] Finished mosdepth for ${sample}" >&2
  #python ${params.scriptdir}/summarize_mosdepth.py \
  #  --prefix $sample \
  #  --summary ${sample}.mosdepth.summary.txt \
  #  --thresholds ${sample}.thresholds.bed.gz \
  #  --out ${sample}_coverage_summary.overall.txt
  """
}

process CoverageGapsAnnotation {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay 
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'

  input:
    tuple val(meta),
      path("${meta.sample}.quantized.bed.gz"),
      path("${meta.sample}.thresholds.bed.gz")
    path bed

  output:
    tuple val(meta),
      path("${meta.sample}.acmg_gaps_lt20.annot.bed"),
      path("${meta.sample}.acmg_gaps_lt30.annot.bed")

  script:
    def sample = meta.sample
  """
  zcat ${sample}.quantized.bed.gz | awk '\$4=="LT20"' \
    | bedtools intersect -wa -a - -b $bed \
    | bedtools sort -i - \
    | bedtools merge -i - > ${sample}.acmg_gaps_lt20.bed

  zcat ${sample}.quantized.bed.gz | awk '\$4=="LT20" || \$4=="GE20_LT30"' \
    | bedtools intersect -wa -a - -b $bed \
    | bedtools sort -i - \
    | bedtools merge -i - > ${sample}.acmg_gaps_lt30.bed

  bedtools intersect -wao -a ${sample}.acmg_gaps_lt20.bed -b $bed \
    | awk 'BEGIN{OFS="\\t"}{ label = (\$7 != "" && \$7 != ".") ? \$7 : \$4 ":" \$5 "-" \$6; print \$1,\$2,\$3,label}' \
    > ${sample}.acmg_gaps_lt20.annot.bed

  bedtools intersect -wao -a ${sample}.acmg_gaps_lt30.bed -b $bed \
    | awk 'BEGIN{OFS="\\t"}{ label = (\$7 != "" && \$7 != ".") ? \$7 : \$4 ":" \$5 "-" \$6; print \$1,\$2,\$3,label}' \
    > ${sample}.acmg_gaps_lt30.annot.bed
  """
}
//process MosdepthCoverage {
//  tag "$sample"
//  publishDir "${params.outdir}/${sample}/qc", mode: 'copy'
//  input:
//    tuple val(sample), path(bam), path(bai)
//    path bed
//  output:
//    tuple val(sample),
//      path("${sample}.mosdepth.summary.txt"),
//      path("${sample}.regions.bed.gz"),
//      path("${sample}.thresholds.bed.gz"),
//      path("${sample}.quantized.bed.gz"),
//      path("${sample}_coverage_summary.overall.txt"),
//      path("${sample}.acmg_gaps_lt20.bed"),
//      path("${sample}.acmg_gaps_lt30.bed"),
//      path("${sample}.acmg_gaps_lt20.annot.bed"),
//      path("${sample}.acmg_gaps_lt30.annot.bed")
//  script:
//  """
//  set -euo pipefail
//  export MOSDEPTH_Q0=LT20
//  export MOSDEPTH_Q1=GE20_LT30
//  export MOSDEPTH_Q2=GE30
//
//  mosdepth --no-per-base --by $bed --thresholds 10,20,30,50,100 --quantize 0:20:30: --fast-mode $sample $bam
//  python ${params.scriptdir}/summarize_mosdepth.py \
//    --prefix $sample \
//    --summary ${sample}.mosdepth.summary.txt \
//    --thresholds ${sample}.thresholds.bed.gz \
//    --out ${sample}_coverage_summary.overall.txt
//
//  zcat ${sample}.quantized.bed.gz | awk '\$4=="LT20"' | bedtools intersect -wa -a - -b $bed | bedtools sort -i - | bedtools merge -i - > ${sample}.acmg_gaps_lt20.bed
//  zcat ${sample}.quantized.bed.gz | awk '\$4=="LT20" || \$4=="GE20_LT30"' | bedtools intersect -wa -a - -b $bed | bedtools sort -i - | bedtools merge -i - > ${sample}.acmg_gaps_lt30.bed
//
//  bedtools intersect -wao -a ${sample}.acmg_gaps_lt20.bed -b $bed | \
//    awk 'BEGIN{OFS="\t"}{ label = (\$7 != "" && \$7 != ".") ? \$7 : \$4 ":" \$5 "-" \$6; print \$1,\$2,\$3,label}' > ${sample}.acmg_gaps_lt20.annot.bed
//
//  bedtools intersect -wao -a ${sample}.acmg_gaps_lt30.bed -b $bed | \
//    awk 'BEGIN{OFS="\t"}{ label = (\$7 != "" && \$7 != ".") ? \$7 : \$4 ":" \$5 "-" \$6; print \$1,\$2,\$3,label}' > ${sample}.acmg_gaps_lt30.annot.bed
//  """
//}

process SexCheck {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay 
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
  input:
    tuple val(meta), path(bam)
  output:
    tuple val(meta), path("${meta.sample}_sex_check.txt")
  script:
    def sample = meta.sample
  """
  x_depth=\$(samtools idxstats $bam | awk '\$1=="X"{print \$3}')
  y_depth=\$(samtools idxstats $bam | awk '\$1=="Y"{print \$3}')
  if [ "\$y_depth" -gt 0 ]; then ratio=\$(echo "scale=3; \$x_depth/\$y_depth" | bc); else ratio="NA"; fi
  if [ "\$ratio" != "NA" ]; then
    if (( \$(echo "\$ratio > 4" | bc -l) )); then sex="Female"; else sex="Male"; fi
  else sex="Unknown"; fi
  {
    echo "Sample: $sample"
    echo "X chromosome depth: \$x_depth"
    echo "Y chromosome depth: \$y_depth"
    echo "X/Y ratio: \$ratio"
    echo "Predicted sex: \$sex"
  } > ${sample}_sex_check.txt
  """
}

process BcftoolsStats {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay 
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
  input:
    tuple val(meta), path(vcf)
  output:
    tuple val(meta), path("${meta.sample}_bcftools_stats.txt")
  script:
    def sample = meta.sample
  """
  bcftools stats $vcf > ${sample}_bcftools_stats.txt
  """
}

process VEP_Annotate {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay 
  publishDir "${params.outdir}/${meta.sample}/vcf", mode: 'copy'
  input:
    tuple val(meta), path(vcf)
    path vep_cache
    path vep_fasta
    path vep_fasta_fai
    path revel_vcf
    path revel_vcf_tbi
    path alpha_missense_vcf
    path alpha_missense_vcf_tbi
    path clinvar_vcf
    path clinvar_vcf_tbi
    path spliceai_snv_vcf
    path spliceai_snv_vcf_tbi
    path spliceai_indel_vcf
    path spliceai_indel_vcf_tbi
    path bayesdel_vcf
    path bayesdel_vcf_tbi
    path vep_plugins
  output:
    tuple val(meta), path("${meta.sample}.vep.vcf")
  script:
    def sample = meta.sample
  """
  set -euo pipefail
  if [[ "$vcf" == *.vcf.gz ]]; then gunzip -c "$vcf" > INPUT_FOR_VEP.vcf; else cp "$vcf" INPUT_FOR_VEP.vcf; fi
  vep \
    -i INPUT_FOR_VEP.vcf \
    -o ${sample}.vep.vcf \
    --offline --cache --dir_cache ${vep_cache} \
    --dir_plugins ${vep_plugins} \
    --fasta ${vep_fasta} \
    --assembly GRCh38 --species homo_sapiens \
    --hgvs --symbol --vcf --everything --canonical --merged \
    --plugin REVEL,${revel_vcf} \
    --plugin AlphaMissense,file=${alpha_missense_vcf},cols=am_pathogenicity:am_class \
    --plugin SpliceAI,snv=${spliceai_snv_vcf},indel=${spliceai_indel_vcf} \
    --plugin BayesDel,file=${bayesdel_vcf} \
    --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,ALLELEID
  """
}

process LeanReport {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay   
  publishDir "${params.outdir}/${meta.sample}/reports", mode: 'copy'
  input:
    tuple val(meta),
          path(vcf), path(exon_cov), path(r1r2), path(frstrand),
          path(flagstat), path(stats),
          path(mosdepth_summary),
          path(sex_check),
          path(gaps20), path(gaps30),
          path(thresholds),
          path(somatic_vcf)
    each path(script)
  output:
    tuple val(meta), path("${meta.sample}_report/${meta.sample}_variants_lean.xlsx")
  script:
    def sample  = meta.sample
    def som_arg = somatic_vcf.name != 'NO_FILE' ? "--somatic-vcf ${somatic_vcf}" : ""
    // Which BAM produced the germline/ACMG calls. Set per-sample in main.nf for
    // somatic runs; a plain germline run has one sample and no ambiguity.
    def gsrc    = meta.germline_source ?: 'germline (single sample)'
  """
  mkdir -p ${sample}_report
  python ${script} \
    $vcf $exon_cov $r1r2 $frstrand ${sample}_report/${sample}_variants_lean.xlsx \
    --sample-id ${sample} --assay ${meta.assay} --build GRCh38 \
    --flagstat ${flagstat} --stats ${stats} \
    --mosdepth-summary ${mosdepth_summary} \
    --acmg-thresholds ${thresholds} \
    --sexcheck ${sex_check} \
    --gaps20 ${gaps20} --gaps30 ${gaps30} \
    --germline-source '${gsrc}' \
    ${som_arg}
  """
}

// ── Somatic-specific VCF processing (no publishDir on intermediates) ──────────

// Restrict the somatic call set to the reporting panel. Runs before
// NormalizeSomatic so the expensive VEP step only annotates in-panel variants.
process BedFilterSomatic {
  tag { "${meta.sample}" }
  input:
    tuple val(meta), path(vcf)
    path bed
  output:
    tuple val(meta), path("${meta.sample}.somatic.panel.vcf.gz")
  script:
    def sample = meta.sample
  """
  tabix -p vcf $vcf || bcftools index -t $vcf
  bcftools view -R $bed $vcf -Oz -o ${sample}.somatic.panel.vcf.gz
  tabix -p vcf ${sample}.somatic.panel.vcf.gz
  """
}

process NormalizeSomatic {
  tag { "${meta.sample}" }
  input:
    tuple val(meta), path(vcf)
    path fasta
    path fai
  output:
    tuple val(meta), path("${meta.sample}.somatic.norm.vcf.gz")
  script:
    def sample = meta.sample
    // Normalization only — no PASS gate and no post-hoc thresholds here.
    // VEP runs once over the full in-panel set, and the report script applies
    // PASS + the tumor-only thresholds when building each sheet. Keeping the
    // filtering in one place (Python) also lets the QC tab explain *why* a
    // variant was excluded instead of it silently disappearing upstream.
  """
  # -f left-aligns indels against the alignment reference. Mutect2 emits a high
  # proportion of indels and this is the path where indel representation drives
  # interpretation -- a non-parsimonious indel matches neither ClinVar nor a
  # position-keyed hotspot list. See NormalizeVCF for the full rationale.
  bcftools norm -m -any -f $fasta $vcf -Oz -o ${sample}.somatic.norm.vcf.gz
  tabix -p vcf ${sample}.somatic.norm.vcf.gz
  """
}


process VEP_Annotate_Somatic {
  tag { "${meta.sample}" }
  publishDir "${params.outdir}/${meta.sample}/somatic", mode: 'copy'
  input:
    tuple val(meta), path(vcf)
    path vep_cache
    path vep_fasta
    path vep_fasta_fai
    path revel_vcf
    path revel_vcf_tbi
    path alpha_missense_vcf
    path alpha_missense_vcf_tbi
    path clinvar_vcf
    path clinvar_vcf_tbi
    path spliceai_snv_vcf
    path spliceai_snv_vcf_tbi
    path spliceai_indel_vcf
    path spliceai_indel_vcf_tbi
    path bayesdel_vcf
    path bayesdel_vcf_tbi
    path vep_plugins
  output:
    tuple val(meta), path("${meta.sample}.somatic.vep.vcf")
  script:
    def sample = meta.sample
  """
  set -euo pipefail
  if [[ "$vcf" == *.vcf.gz ]]; then gunzip -c "$vcf" > INPUT_FOR_VEP.vcf; else cp "$vcf" INPUT_FOR_VEP.vcf; fi
  vep \
    -i INPUT_FOR_VEP.vcf \
    -o ${sample}.somatic.vep.vcf \
    --offline --cache --dir_cache ${vep_cache} \
    --dir_plugins ${vep_plugins} \
    --fasta ${vep_fasta} \
    --assembly GRCh38 --species homo_sapiens \
    --hgvs --symbol --vcf --everything --canonical --merged \
    --plugin REVEL,${revel_vcf} \
    --plugin AlphaMissense,file=${alpha_missense_vcf},cols=am_pathogenicity:am_class \
    --plugin SpliceAI,snv=${spliceai_snv_vcf},indel=${spliceai_indel_vcf} \
    --plugin BayesDel,file=${bayesdel_vcf} \
    --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,ALLELEID
  """
}

process GENERATE_ACMG_REPORT {
  tag { "${meta.sample} (${meta.assay})" }
  publishDir "${params.outdir}/${meta.sample}/reports", mode: 'copy'
  input:
    tuple val(meta), path(excel_file)
    each path(python_skeleton)
    each path(template_dir)
  output:
    tuple val(meta), path("${meta.sample}_report/${meta.sample}_clinical_report.html")
  script:
    def sample = meta.sample
    def assay = meta.assay
  """
  python ${python_skeleton}/generate_report.py \
    ${excel_file} ${sample}_report \
    --sample-id ${sample} \
    --assay ${assay} \
    --template-dir ${template_dir} \
    --format html
  """
}


/******************  SUBWORKFLOW: consumes Sarek outputs  ********************/

workflow POST_SAREK {
  take:
    vcf_ch          // tuple(meta, vcf)
    bam_ch          // tuple(meta, bam, bai)
    bed_ch          // value channel with BED
    somatic_vcf_ch  // tuple(meta, somatic_vep_vcf) — pass Channel.empty() if not somatic mode

  main:
    // join per-sample → (s//ample, vcf, bam, bai)
    sample_inputs = vcf_ch.join(bam_ch)
    script_ch = Channel.fromPath("${params.scriptdir}/generate_lean_report_org.py").first()
    report_script_ch = Channel.fromPath("${params.scriptdir}/python-skeleton/", type: 'dir').first()
    template_dir_ch = Channel.fromPath("${params.template_dir}", type: 'dir').first()
    // VCF path
    BedFilterVCF(sample_inputs.map { s, vcf, bam, bai -> tuple(s, vcf) }, bed_ch)
    def _normFasta = resolveNormFasta()
    norm_fasta_ch = Channel.value(file(_normFasta))
    norm_fai_ch   = Channel.value(file("${_normFasta}.fai"))

    NormalizeVCF(BedFilterVCF.out, norm_fasta_ch, norm_fai_ch)
    FilterVCF(NormalizeVCF.out)
    // No FORMAT/VAF tag is written. `bcftools +fill-tags` used to run here, but
    // NormalizeVCF already split multiallelics, so fill-tags only ever saw
    // biallelic records with AD=[site_ref, this_alt] and computed
    // alt/(ref+alt). At a 1/2 site the sibling allele's reads are not in the
    // record, so there is no denominator to compute against. Reordering cannot
    // fix it: the split has to happen before VEP.
    //
    // The report recomputes VAF from AD against site-level depth
    // (AD_ref + sum of every ALT's AD) and overwrites every row.
    vep_ch = params.run_vep ? VEP_Annotate(
      FilterVCF.out, 
      file(params.vep_cache), 
      file(params.vep_fasta), 
      file(params.vep_fasta + ".fai"), 
      file(params.revel_vcf), 
      file(params.revel_vcf + ".tbi"), 
      file(params.alpha_missense_vcf), 
      file(params.alpha_missense_vcf + ".tbi"), 
      file(params.clinvar_vcf), 
      file(params.clinvar_vcf + ".tbi"), 
      file(params.spliceai_snv_vcf), 
      file(params.spliceai_snv_vcf + ".tbi"),
      file(params.spliceai_indel_vcf), 
      file(params.spliceai_indel_vcf + ".tbi"),
      file(params.bayesdel_vcf), 
      file(params.bayesdel_vcf + ".tbi"),
      file(params.vep_plugins)
      ) : FilterVCF.out  // (sample, vcf)

    // BAM path
    BedFilterBAM(sample_inputs.map { s, vcf, bam, bai -> tuple(s, vcf, bam) }, bed_ch)
    bam_sample_ch = BedFilterBAM.out.map { s, bam, bai -> tuple(s, bam, bai) }

    CoverageSummary(bam_sample_ch.map { s, bam, bai -> tuple(s, bam) }, bed_ch)
    R1R2Ratio(bam_sample_ch, bed_ch)
    ForwardReverseRatio(bam_sample_ch, bed_ch)
    SamtoolsFlagstat(bam_sample_ch.map { s, bam, bai -> tuple(s, bam) })
    SamtoolsStats(bam_sample_ch.map { s, bam, bai -> tuple(s, bam) })
    MosdepthRun(bam_sample_ch, bed_ch)
    CoverageGapsAnnotation(MosdepthRun.out.map { s, summary, thresholds, quantized -> tuple(s, quantized, thresholds) }, bed_ch)
    SexCheck(bam_sample_ch.map { s, bam, bai -> tuple(s, bam) })
    BcftoolsStats(vep_ch.map { s, vcf -> tuple(s, vcf) })

    // prepare joins keyed by sample
    exon_cov_ch         = CoverageSummary.out.map { s, summary -> tuple(s, summary) }
    gaps20_ch           = CoverageGapsAnnotation.out.map { s, a20, a30 -> tuple(s, a20) }
    gaps30_ch           = CoverageGapsAnnotation.out.map { s, a20, a30 -> tuple(s, a30) }
    mosdepth_summary_ch = MosdepthRun.out.map { s, summary, thresholds, quantized -> tuple(s, summary) }
    thresholds_ch       = MosdepthRun.out.map { s, summary, thresholds, quantized -> tuple(s, thresholds) }

    // join all for LeanReport. Every channel above is keyed by the SAME meta map
    // (they all descend from sample_inputs), so joining on meta is safe here.
    qc_joined_ch = vep_ch
      .join(exon_cov_ch)
      .join(R1R2Ratio.out)
      .join(ForwardReverseRatio.out)
      .join(SamtoolsFlagstat.out)
      .join(SamtoolsStats.out)
      .join(mosdepth_summary_ch)
      .join(SexCheck.out)
      .join(gaps20_ch)
      .join(gaps30_ch)
      .join(thresholds_ch)

    // somatic_vcf_ch is OPTIONAL and comes from a different subworkflow, so its
    // meta is built independently -- it carries no germline_source, while the
    // chain above does. Joining those two on the meta MAP matched nothing and
    // `remainder: true` then emitted the unmatched somatic entry as a 3-element
    // tuple [meta, null, vcf] into a closure expecting 13, which is a
    // MissingMethodException at runtime, not a parse error.
    //
    // Key this join on the sample NAME instead. It is the one field both sides
    // agree on by construction, and it stays correct if either meta gains a
    // field later.
    som_keyed_ch = somatic_vcf_ch.map { meta, somatic_vcf -> tuple(meta.sample, somatic_vcf) }

    lean_input_ch = qc_joined_ch
      .map { items -> tuple(items[0].sample, items) }
      .join(som_keyed_ch, remainder: true)
      // A somatic VCF with no matching QC chain would arrive with items == null.
      // That should be impossible (both derive from the same sample set) but
      // dropping it is cheaper than a confusing failure 15 hours in.
      .filter { sample, items, somatic_vcf -> items != null }
      .map { sample, items, somatic_vcf ->
          items + [ somatic_vcf != null ? somatic_vcf : file("${projectDir}/assets/NO_FILE") ]
      }
    LeanReport(lean_input_ch, script_ch)
    GENERATE_ACMG_REPORT(LeanReport.out, report_script_ch, template_dir_ch)
}

// ── Somatic VEP annotation subworkflow ────────────────────────────────────────

workflow POST_SAREK_SOMATIC {
  take:
    vcf_ch          // tuple(meta, rescued_vcf) — output of MUTECT2_RESCUE (with meta)
    somatic_bed_ch  // value channel with the somatic reporting panel BED

  main:
    // Restrict to the reporting panel first — everything downstream, VEP
    // included, only ever sees in-panel variants.
    BedFilterSomatic(vcf_ch, somatic_bed_ch)

    // Annotate the FULL in-panel set once. PASS and the tumor-only thresholds
    // are applied downstream in the report script, so VEP runs a single time
    // per sample regardless of how many views of the data the report needs.
    def _somNormFasta = resolveNormFasta()
    som_fasta_ch = Channel.value(file(_somNormFasta))
    som_fai_ch   = Channel.value(file("${_somNormFasta}.fai"))

    NormalizeSomatic(BedFilterSomatic.out, som_fasta_ch, som_fai_ch)
    // No FORMAT/VAF tag. `bcftools +fill-tags` used to run here and nothing read
    // its output: the report uses Mutect2's own FORMAT/AF, a model posterior
    // that accounts for read orientation and filtered depth. The AD ratio
    // fill-tags computes reads 1.000 whenever there are zero ref reads -- on as
    // few as 6 reads -- where AF reports 0.833. AF is also multiallelic-aware,
    // which the post-split tag was not.
    somatic_vep_ch = params.run_vep ? VEP_Annotate_Somatic(
      NormalizeSomatic.out,
      file(params.vep_cache),
      file(params.vep_fasta),
      file(params.vep_fasta + ".fai"),
      file(params.revel_vcf),
      file(params.revel_vcf + ".tbi"),
      file(params.alpha_missense_vcf),
      file(params.alpha_missense_vcf + ".tbi"),
      file(params.clinvar_vcf),
      file(params.clinvar_vcf + ".tbi"),
      file(params.spliceai_snv_vcf),
      file(params.spliceai_snv_vcf + ".tbi"),
      file(params.spliceai_indel_vcf),
      file(params.spliceai_indel_vcf + ".tbi"),
      file(params.bayesdel_vcf),
      file(params.bayesdel_vcf + ".tbi"),
      file(params.vep_plugins)
    ) : NormalizeSomatic.out

  emit:
    // tuple(meta, <sample>.somatic.vep.vcf) — every in-panel call, annotated,
    // FILTER preserved. The report script derives the reported and QC views.
    somatic_vep = somatic_vep_ch
}
