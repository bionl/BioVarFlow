// modules/vep.nf
nextflow.enable.dsl=2

// -------- Parameters (used by processes) --------
params.bed         = params.bed         ?: "${workflow.projectDir}/data/annotated_merged_MANE_deduped.bed"
params.outdir      = params.outdir      ?: "results"
params.scriptdir   = params.scriptdir   ?: "${workflow.projectDir}/scripts"
params.template_dir= params.template_dir?: "${workflow.projectDir}/scripts/template-files"

params.run_vep     = params.run_vep     ?: true
params.min_dp   = params.min_dp   ?: 10
params.min_qual = params.min_qual ?: 10
// VEP resource params expected from main/config:
// params.vep_fasta, params.revel_vcf, params.alpha_missense_vcf, params.clinvar_vcf


/********************  PROCESSES (unchanged logic, publish to per-sample dirs)  ********************/

process BedFilterVCF {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  publishDir "${params.outdir}/${meta.sample}/vcf", mode: 'copy'
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
  publishDir "${params.outdir}/${meta.sample}/vcf", mode: 'copy'
  input:
    tuple val(meta), path(vcf)
  output:
    tuple val(meta), path("${meta.sample}.normalized.vcf.gz")
  script:
    def sample = meta.sample
  """
  bcftools norm -m -any $vcf -Oz -o ${sample}.normalized.vcf.gz
  tabix -p vcf ${sample}.normalized.vcf.gz
  """
}

process FilterVCF {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  publishDir "${params.outdir}/${meta.sample}/vcf", mode: 'copy'
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

process AddVAF {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay 
  publishDir "${params.outdir}/${meta.sample}/vcf", mode: 'copy'
  input:
    tuple val(meta), path(vcf)
  output:
    tuple val(meta), path("${meta.sample}.vaf_added.vcf.gz")
  script:
    def sample = meta.sample
  """
  bcftools +fill-tags $vcf -Oz -o ${sample}.vaf_added.vcf.gz -- -t FORMAT/VAF
  tabix -p vcf ${sample}.vaf_added.vcf.gz
  """
}

process BedFilterBAM {
  tag { "${meta.sample} (${meta.assay})" } // meta is a map containing sample and assay
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
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
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
  input:
    tuple val(meta), path(bam)
    path bed
  output:
    tuple val(meta), path("${meta.sample}_coverage_summary.sorted.txt"), path("${meta.sample}_coverage_per_base.txt")
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
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
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
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
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
      path("${meta.sample}.acmg_gaps_lt20.bed"),
      path("${meta.sample}.acmg_gaps_lt30.bed"),
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
          path(thresholds)
    each path(script)
  output:
    tuple val(meta), path("${meta.sample}_report/${meta.sample}_variants_lean.xlsx")
  script:
    def sample = meta.sample
  """
  mkdir -p ${sample}_report
  python ${script} \
    $vcf $exon_cov $r1r2 $frstrand ${sample}_report/${sample}_variants_lean.xlsx \
    --sample-id ${sample} --assay ${meta.assay} --build GRCh38 \
    --flagstat ${flagstat} --stats ${stats} \
    --mosdepth-summary ${mosdepth_summary} \
    --acmg-thresholds ${thresholds} \
    --sexcheck ${sex_check} \
    --gaps20 ${gaps20} --gaps30 ${gaps30}
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
    vcf_ch   // (sample, vcf)
    bam_ch//  // // (samp//le, bam, bai)
    bed_ch   // value channel with //BED

  main:
    // join per-sample → (s//ample, vcf, bam, bai)
    sample_inputs = vcf_ch.join(bam_ch)
    script_ch = Channel.fromPath("${params.scriptdir}/generate_lean_report_org.py").first()
    report_script_ch = Channel.fromPath("${params.scriptdir}/python-skeleton/", type: 'dir').first()
    template_dir_ch = Channel.fromPath("${params.template_dir}", type: 'dir').first()
    // VCF path
    BedFilterVCF(sample_inputs.map { s, vcf, bam, bai -> tuple(s, vcf) }, bed_ch)
    NormalizeVCF(BedFilterVCF.out)
    FilterVCF(NormalizeVCF.out)
    AddVAF(FilterVCF.out)
    vep_ch = params.run_vep ? VEP_Annotate(
      AddVAF.out, 
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
      ) : AddVAF.out  // (sample, vcf)

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
    exon_cov_ch         = CoverageSummary.out.map { s, summary, per_base -> tuple(s, summary) }
    gaps20_ch           = CoverageGapsAnnotation.out.map { s, g20, g30, a20, a30 -> tuple(s, a20) }
    gaps30_ch           = CoverageGapsAnnotation.out.map { s, g20, g30, a20, a30 -> tuple(s, a30) }
    mosdepth_summary_ch = MosdepthRun.out.map { s, summary, thresholds, quantized -> tuple(s, summary) }
    thresholds_ch       = MosdepthRun.out.map { s, summary, thresholds, quantized -> tuple(s, thresholds) }

    // join all for LeanReport
    lean_input_ch = vep_ch
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
    LeanReport(lean_input_ch, script_ch)
    GENERATE_ACMG_REPORT(LeanReport.out, report_script_ch, template_dir_ch)
}
