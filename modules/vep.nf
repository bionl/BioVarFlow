// modules/vep.nf
nextflow.enable.dsl=2

// -------- Parameters (used by processes) --------
params.bed         = params.bed         ?: "${workflow.projectDir}/data/ACMG_SF_plus_HemOnc_MANE_exons_50bp.refseq.bed"
params.sf_genes    = params.sf_genes    ?: "${workflow.projectDir}/data/acmg_sf_gene_list.txt"
params.hemonc_genes= params.hemonc_genes?: "${workflow.projectDir}/data/hemonc.txt"
params.outdir      = params.outdir      ?: "results"
params.scriptdir   = params.scriptdir   ?: "${workflow.projectDir}/scripts"
params.template_dir= params.template_dir?: "${workflow.projectDir}/scripts/template-files"

params.run_vep     = params.run_vep     ?: true
params.min_dp   = params.min_dp   ?: 10
params.min_qual = params.min_qual ?: 10
// Rows kept on the Exomiser tab, by rank. 0 = keep everything Exomiser emitted
// (745 for IQMM); set a positive number to cap it.
params.exomiser_top = params.exomiser_top ?: 0
// Cross-case table: 'all' keeps every row of the source sheet, 'coding'
// restricts to protein-altering and splice-affecting consequences.
params.merge_filter = params.merge_filter ?: 'all'

// Reference used by NormalizeVCF for left-alignment/trimming: see the lazy
// fallback to params.fasta inside POST_SAREK. Deliberately NOT resolved here at
// module script level — params.fasta is assigned in external/sarek/main.nf
// (line 41, getGenomeAttribute), so its value depends on include order and
// would be null if this module were ever included first.
// Override with --norm_fasta to use a local/GCS mirror and avoid igenomes egress.
params.norm_fasta  = params.norm_fasta  ?: null
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
  #          breaks position-keyed annotation lookups — SpliceAI silently returns
  #          nothing for the affected records.
  #
  # The reference MUST be the one the BAMs were aligned to (GATK
  # Homo_sapiens_assembly38.fasta, chr-prefixed contigs). Passing the Ensembl
  # VEP fasta here fails on every record: its contigs are named 1/2/…/MT.
  #
  # No -c override: bcftools exits on a REF mismatch by default, which is what we
  # want — a mismatch means the wrong reference, and continuing would silently
  # corrupt allele representations.
  bcftools norm -m -any -f $fasta $vcf -Oz -o ${sample}.normalized.vcf.gz
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

  # Label each gap by the EXON it falls in (chrom:start-end from the -b BED),
  # not by the BED's 4th column. That column is GENE|TRANSCRIPT and repeats
  # across every exon of a gene, so using it collapsed all gaps to gene level and
  # stamped one gene-wide total onto every exon row of the report.
  #
  # Also CLIP each gap to the exon: bedtools merge upstream can join adjacent
  # low-coverage runs that span introns, so an unclipped gap routinely reported
  # more bp than the exon it was attached to.
  #   -a cols: \$1,\$2,\$3 = gap    -b cols: \$4,\$5,\$6 = exon, \$7 = GENE|TRANSCRIPT
  # -wo (not -wao) emits only genuinely overlapping pairs.
  bedtools intersect -wo -a ${sample}.acmg_gaps_lt20.bed -b $bed \
    | awk 'BEGIN{OFS="\\t"}{ s=(\$2>\$5?\$2:\$5); e=(\$3<\$6?\$3:\$6); if(e>s) print \$1,s,e,\$4":"\$5"-"\$6 }' \
    > ${sample}.acmg_gaps_lt20.annot.bed

  bedtools intersect -wo -a ${sample}.acmg_gaps_lt30.bed -b $bed \
    | awk 'BEGIN{OFS="\\t"}{ s=(\$2>\$5?\$2:\$5); e=(\$3<\$6?\$3:\$6); if(e>s) print \$1,s,e,\$4":"\$5"-"\$6 }' \
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

// ── Per-variant, per-allele strand bias from the RAW alignment ──────────────
//
// The existing strand QC (ForwardReverseRatio) is per-EXON and counts every
// read regardless of allele, so it cannot see a skew confined to the alt reads.
// GATK's INFO/FS can, but it is computed on post-reassembly allele depths: for
// MSH2 chr2:47414420 T>G it reported FS=43.4 (under the standard >60 filter)
// while the raw pileup showed 0 forward / 11 reverse alt reads against 49/8
// reference -- Fisher p = 4.9e-08, i.e. FS ~73.
//
// bcftools mpileup recomputes allele depths from the alignment itself. ADF/ADR
// are the per-allele forward/reverse counts, giving the 2x2 table directly, and
// unlike a plain samtools pileup they cover indels as well as SNVs. --no-BAQ
// keeps indel-adjacent bases from being down-weighted.
process StrandBiasPileup {
  tag { "${meta.sample} (${meta.assay})" }
  input:
    // raw_bam is the full markduplicates BAM, NOT the BED-filtered one. The
    // panel BAM cannot answer anything about off-panel variants -- those reads
    // were removed by `samtools view -L` -- and that left every EXOMISER_ONLY
    // row with a blank StrandBias, which reads as "clean" to a reviewer.
    // At a panel position the two BAMs give identical counts anyway: a read
    // overlapping a panel base overlaps the BED by definition, so it survives
    // the filter. Using the raw BAM for both therefore costs nothing in
    // accuracy and removes the blind spot.
    tuple val(meta), path(vep_vcf), path(cons_vcf), path(cons_tbi),
          path(exomiser_tsv), path(raw_bam), path(raw_bai)
    path fasta
    path fai
  output:
    tuple val(meta), path("${meta.sample}_mpileup_adf_adr.tsv"),
                     path("${meta.sample}_sb_sites.vcf.gz")
  script:
    def sample = meta.sample
  """
  set -euo pipefail

  # Sites to test = panel variants (the VEP VCF) + everything Exomiser ranked.
  # Exomiser writes CONTIG without the 'chr' prefix, so add it back.
  bcftools query -f '%CHROM\\t%POS\\n' $vep_vcf > sites_raw.txt
  if [ -s "$exomiser_tsv" ]; then
    awk -F'\\t' 'NR>1 && \$15 != "" {
        c = \$15; if (c !~ /^chr/) c = "chr" c; print c "\\t" \$16
      }' $exomiser_tsv >> sites_raw.txt
  fi
  sort -u sites_raw.txt > sites.txt

  # Subset the consensus VCF to those positions with awk rather than
  # `bcftools view -R/-T`: both expect the positions file in the VCF's own
  # contig order, and ours is a lexicographic sort -u (chr10 before chr2). A
  # hash join sidesteps the ordering question entirely, and the output keeps
  # the consensus VCF's order -- which is what the regions file below needs.
  #
  # The consensus VCF is used, not the VEP VCF, because only it contains the
  # off-panel records, with the REF/ALT and FORMAT/AD that strand_bias.py needs
  # for allele matching and the indel VAF-reconciliation gate.
  bcftools view $cons_vcf \\
  | awk -F'\\t' 'NR==FNR { keep[\$1"\\t"\$2]; next }
                 /^#/ { print; next }
                 (\$1"\\t"\$2) in keep' sites.txt - \\
  | bgzip -c > ${sample}_sb_sites.vcf.gz

  # Regions for mpileup, derived from the subset VCF so they are already in
  # reference order. -R seeks via the BAM index instead of streaming the whole
  # file, which matters now that this runs on the full exome BAM.
  bcftools query -f '%CHROM\\t%POS\\n' ${sample}_sb_sites.vcf.gz | uniq > regions.txt

  # Deliberately NOT `bcftools call -C alleles`. Constraining the pileup to the
  # caller's alleles looks like the right way to fix indel counting, and it does
  # match more indels -- but `call` makes a GENOTYPE decision, and when alt
  # support is weak and one-sided it calls the site hom-ref and drops the ALT
  # entirely. That is exactly the artifact class this check exists to find: on
  # IQMM it returned ALT='.' for MSH2, RUNX1 and DSG2 alike and flagged nothing.
  # Plain mpileup reports what it sees and lets the Fisher test decide.
  bcftools mpileup \\
    --regions-file regions.txt \\
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR \\
    --fasta-ref $fasta \\
    --min-BQ 13 --min-MQ 0 --no-BAQ --max-depth 8000 \\
    -Ou $raw_bam \\
  | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%ADF\\t%ADR]\\n' \\
  > ${sample}_mpileup_adf_adr.tsv
  """
}

// Fisher test over the pileup counts. Split from StrandBiasPileup only because
// of containers: the bcftools biocontainer has no python3, and this half runs in
// the same python-bionl image as the other report scripts.
process StrandBiasTest {
  tag { "${meta.sample} (${meta.assay})" }
  publishDir "${params.outdir}/${meta.sample}/qc", mode: 'copy'
  input:
    tuple val(meta), path(pileup), path(vcf)
    each path(script)
  output:
    tuple val(meta), path("${meta.sample}_strand_bias.tsv")
  script:
    def sample = meta.sample
  """
  python3 ${script} \\
    --mpileup ${pileup} \\
    --vcf $vcf \\
    --out ${sample}_strand_bias.tsv
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
    --pick --pick_order mane_select,mane_plus_clinical,canonical,tsl,biotype,ccds,rank,length \
    --plugin REVEL,${revel_vcf} \
    --plugin AlphaMissense,file=${alpha_missense_vcf},cols=am_pathogenicity:am_class \
    --plugin SpliceAI,snv=${spliceai_snv_vcf},indel=${spliceai_indel_vcf} \
    --plugin BayesDel,file=${bayesdel_vcf} \
    --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,ALLELEID
  """
}

// ── Cross-case table ────────────────────────────────────────────────────────
//
// One row per (case, variant) across every sample in the run, in the agreed
// deliverable format. Runs once, after every per-sample workbook exists --
// hence the .collect() on LeanReport.out at the call site.
//
// Sourcing is per case: the Prioritised sheet where Exomiser ran, the panel
// Reportable sheet where it could not (no HPO terms), recorded in a Ranked_By
// column. See the script's module docstring for why there is no middle option.
process MergeVariantsTable {
  tag { "${xlsx instanceof List ? xlsx.size() : 1} case(s)" }
  publishDir "${params.outdir}/merged", mode: 'copy'
  input:
    path xlsx
    each path(script)
  output:
    path "genetic_variants_all_samples.xlsx", emit: table
  script:
  """
  python ${script} ${xlsx} \\
    --out genetic_variants_all_samples.xlsx \\
    --filter ${params.merge_filter}
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
          path(strand_bias),
          path(exomiser_tsv)
    each path(script)
    each path(sf_genes_file)
    each path(hemonc_genes_file)
  output:
    tuple val(meta), path("${meta.sample}_report/${meta.sample}_variants.xlsx")
  script:
    def sample = meta.sample
  """
  mkdir -p ${sample}_report
  python ${script} \
    $vcf $exon_cov $r1r2 $frstrand ${sample}_report/${sample}_variants.xlsx \
    --sample-id ${sample} --assay ${meta.assay} --build GRCh38 \
    --flagstat ${flagstat} --stats ${stats} \
    --mosdepth-summary ${mosdepth_summary} \
    --acmg-thresholds ${thresholds} \
    --sexcheck ${sex_check} \
    --sf-genes ${sf_genes_file} \
    --hemonc-genes ${hemonc_genes_file} \
    --gaps20 ${gaps20} --gaps30 ${gaps30} \
    --strand-bias ${strand_bias} \
    --exomiser ${exomiser_tsv} \
    --exomiser-top ${params.exomiser_top}
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
    exomiser_ch // (sample, <sample>_exomiser.variants.tsv) -- may be empty
    raw_cons_ch // (sample, consensus vcf, tbi) -- unfiltered; may be empty
    vcf_ch   // (sample, vcf)
    bam_ch//  // // (samp//le, bam, bai)
    bed_ch   // value channel with //BED

  main:
    // join per-sample → (s//ample, vcf, bam, bai)
    sample_inputs = vcf_ch.join(bam_ch)
    script_ch = Channel.fromPath("${params.scriptdir}/generate_lean_report_org.py").first()
    report_script_ch = Channel.fromPath("${params.scriptdir}/python-skeleton/", type: 'dir').first()
    template_dir_ch = Channel.fromPath("${params.template_dir}", type: 'dir').first()
    // Gene-set lists for the lean report. Staged as value channels so Nextflow
    // stages them into every LeanReport task (works both locally and on cloud).
    // Overridable via --sf_genes / --hemonc_genes (defaults in params block above).
    sf_genes_ch     = Channel.fromPath(params.sf_genes).first()
    hemonc_genes_ch = Channel.fromPath(params.hemonc_genes).first()
    // VCF path
    BedFilterVCF(sample_inputs.map { s, vcf, bam, bai -> tuple(s, vcf) }, bed_ch)
    // Reference for left-alignment/trimming. params.fasta is the GATK
    // GRCh38 assembly Sarek aligned against — NOT params.vep_fasta (Ensembl,
    // differently-named contigs). Overridable via --norm_fasta for sites that
    // keep a local mirror instead of pulling from igenomes.
    // Resolved here, not at module load: params.fasta is set by
    // external/sarek/main.nf when it is included, so reading it inside the
    // workflow body guarantees it is populated regardless of include order.
    // Same fallback chain as main.nf's resolveRefFasta(): params.fasta is set by
    // external/sarek/main.nf and is not reliably visible from another module's
    // binding, so fall through to the igenomes map, which config-parse time
    // always populates.
    def _normFasta = params.norm_fasta ?: params.fasta
    if (!_normFasta && params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        _normFasta = params.genomes[params.genome].fasta
    }
    if (!_normFasta) {
        error "❌ No reference for NormalizeVCF.\n" +
              "   Tried --norm_fasta, params.fasta, params.genomes[${params.genome}].fasta.\n" +
              "   genome=${params.genome}  genomes_loaded=${params.genomes ? params.genomes.size() : 0}\n" +
              "   Set --norm_fasta to the GATK assembly the BAMs were aligned against."
    }
    norm_fasta_ch = Channel.value(file(_normFasta))
    norm_fai_ch   = Channel.value(file("${_normFasta}.fai"))

    // No FORMAT/VAF tag is written. `bcftools +fill-tags` used to run here, but
    // CONSENSUS_CALLING already splits multiallelics (NormalizeDV/NormalizeHC
    // run `norm -m -any`), so fill-tags only ever saw biallelic records with
    // AD=[site_ref, this_alt] and computed alt/(ref+alt). At a 1/2 site the
    // sibling allele's reads are not in the record, so there is no denominator
    // to compute against: ACTC1 chr15:34791307 came out 1.000 on BOTH alleles
    // instead of 0.585/0.415, and MSH2 chr2:47414420 T>G read 0.79 instead of
    // 0.34. Reordering cannot fix it -- the split happens upstream.
    //
    // Verified on IQMM: of 915 records the tag was right on the 907 biallelic
    // rows and wrong on all 8 multiallelic ones, and nothing read it. The
    // report recomputes VAF from AD against site-level depth
    // (AD_ref + sum of every ALT's AD) and overwrites all 915 rows, so a wrong
    // value in a published clinical VCF was the tag's only remaining effect.
    NormalizeVCF(BedFilterVCF.out, norm_fasta_ch, norm_fai_ch)
    FilterVCF(NormalizeVCF.out)
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

    // Per-variant strand bias, recomputed from the raw alignment (see the
    // StrandBiasPileup header for why the caller's own FS is not enough).
    strand_bias_script_ch = Channel.fromPath("${params.scriptdir}/strand_bias.py").first()

    // Sites come from the panel VCF AND from Exomiser; alleles and AD come from
    // the unfiltered consensus, which is the only VCF holding the off-panel
    // records. Reads come from the RAW BAM (bam_ch), not BedFilterBAM.out --
    // see the StrandBiasPileup header.
    //
    // remainder:true on both joins so a sample with no Exomiser result, or a
    // run with no consensus at all, still gets its panel variants tested.
    no_file_sb = file("${workflow.projectDir}/assets/NO_FILE")
    sb_input_ch = vep_ch
      .join(raw_cons_ch,  remainder: true)
      .join(exomiser_ch,  remainder: true)
      .join(bam_ch,       remainder: true)
      .filter { it[1] != null && it[2] != null && it[5] != null }
      .map { s, vep, cons, tbi, exo, bam, bai ->
             tuple(s, vep, cons, tbi, exo ?: no_file_sb, bam, bai) }

    StrandBiasPileup(sb_input_ch, norm_fasta_ch, norm_fai_ch)
    StrandBiasTest(StrandBiasPileup.out, strand_bias_script_ch)

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
      .join(StrandBiasTest.out)

    // Exomiser runs as ONE batch for the whole run, so its per-sample TSVs are
    // joined back here by sample name. remainder:true plus the NO_FILE stand-in
    // keeps LeanReport running for samples Exomiser skipped, and -- with
    // errorStrategy 'ignore' on EXOMISER_BATCH -- keeps a failed Exomiser from
    // costing the reports entirely. The workbook simply omits the two tabs.
    no_exomiser = file("${workflow.projectDir}/assets/NO_FILE")
    lean_with_exo_ch = lean_input_ch
      .map { tup -> tuple(tup[0].sample, tup) }
      .join(exomiser_ch, remainder: true)
      .filter { s, tup, exo -> tup != null }
      .map    { s, tup, exo -> tup + [ exo ?: no_exomiser ] }
    LeanReport(lean_with_exo_ch, script_ch, sf_genes_ch, hemonc_genes_ch)
    GENERATE_ACMG_REPORT(LeanReport.out, report_script_ch, template_dir_ch)

    // collect() so this fires once, with every sample's workbook staged.
    merge_script_ch = Channel.fromPath("${params.scriptdir}/merge_reportable_variants.py").first()
    MergeVariantsTable(LeanReport.out.map { meta, xlsx -> xlsx }.collect(), merge_script_ch)
}
