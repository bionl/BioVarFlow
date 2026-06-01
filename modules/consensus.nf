// modules/consensus.nf
nextflow.enable.dsl=2

params.outdir = params.outdir ?: "results"
params.ref_fasta = params.ref_fasta ?: params.vep_fasta

/********************  CONSENSUS CALLING PROCESSES  ********************/
process CONS_REHEADER_VCF {
  tag "${sample}:${caller}"

  input:
    tuple val(sample), val(caller), path(vcf)

  output:
    tuple val(sample), val(caller), path("${sample}.${caller}.rehead.vcf.gz"), path("${sample}.${caller}.rehead.vcf.gz.tbi")

  script:
  """
  set -euo pipefail
  
  n=\$(bcftools query -l $vcf | wc -l)
  if [ "\$n" -ne 1 ]; then
    echo "ERROR: expected single-sample VCF, found \$n samples in $vcf" >&2
    exit 2
  fi

  old=\$(bcftools query -l $vcf | head -n 1)

  if [ "\$old" = "${sample}" ]; then
    cp $vcf ${sample}.${caller}.rehead.vcf.gz
  else
    echo -e "\${old}\\t${sample}" > rename.txt
    bcftools reheader -s rename.txt -o ${sample}.${caller}.rehead.vcf.gz $vcf
  fi

  bcftools index -t ${sample}.${caller}.rehead.vcf.gz
  """
}

process NormalizeDV {
  tag "$sample"
  //publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(vcf)
    path ref_fasta
    path ref_fai
  
  output:
    tuple val(sample), path("${sample}.dv.norm.pass.vcf.gz"), path("${sample}.dv.norm.pass.vcf.gz.tbi")
  
  script:
  """
  #bcftools norm -f ${ref_fasta} ${vcf} -Oz -o ${sample}.dv.norm.vcf.gz
  bcftools norm -m -any ${vcf} -Oz -o ${sample}.dv.norm.vcf.gz
  bcftools index -t ${sample}.dv.norm.vcf.gz
  bcftools view -f PASS,. ${sample}.dv.norm.vcf.gz -Oz -o ${sample}.dv.norm.pass.vcf.gz
  bcftools index -t ${sample}.dv.norm.pass.vcf.gz
  """
}

process NormalizeHC {
  tag "$sample"
  //publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(vcf)
    path ref_fasta
    path ref_fai
  
  output:
    tuple val(sample), path("${sample}.hc.norm.pass.vcf.gz"), path("${sample}.hc.norm.pass.vcf.gz.tbi")
  
  script:
  """
  #bcftools norm -f ${ref_fasta} ${vcf} -Oz -o ${sample}.hc.norm.vcf.gz
  bcftools norm -m -any ${vcf} -Oz -o ${sample}.hc.norm.vcf.gz
  bcftools index -t ${sample}.hc.norm.vcf.gz
  bcftools view -f PASS,. ${sample}.hc.norm.vcf.gz -Oz -o ${sample}.hc.norm.pass.vcf.gz
  bcftools index -t ${sample}.hc.norm.pass.vcf.gz
  """
}

process SplitByCallerOverlap {
  tag "$sample"
  //publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(dv_vcf), path(dv_tbi), path(hc_vcf), path(hc_tbi)
  
  output:
    tuple val(sample), 
      path("isec_${sample}/0000.vcf"), 
      path("isec_${sample}/0001.vcf"), 
      path("isec_${sample}/0002.vcf")
  
  script:
  """
  mkdir -p isec_${sample}
  bcftools isec -p isec_${sample} ${dv_vcf} ${hc_vcf}
  """
}

process AnnotateDVOnly {
  tag "$sample"
  //publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(dv_only_vcf)
  
  output:
    tuple val(sample), path("${sample}.dv_only.tag.vcf.gz"), path("${sample}.dv_only.tag.vcf.gz.tbi")
  
  script:
  """
  echo '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Which callers support this variant (DV,HC)">' > callers.hdr
  
  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\tDV\\n' ${dv_only_vcf} > dv_callers.tsv
  bgzip -f dv_callers.tsv
  tabix -s 1 -b 2 -e 2 dv_callers.tsv.gz
  
  bcftools annotate \
    -h callers.hdr \
    -a dv_callers.tsv.gz \
    -c CHROM,POS,REF,ALT,CALLERS \
    ${dv_only_vcf} \
    -Oz -o ${sample}.dv_only.tag.vcf.gz
  
  bcftools index -t ${sample}.dv_only.tag.vcf.gz
  """
}

process AnnotateHCOnly {
  tag "$sample"
  //publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(hc_only_vcf)
  
  output:
    tuple val(sample), path("${sample}.hc_only.tag.vcf.gz"), path("${sample}.hc_only.tag.vcf.gz.tbi")
  
  script:
  """
  echo '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Which callers support this variant (DV,HC)">' > callers.hdr
  
  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\tHC\\n' ${hc_only_vcf} > hc_callers.tsv
  bgzip -f hc_callers.tsv
  tabix -s 1 -b 2 -e 2 hc_callers.tsv.gz
  
  bcftools annotate \
    -h callers.hdr \
    -a hc_callers.tsv.gz \
    -c CHROM,POS,REF,ALT,CALLERS \
    ${hc_only_vcf} \
    -Oz -o ${sample}.hc_only.tag.vcf.gz
  
  bcftools index -t ${sample}.hc_only.tag.vcf.gz
  """
}

process AnnotateShared {
  tag "$sample"
  //publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(shared_vcf)
  
  output:
    tuple val(sample), path("${sample}.both.tag.vcf.gz"), path("${sample}.both.tag.vcf.gz.tbi")
  
  script:
  """
  echo '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Which callers support this variant (DV,HC)">' > callers.hdr
  
  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\tDV,HC\\n' ${shared_vcf} > both_callers.tsv
  bgzip -f both_callers.tsv
  tabix -s 1 -b 2 -e 2 both_callers.tsv.gz
  
  bcftools annotate \
    -h callers.hdr \
    -a both_callers.tsv.gz \
    -c CHROM,POS,REF,ALT,CALLERS \
    ${shared_vcf} \
    -Oz -o ${sample}.both.tag.vcf.gz
  
  bcftools index -t ${sample}.both.tag.vcf.gz
  """
}

process BuildConsensusVCF {
  tag "$sample"
  publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), 
      path(dv_only_vcf), path(dv_only_tbi),
      path(hc_only_vcf), path(hc_only_tbi),
      path(both_vcf), path(both_tbi)
  
  output:
    tuple val(sample), path("${sample}.consensus.vcf.gz"), path("${sample}.consensus.vcf.gz.tbi")
  
  script:
  """
  bcftools concat -a \
    ${both_vcf} \
    ${dv_only_vcf} \
    ${hc_only_vcf} \
    -Oz -o ${sample}.consensus.vcf.gz
  
  bcftools index -t ${sample}.consensus.vcf.gz
  """
}

process FilterConsensusVCF {
  tag "$sample"
  publishDir "${params.outdir}/${sample}/consensus", mode: 'copy'
  
  input:
    tuple val(sample), path(consensus_vcf), path(consensus_tbi)
  
  output:
    tuple val(sample), path("${sample}.consensus.filtered.vcf.gz")
  
  script:
  """
  bcftools view \
    -i '((INFO/CALLERS="DV,HC") && FORMAT/DP>=10 && FORMAT/GQ>=10) || \
        ((INFO/CALLERS="HC") && QUAL>=30 && FORMAT/DP>=10 && FORMAT/GQ>=20) || \
        ((INFO/CALLERS="DV") && QUAL>=10 && FORMAT/DP>=10 && FORMAT/GQ>=20)' \
    ${consensus_vcf} \
    -Oz -o ${sample}.consensus.filtered.vcf.gz
  
  bcftools index -t ${sample}.consensus.filtered.vcf.gz
  """
}

/********************  CONSENSUS WORKFLOW  ********************/

workflow CONSENSUS_CALLING {
  take:
    dv_vcf_ch    // (sample, dv_vcf)
    hc_vcf_ch    // (sample, hc_vcf)
    ref_fasta_ch // reference fasta
    ref_fai_ch   // reference fasta index
  
  main:
    dv_labeled = dv_vcf_ch.map { s, v -> tuple(s, 'DV', v) }
    hc_labeled = hc_vcf_ch.map { s, v -> tuple(s, 'HC', v) }
    both_callers_ch = dv_labeled.mix(hc_labeled)

    // Reheader ONCE for both callers
    CONS_REHEADER_VCF(both_callers_ch)

    // Split back to DV and HC channels (drop caller and tbi for NormalizeDV/HC inputs)
    dv_fix_ch = CONS_REHEADER_VCF.out
      .filter { s, caller, vcf, tbi -> caller == 'DV' }
      .map    { s, caller, vcf, tbi -> tuple(s, vcf) }

    hc_fix_ch = CONS_REHEADER_VCF.out
      .filter { s, caller, vcf, tbi -> caller == 'HC' }
      .map    { s, caller, vcf, tbi -> tuple(s, vcf) }

    // Normalize both VCFs
    NormalizeDV(dv_fix_ch, ref_fasta_ch, ref_fai_ch)
    NormalizeHC(hc_fix_ch, ref_fasta_ch, ref_fai_ch)
    
    // Join DV and HC normalized VCFs by sample
    combined_ch = NormalizeDV.out
      .join(NormalizeHC.out)
      .map { sample, dv_vcf, dv_tbi, hc_vcf, hc_tbi -> 
        tuple(sample, dv_vcf, dv_tbi, hc_vcf, hc_tbi) 
      }
    
    // Split by caller overlap
    SplitByCallerOverlap(combined_ch)
    
    // Annotate each category
    dv_only_ch = SplitByCallerOverlap.out.map { s, dv, hc, both -> tuple(s, dv) }
    hc_only_ch = SplitByCallerOverlap.out.map { s, dv, hc, both -> tuple(s, hc) }
    shared_ch  = SplitByCallerOverlap.out.map { s, dv, hc, both -> tuple(s, both) }
    
    AnnotateDVOnly(dv_only_ch)
    AnnotateHCOnly(hc_only_ch)
    AnnotateShared(shared_ch)
    
    // Combine all annotated VCFs
    consensus_input_ch = AnnotateDVOnly.out
      .join(AnnotateHCOnly.out)
      .join(AnnotateShared.out)
      .map { sample, dv_vcf, dv_tbi, hc_vcf, hc_tbi, both_vcf, both_tbi ->
        tuple(sample, dv_vcf, dv_tbi, hc_vcf, hc_tbi, both_vcf, both_tbi)
      }
    
    // Build consensus
    BuildConsensusVCF(consensus_input_ch)
    
    // Filter consensus
    FilterConsensusVCF(BuildConsensusVCF.out)
  
  emit:
    consensus_vcf = FilterConsensusVCF.out
}
