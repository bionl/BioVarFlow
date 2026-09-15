// modules/exomiser.nf
nextflow.enable.dsl=2

// ── Exomiser: phenotype-driven ranking over the WHOLE callset ────────────────
//
// This runs on the UNFILTERED consensus VCF, not consensus.filtered, and not
// the panel-filtered VCF the reports are built from. Those two filters exist to
// make the reporting sheets trustworthy -- BedFilterVCF restricts to the 285
// panel genes and drops 99% of calls, and the tiered DP/GQ/QUAL filter halves
// what is left. Exomiser's whole purpose here is to surface findings OUTSIDE
// the panel, and it applies its own frequency and pathogenicity filters before
// ranking, so pre-filtering only removes candidates before it can see them.
//
// One batch task for the whole run rather than one per sample: the Exomiser
// data directory is the expensive part, and staging it once amortises across
// every sample.
params.exomiser_data_gs  = params.exomiser_data_gs  ?: "gs://vep-data-vaic/exomiser"
params.exomiser_version  = params.exomiser_version  ?: "2512"
params.phenotypes        = params.phenotypes        ?: null
// Not `?: true` -- Groovy elvis treats false as unset, so
// `--run_exomiser false` would be silently flipped back on.
if (!params.containsKey("run_exomiser")) params.run_exomiser = true

process EXOMISER_BATCH {
  tag { "${samples.size()} sample(s)" }
  publishDir "${params.outdir}/exomiser", mode: 'copy'

  input:
    val  samples        // [[sample:, pheno_sex:, hpo:], ...] -- order matches vcfs
    path vcfs           // staged consensus VCFs, named <sample>.consensus.vcf.gz
    path tbis

  output:
    path "results/*", emit: reports
    path "batch.txt",  emit: batch_file

  script:
    // Phenopacket v2. Only subject.sex and phenotypicFeatures affect scoring --
    // labels and age do not, so the CSV carries neither.
    def stamp = new Date().format("yyyy-MM-dd'T'HH:mm:ss'Z'", TimeZone.getTimeZone('UTC'))
    def writes = samples.collect { m ->
        def features = m.hpo.split(';')
                            .findAll { it?.trim() }
                            .collect { "  - type: { id: \"${it.trim()}\" }" }
                            .join('\n')
        def ppkt = """id: ${m.sample}
subject:
  id: ${m.sample}
  sex: ${m.pheno_sex}
phenotypicFeatures:
${features}
metaData:
  created: "${stamp}"
  createdBy: "biovarflow"
  phenopacketSchemaVersion: "2.0"
"""
        "cat > ${m.sample}.phenopacket.yml <<'PPKT'\n${ppkt}PPKT"
    }.join('\n')

    def batch = samples.collect { m ->
        "--sample ${m.sample}.phenopacket.yml --vcf ${m.sample}.consensus.vcf.gz " +
        "--assembly GRCh38 --preset exome --output-directory results " +
        "--output-filename ${m.sample}_exomiser --output-format HTML,JSON,TSV_GENE,TSV_VARIANT"
    }.join('\n')
  """
  set -euo pipefail
  mkdir -p ${params.exomiser_data_dir} results

  # The data directory is pulled rather than staged as a Nextflow `path` input:
  # 2512_hg38 plus 2512_phenotype is far too large to localise per task, and one
  # batch task means this happens once per run.
  gcloud storage rsync -r ${params.exomiser_data_gs} ${params.exomiser_data_dir}

  export EXOMISER_DATA_DIRECTORY=${params.exomiser_data_dir}
  export EXOMISER_HG38_DATA_VERSION=${params.exomiser_version}
  export EXOMISER_PHENOTYPE_DATA_VERSION=${params.exomiser_version}
  export JAVA_TOOL_OPTIONS=-Xmx${task.memory.toGiga().intdiv(2)}g

${writes}

  cat > batch.txt <<'BATCH'
${batch}
BATCH

  # --dry-run validates every line before any analysis starts, so a malformed
  # phenopacket fails in seconds instead of after the first sample completes.
  exomiser batch --dry-run batch.txt
  exomiser batch batch.txt
  """
}


workflow EXOMISER {
  take:
    raw_consensus_ch   // (meta, vcf, tbi) -- UNFILTERED consensus

  main:
    if (!params.phenotypes) {
      error "❌ --run_exomiser requires --phenotypes <csv>\n" +
            "   Columns: sample,sex,hpo\n" +
            "   sex is MALE / FEMALE / UNKNOWN_SEX; hpo is ';'-separated HPO IDs.\n" +
            "   e.g.  IQMM,MALE,HP:0004808;HP:0002863"
    }

    pheno_ch = Channel.fromPath(params.phenotypes, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        def id = (row.sample ?: '').trim()
        if (!id) error "❌ --phenotypes row with no sample id: ${row}"
        tuple(id, (row.sex ?: 'UNKNOWN_SEX').trim(), (row.hpo ?: '').trim())
      }

    // remainder:true so a sample missing from the CSV surfaces as a warning
    // rather than vanishing from the run.
    joined_ch = raw_consensus_ch
      .map { meta, vcf, tbi -> tuple(meta.sample, vcf, tbi) }
      .join(pheno_ch, remainder: true)
      .branch {
        ready:   it[1] != null && it[3] != null && it[4]
        no_vcf:  it[1] == null
        no_pheno: true
      }

    joined_ch.no_pheno.view { "⚠️  Exomiser: no HPO terms for ${it[0]} — sample skipped." }
    joined_ch.no_vcf.view   { "⚠️  Exomiser: ${it[0]} in --phenotypes but not in the run — ignored." }

    // Collect once so `samples` and `vcfs` stay index-aligned.
    rows_ch = joined_ch.ready.toList().filter { it.size() > 0 }

    EXOMISER_BATCH(
      rows_ch.map { rows -> rows.collect { s, vcf, tbi, sex, hpo ->
                              [ sample: s, pheno_sex: sex, hpo: hpo ] } },
      rows_ch.map { rows -> rows.collect { s, vcf, tbi, sex, hpo -> vcf } },
      rows_ch.map { rows -> rows.collect { s, vcf, tbi, sex, hpo -> tbi } }
    )

  emit:
    reports = EXOMISER_BATCH.out.reports
}
