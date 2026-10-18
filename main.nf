nextflow.enable.dsl=2

// ═══════════════════════════════════════════════════════════════════════════
// IMPORTS
// ═══════════════════════════════════════════════════════════════════════════

include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION } \
  from './external/sarek/subworkflows/local/utils_nfcore_sarek_pipeline'

include { NFCORE_SAREK } \
  from './external/sarek/main.nf'

include { POST_SAREK } \
  from './modules/vep.nf'

include { FilterConsensusVCF } \
    from './modules/consensus.nf'
include { CONSENSUS_CALLING } \
  from './modules/consensus.nf'

include { EXOMISER } \
  from './modules/exomiser.nf'

include { DB_QC_EXPORT } \
  from './modules/db_qc_export.nf'

include { MANIFEST } \
  from './modules/manifest.nf'

include { MUTECT2_RESCUE } \
  from './modules/mutect2_rescue.nf'

// ═══════════════════════════════════════════════════════════════════════════
// PARAMETERS & VALIDATION
// ═══════════════════════════════════════════════════════════════════════════

// Default parameters
params.input                  = params.input ?: params.samplesheet
params.outdir                 = params.outdir ?: params.output
params.bed                    = params.bed ?: "${workflow.projectDir}/data/ACMG_SF_MANE_exons_50bp.bed"
//params.run_variant_calling    = params.run_variant_calling instanceof Boolean ? params.run_variant_calling : true
params.create_consensus       = params.create_consensus instanceof Boolean ? params.create_consensus : true
params.run_db_qc              = params.run_db_qc instanceof Boolean ? params.run_db_qc : true
// Run identifier propagated to the run-output manifest (consumed by the
// downstream DB ingestion pipeline). Falls back to workflow.runName at
// manifest-build time if left null.
params.run_id                 = params.run_id ?: null
// Reference for consensus normalisation. Deliberately NOT defaulted here:
// params.fasta is assigned by external/sarek/main.nf at include time and is not
// reliably populated while this script-level params block runs (and an
// externally-supplied null cannot be overwritten). Resolved instead by
// resolveRefFasta() at the point of use, which errors clearly if unset.
params.ref_fasta              = params.ref_fasta ?: null

// Which half of the pipeline to run.
//
//   calling    FASTQ/BAM -> Sarek -> DeepVariant + HaplotypeCaller -> consensus
//              VCF. Stops there. Needs NO panel inputs at all: no --bed, no
//              --sf_genes, no --phenotypes, no VEP resources. This is the
//              deliverable for a customer who wants variants, not a report.
//
//   reporting  consensus VCF + RAW BAM -> annotation, Exomiser and the Excel
//              report. Takes --reporting_input, a CSV of
//              sample,raw_vcf,raw_tbi,bam,bai.
//
//              The RAW consensus is the handover, not the filtered one: it is
//              what Exomiser needs (it ranks genome-wide, and the filtered VCF
//              is panel-restricted) and what strand_bias.py needs for off-panel
//              REF/ALT and AD. The filtered VCF is re-derived here by
//              FilterConsensusVCF, so the tiered CALLERS thresholds live in one
//              place and the customer hands over one VCF rather than two.
//
//              The BAM must be the RAW alignment, not a panel-filtered one:
//              BedFilterBAM does the panel restriction itself, and
//              StrandBiasPileup deliberately reads the unfiltered BAM so
//              off-panel variants are not left blank.
//
//   full       both, as a single run. The default, so existing invocations are
//              unaffected.
params.stage                  = params.stage ?: 'full'
params.reporting_input        = params.reporting_input ?: null

if (!(params.stage in ['full', 'calling', 'reporting'])) {
    error "❌ --stage must be one of: full, calling, reporting (got '${params.stage}')"
}
if (params.stage == 'reporting' && !params.reporting_input) {
    error "❌ --stage reporting requires --reporting_input (CSV: sample,raw_vcf,raw_tbi,bam,bai[,assay])"
}
//params.vep_fasta              = params.vep_fasta ?: params.vep_fasta

// Validate required parameters
//if (params.run_variant_calling) {
//    if (!params.input)  error "❌ Missing --input (samplesheet CSV) when run_variant_calling=true"
//    if (!params.outdir) error "❌ Missing --outdir when run_variant_calling=true"
//    if (params.create_consensus && !params.ref_fasta) {
//        error "❌ Missing --vep_fasta when create_consensus=true"
//    }
//} else {
//    if (!params.post_samplesheet && !params.variant_calling_outdir)
//        error "❌ When run_variant_calling=false provide either --post_samplesheet or --variant_calling_outdir"
//    
//    if (params.post_samplesheet && params.variant_calling_outdir)
//        error "❌ Cannot provide both --post_samplesheet and --variant_calling_outdir. Choose one."
//}

// ═══════════════════════════════════════════════════════════════════════════
// HELPER FUNCTIONS
// ═══════════════════════════════════════════════════════════════════════════

def isGcsPath(path) {
    return path.toString().startsWith('gs://')
}

// Resolve the reference used for bcftools norm -f in the consensus subworkflow.
// Must be the assembly the BAMs were aligned to (GATK, chr-prefixed contigs) --
// NOT params.vep_fasta, which is the Ensembl build named 1/2/.../MT and fails on
// every record. Called from inside a workflow body so params.fasta is populated.
//
// Fallback chain, most-explicit first:
//   1. --ref_fasta                   explicit override
//   2. params.fasta                  set by external/sarek/main.nf; not reliably
//                                    visible from this binding, never relied on alone
//   3. params.genomes[genome].fasta  the igenomes map, loaded via includeConfig
//                                    at config-parse time and always available here
def resolveRefFasta() {
    def f = params.ref_fasta ?: params.fasta
    if (!f && params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        f = params.genomes[params.genome].fasta
    }
    if (!f) {
        error "❌ No reference genome for consensus normalisation.\n" +
              "   Tried --ref_fasta, params.fasta, and params.genomes[${params.genome}].fasta.\n" +
              "   genome=${params.genome}  igenomes_ignore=${params.igenomes_ignore}  " +
              "genomes_loaded=${params.genomes ? params.genomes.size() : 0}\n" +
              "   Set --ref_fasta to the GATK assembly the BAMs were aligned against " +
              "(chr-prefixed contigs)."
    }
    return f
}

def validateBedFile() {
    def bedFile = params.bed ? file(params.bed) : null
    if (!bedFile?.exists()) {
        error "❌ BED file not found: ${params.bed}"
    }
    return bedFile
}

// ═══════════════════════════════════════════════════════════════════════════
// WORKFLOWS
// ═══════════════════════════════════════════════════════════════════════════

workflow COLLECT_VARIANT_CALLING_OUTPUTS {
    take:
        trigger    // Completion signal channel
        outdir     // Output directory to search

    main:
        def isGCS = isGcsPath(outdir)
        
        // Collect DeepVariant VCFs
        dv_vcf_ch = trigger
            .flatMap { 
                file("${outdir}/variant_calling/deepvariant/*/*.vcf.gz", checkIfExists: !isGCS)
            }
            .filter { vcf -> 
                vcf.name.endsWith('.vcf.gz') && 
                !vcf.name.contains('.g.vcf.gz') && 
                !vcf.name.endsWith('.tbi') 
            }
            .map { vcf -> tuple(vcf.parent.name, vcf) }
        
        // Collect HaplotypeCaller VCFs
        hc_vcf_ch = trigger
            .flatMap { 
                file("${outdir}/variant_calling/haplotypecaller/*/*filtered.vcf.gz", checkIfExists: !isGCS)
            }
            .filter { vcf -> 
                vcf.name.endsWith('.vcf.gz') && 
                !vcf.name.contains('.g.vcf.gz') && 
                !vcf.name.endsWith('.tbi') 
            }
            .map { vcf -> tuple(vcf.parent.name, vcf) }

        // Collect BAMs with BAI
        bam_ch = trigger
            .flatMap { 
                file("${outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: !isGCS)
            }
            .map { bam -> 
                def sample = bam.parent.name
                def bamPath = bam.toString()
                def baiPath = "${bamPath}.bai"
                
                def bai
                if (isGCS) {
                    bai = file(baiPath, checkIfExists: false)
                } else {
                    bai = file(baiPath)
                    if (!bai.exists()) {
                        bai = file("${bam.parent}/${bam.baseName}.bai")
                        if (!bai.exists()) {
                            error "❌ BAM index not found for ${bam}"
                        }
                    }
                }
                
                tuple(sample, bam, bai)
            }

        // Collect Sarek QC reports (reused by DB_QC_EXPORT to avoid re-running tools)

        // Alignment QC: prefer samtools stats, fall back to markdup metrics
        samtools_stats_raw = trigger
            .flatMap {
                file("${outdir}/reports/samtools/*/*.md.cram.stats", checkIfExists: !isGCS)
            }
            .map { f -> tuple(f.parent.name, f, 'samtools_stats') }

        markdup_raw = trigger
            .flatMap {
                file("${outdir}/reports/markduplicates/*/*.md.cram.metrics", checkIfExists: !isGCS)
            }
            .map { f -> tuple(f.parent.name, f, 'markdup_metrics') }

        // Merge both sources per sample; prefer samtools stats when available
        alignment_qc_ch = samtools_stats_raw
            .mix(markdup_raw)
            .groupTuple(by: 0)
            .map { sample, files, types ->
                def idx = types.indexOf('samtools_stats')
                if (idx >= 0) {
                    tuple(sample, files[idx], 'samtools_stats')
                } else {
                    tuple(sample, files[0], 'markdup_metrics')
                }
            }

        mosdepth_summary_ch = trigger
            .flatMap {
                file("${outdir}/reports/mosdepth/*/*.md.mosdepth.summary.txt", checkIfExists: !isGCS)
            }
            .map { summary -> tuple(summary.parent.name, summary) }

        mosdepth_dist_ch = trigger
            .flatMap {
                file("${outdir}/reports/mosdepth/*/*.md.mosdepth.global.dist.txt", checkIfExists: !isGCS)
            }
            .map { dist -> tuple(dist.parent.name, dist) }

    emit:
        dv_vcf           = dv_vcf_ch
        hc_vcf           = hc_vcf_ch
        bam              = bam_ch
        alignment_qc     = alignment_qc_ch       // tuple(sample, file, type)
        mosdepth_summary = mosdepth_summary_ch
        mosdepth_dist    = mosdepth_dist_ch
}

workflow RUN_FROM_VARIANT_CALLING_OUTDIR {
    take:
        variant_calling_outdir
        bed_ch

    main:
        def isGCS = isGcsPath(variant_calling_outdir)
        
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Using existing variant calling results                    ║
        ║  Location: ${variant_calling_outdir}
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

        // Collect VCFs (use consensus or single caller depending on params)
        if (params.create_consensus) {
            // Collect DeepVariant VCFs
            dv_vcf_ch = Channel
                .fromPath("${variant_calling_outdir}/variant_calling/deepvariant/*/*.vcf.gz", checkIfExists: !isGCS)
                .filter { vcf -> 
                    vcf.name.endsWith('.vcf.gz') && 
                    !vcf.name.contains('.g.vcf.gz') && 
                    !vcf.name.endsWith('.tbi') 
                }
                .map { vcf -> 
                    def sample = vcf.parent.name
                    tuple(sample, vcf) 
                }
            
            // Collect HaplotypeCaller VCFs
            hc_vcf_ch = Channel
                .fromPath("${variant_calling_outdir}/variant_calling/haplotypecaller/*/*.vcf.gz", checkIfExists: !isGCS)
                .filter { vcf -> 
                    vcf.name.endsWith('.vcf.gz') && 
                    !vcf.name.contains('.g.vcf.gz') && 
                    !vcf.name.endsWith('.tbi') 
                }
                .map { vcf -> 
                    def sample = vcf.parent.name
                    tuple(sample, vcf) 
                }
            
            // Create reference channels
            def _refFasta = resolveRefFasta()
            ref_fasta_ch = Channel.value(file(_refFasta))
            ref_fai_ch = Channel.value(file(_refFasta + ".fai"))
            
            // Run consensus calling
            CONSENSUS_CALLING(dv_vcf_ch, hc_vcf_ch, ref_fasta_ch, ref_fai_ch)
            vcf_ch = CONSENSUS_CALLING.out.consensus_vcf
            
        } else {
            // Use single caller VCF (default to DeepVariant or configurable)
            vcf_ch = Channel
                .fromPath("${variant_calling_outdir}/variant_calling/*/*/*.vcf.gz", checkIfExists: !isGCS)
                .filter { vcf -> 
                    vcf.name.endsWith('.vcf.gz') && 
                    !vcf.name.contains('.g.vcf.gz') && 
                    !vcf.name.endsWith('.tbi') 
                }
                .map { vcf -> 
                    def sample = vcf.parent.name
                    tuple(sample, vcf) 
                }
        }

        // Collect BAMs with BAI
        bam_ch = Channel
            .fromPath("${variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam", checkIfExists: !isGCS)
            .map { bam -> 
                def sample = bam.parent.name
                def bamPath = bam.toString()
                def baiPath = "${bamPath}.bai"
                
                def bai
                if (isGCS) {
                    bai = file(baiPath, checkIfExists: false)
                } else {
                    bai = file(baiPath)
                    if (!bai.exists()) {
                        bai = file("${bam.parent}/${bam.baseName}.bai")
                        if (!bai.exists()) {
                            error "❌ BAM index not found for ${bam}"
                        }
                    }
                }
                
                tuple(sample, bam, bai) 
            }
        
        // Debug output
        vcf_ch.view { s, v -> "📄 VCF -> ${s} :: ${v.name}" }
        bam_ch.view { s, a, i -> "🧬 BAM -> ${s} :: ${a.name}" }
        
        // Safety checks
        vcf_ch
            .count()
            .subscribe { count ->
                if (count == 0) {
                    error "❌ No VCFs found"
                }
                log.info "✓ Found ${count} VCF file(s)"
            }
        
        bam_ch
            .count()
            .subscribe { count ->
                if (count == 0) {
                    error "❌ No BAMs found in ${variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam"
                }
                log.info "✓ Found ${count} BAM file(s)"
            }

        // Run post-processing
        POST_SAREK(Channel.empty(), Channel.empty(), vcf_ch, bam_ch, bed_ch)
}

workflow RUN_FROM_POST_SAMPLESHEET {
    take:
        post_samplesheet
        bed_ch

    main:
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Using custom post-samplesheet                             ║
        ║  File: ${post_samplesheet}
        ║  Note: Consensus calling is skipped when using             ║
        ║        post-samplesheet (single VCF per sample expected)   ║
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

        // Parse samplesheet
        Channel
            .fromPath(post_samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def v = file(row.vcf)
                def b = file(row.bam)
                def bi = file(row.bai ?: "${b}.bai")
                
                // Only validate for non-GCS paths
                def isGCS = row.vcf.startsWith('gs://')
                if (!isGCS) {
                    if (!v.exists()) error "❌ VCF not found: ${v}"
                    if (!b.exists()) error "❌ BAM not found: ${b}"
                    if (!bi.exists()) {
                        bi = file("${b.parent}/${b.baseName}.bai")
                        if (!bi.exists()) error "❌ BAI not found for ${b}"
                    }
                }
                
                tuple(row.sample, v, b, bi)
            }
            .multiMap { sample, vcf, bam, bai ->
                vcf: tuple(sample, vcf)
                bam: tuple(sample, bam, bai)
            }
            .set { result }

        vcf_ch = result.vcf
        bam_ch = result.bam
        
        // Debug output
        vcf_ch.view { s, v -> "📄 VCF -> ${s} :: ${v.name}" }
        bam_ch.view { s, a, i -> "🧬 BAM -> ${s} :: ${a.name}" }

        // Run post-processing (no consensus for post-samplesheet)
        POST_SAREK(Channel.empty(), Channel.empty(), vcf_ch, bam_ch, bed_ch)
}

// Second half: consensus VCF + raw BAM -> annotation, Exomiser, report.
// Mirrors the tail of RUN_FULL_VARIANT_CALLING rather than re-deriving it, so
// the two entry points produce the same sheets from the same code.
workflow RUN_REPORTING {
    take:
        bed_ch

    main:
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Stage: REPORTING only                                     ║
        ║  Input: ${params.reporting_input}
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

        // CSV: sample,raw_vcf,raw_tbi,bam,bai,assay
        // raw_tbi and bai default to <file>.tbi / <file>.bai. Both are verified
        // for local paths; on GCS `file()` cannot cheaply confirm existence, so
        // a wrong path there fails at stage-in with a clear Nextflow error
        // rather than silently producing an empty channel.
        def rows_ch = Channel
            .fromPath(params.reporting_input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                if (!row.sample)  error "❌ --reporting_input row missing 'sample': ${row}"
                if (!row.raw_vcf) error "❌ --reporting_input row missing 'raw_vcf': ${row}"
                if (!row.bam)     error "❌ --reporting_input row missing 'bam': ${row}"

                def v  = file(row.raw_vcf)
                def vt = file(row.raw_tbi ?: "${row.raw_vcf}.tbi")
                def b  = file(row.bam)
                def bi = file(row.bai ?: "${row.bam}.bai")

                if (!row.raw_vcf.startsWith('gs://')) {
                    if (!v.exists())  error "❌ VCF not found: ${v}"
                    if (!vt.exists()) error "❌ VCF index not found: ${vt}"
                    if (!b.exists())  error "❌ BAM not found: ${b}"
                    if (!bi.exists()) error "❌ BAI not found: ${bi}"
                }
                // assay comes from the CSV, not params.input: assayMap is built
                // from the calling samplesheet, which a reporting-only run does
                // not have. Optional, defaults to NA like the full-run fallback.
                tuple(row.sample, v, vt, b, bi, (row.assay ?: 'NA'))
            }

        // The handover is the RAW consensus; re-derive the filtered VCF here so
        // the CALLERS-tiered thresholds are applied by the same process that
        // applies them in a full run.
        FilterConsensusVCF(rows_ch.map { s, v, vt, b, bi, assay -> tuple(s, v, vt) })

        // assay per sample, so vcf_with_meta_ch can carry the same meta MAP the
        // other channels do -- POST_SAREK joins on that map, and a mismatched
        // key silently matches nothing.
        def assay_ch = rows_ch.map { sample, v, vt, b, bi, assay -> tuple(sample, assay) }

        def vcf_with_meta_ch = FilterConsensusVCF.out
            .join(assay_ch)
            .map { sample, vcf, assay -> tuple([ sample: sample, assay: assay ], vcf) }

        def bam_with_meta_ch = rows_ch
            .map { sample, v, vt, b, bi, assay ->
                tuple([ sample: sample, assay: assay ], b, bi) }

        def raw_cons_keyed_ch = rows_ch
            .map { sample, v, vt, b, bi, assay ->
                tuple([ sample: sample, assay: assay ], v, vt) }

        exomiser_tsv_ch = Channel.empty()
        if (params.run_exomiser && params.phenotypes) {
            EXOMISER(rows_ch.map { sample, v, vt, b, bi, assay ->
                tuple([ sample: sample, assay: assay ], v, vt) })
            exomiser_tsv_ch = EXOMISER.out.reports
                .flatten()
                .filter { f -> f.name.endsWith('_exomiser.variants.tsv') }
                .map    { f -> tuple(f.name.replaceFirst(/_exomiser\.variants\.tsv$/, ''), f) }
        } else if (params.run_exomiser) {
            log.warn "Exomiser skipped: --stage reporting needs --phenotypes."
        }

        POST_SAREK(exomiser_tsv_ch, raw_cons_keyed_ch, vcf_with_meta_ch, bam_with_meta_ch, bed_ch)
}

workflow RUN_FULL_VARIANT_CALLING {
    take:
        bed_ch

    main:
        // Only the consensus branch produces one; declared up front so the
        // other branches leave a defined value.
        raw_consensus_ch = null

        if (params.somatic_mode) {
            log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Running somatic variant calling pipeline (Sarek/Mutect2)  ║
        ║  Tumor-only and tumor+normal pairs supported               ║
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()
        } else {
            log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Running full variant calling pipeline (Sarek)             ║
        ║  Consensus calling: ${params.create_consensus ? 'ENABLED' : 'DISABLED'}
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()
        }

        PIPELINE_INITIALISATION(
            params.version,
            params.validate_params,
            args,
            params.outdir,
            params.input,
            params.help,
            params.help_full,
            params.show_hidden,
        )

        def assayMap = [:]

        PIPELINE_INITIALISATION.out.samplesheet
            .map { row ->
                def meta = row[0]
                tuple(meta.sample.toString(), (meta.assay ?: 'NA').toString())
            }
            .toList()
            .subscribe { pairs ->
                assayMap = pairs.collectEntries { s, a -> [(s): a] }
                log.info "✓ Loaded assay metadata for ${assayMap.size()} sample(s)"
                assayMap.each { k, v -> log.info "  ${k} -> ${v}" }
            }

        NFCORE_SAREK(PIPELINE_INITIALISATION.out.samplesheet)

        PIPELINE_COMPLETION(
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir,
            params.monochrome_logs,
            params.hook_url,
            NFCORE_SAREK.out.multiqc_report
        )

        COLLECT_VARIANT_CALLING_OUTPUTS(
            NFCORE_SAREK.out.multiqc_report,
            params.outdir
        )

        if (params.somatic_mode) {
            def isGCS = isGcsPath(params.outdir)
            def mutect2_vcf_ch = NFCORE_SAREK.out.multiqc_report
                .flatMap {
                    file("${params.outdir}/variant_calling/mutect2/*/*.mutect2.filtered.vcf.gz", checkIfExists: !isGCS)
                }
                .filter { vcf -> vcf.name.endsWith('.vcf.gz') && !vcf.name.endsWith('.tbi') }
                .map { vcf -> tuple(vcf.parent.name, vcf) }

            def rescue_script_ch = Channel.value(file("${params.scriptdir}/mutect2_rescue.py"))
            MUTECT2_RESCUE(mutect2_vcf_ch, rescue_script_ch)
            final_vcf_ch = MUTECT2_RESCUE.out.map { sample, vcf, tbi -> tuple(sample, vcf) }
        } else if (params.create_consensus) {
            def _refFasta = resolveRefFasta()
            ref_fasta_ch = Channel.value(file(_refFasta))
            ref_fai_ch   = Channel.value(file(_refFasta + ".fai"))

            CONSENSUS_CALLING(
                COLLECT_VARIANT_CALLING_OUTPUTS.out.dv_vcf,
                COLLECT_VARIANT_CALLING_OUTPUTS.out.hc_vcf,
                ref_fasta_ch,
                ref_fai_ch
            )

            final_vcf_ch = CONSENSUS_CALLING.out.consensus_vcf
            // Exomiser ranks genome-wide, so it takes the UNFILTERED consensus.
            // Everything downstream is panel-restricted and would hide exactly
            // the off-panel findings Exomiser is run to surface.
            raw_consensus_ch = CONSENSUS_CALLING.out.consensus_vcf_raw
        } else {
            final_vcf_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.dv_vcf
        }

        def vcf_with_meta_ch = final_vcf_ch.map { sample, vcf ->
            def meta = [ sample: sample, assay: assayMap.get(sample, 'NA') ]
            tuple(meta, vcf)
        }

        // When starting from BAM files, Sarek writes markdup outputs to
        // preprocessing/markduplicates/ not preprocessing/mapped/, so we read
        // BAMs directly from the input samplesheet instead.
        def bam_with_meta_ch
        if (params.input_bam) {
            bam_with_meta_ch = Channel
                .fromPath(params.input_bam, checkIfExists: true)
                .splitCsv(header: true)
                .map { row ->
                    def meta = [ sample: row.sample, assay: assayMap.get(row.sample, 'NA') ]
                    tuple(meta, file(row.bam), file(row.bai))
                }
        } else {
            bam_with_meta_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.bam.map { sample, bam, bai ->
                def meta = [ sample: sample, assay: assayMap.get(sample, 'NA') ]
                tuple(meta, bam, bai)
            }
        }

        if (params.run_db_qc) {
            def alignment_qc_meta_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.alignment_qc
                .map { sample, f, type ->
                    tuple([ sample: sample, assay: assayMap.get(sample, 'NA') ], f, type)
                }
            def mosdepth_summary_meta_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.mosdepth_summary
                .map { sample, summary ->
                    tuple([ sample: sample, assay: assayMap.get(sample, 'NA') ], summary)
                }
            def mosdepth_dist_meta_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.mosdepth_dist
                .map { sample, dist ->
                    tuple([ sample: sample, assay: assayMap.get(sample, 'NA') ], dist)
                }

            DB_QC_EXPORT(
                vcf_with_meta_ch,
                alignment_qc_meta_ch,
                mosdepth_summary_meta_ch,
                mosdepth_dist_meta_ch
            )

            MANIFEST(
                vcf_with_meta_ch,
                bam_with_meta_ch,
                DB_QC_EXPORT.out.qc_json
            )
        } else {
            log.warn "params.run_db_qc=false → skipping run_output_manifest.tsv (QC Gate JSON is required to populate qc_status / qc_recommendation)."
        }

        // --stage calling stops here: the consensus VCF, the QC export and the
        // run manifest are the deliverable. Everything below is panel-scoped
        // annotation and reporting.
        if (params.stage == 'calling') {
            log.info """
            ╔════════════════════════════════════════════════════════════╗
            ║  Stage: CALLING only                                       ║
            ║  Consensus VCFs are in <outdir>/<sample>/consensus/        ║
            ║  To report on them later, rerun with:                      ║
            ║    --stage reporting --reporting_input <csv>               ║
            ╚════════════════════════════════════════════════════════════╝
            """.stripIndent()
            return
        }

        // Exomiser runs BEFORE POST_SAREK so its per-sample TSVs can become the
        // 'Exomiser' and 'Prioritised' tabs. EXOMISER_BATCH carries
        // errorStrategy 'ignore' and POST_SAREK substitutes a NO_FILE
        // placeholder, so a failed or skipped Exomiser costs those two tabs and
        // nothing else.
        exomiser_tsv_ch = Channel.empty()
        if (params.run_exomiser && raw_consensus_ch != null && params.phenotypes) {
            EXOMISER(
                raw_consensus_ch.map { sample, vcf, tbi ->
                    tuple([ sample: sample, assay: assayMap.get(sample, 'NA') ], vcf, tbi)
                }
            )
            exomiser_tsv_ch = EXOMISER.out.reports
                .flatten()
                .filter { f -> f.name.endsWith('_exomiser.variants.tsv') }
                .map    { f -> tuple(f.name.replaceFirst(/_exomiser\.variants\.tsv$/, ''), f) }
        } else if (params.run_exomiser) {
            log.warn "Exomiser skipped: needs --create_consensus and --phenotypes."
        }

        // Strand bias tests off-panel variants too and needs the unfiltered
        // consensus for their REF/ALT and AD. Keyed by meta to match vcf_ch.
        raw_cons_keyed_ch = raw_consensus_ch == null ? Channel.empty()
            : raw_consensus_ch.map { sample, vcf, tbi ->
                  tuple([ sample: sample, assay: assayMap.get(sample, 'NA') ], vcf, tbi) }

        // POST_SAREK (VEP annotation) is germline-only — skip in somatic mode
        if (!params.somatic_mode) {
            POST_SAREK(exomiser_tsv_ch, raw_cons_keyed_ch, vcf_with_meta_ch, bam_with_meta_ch, bed_ch)
        }
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow {

    // --stage calling needs no panel BED at all, so do not demand one. Both
    // other stages do: validateBedFile() errors clearly if it is missing.
    if (params.stage == 'calling') {
        RUN_FULL_VARIANT_CALLING(Channel.empty())
    }
    else if (params.stage == 'reporting') {
        RUN_REPORTING(Channel.value(validateBedFile()))
    }
    else {
        RUN_FULL_VARIANT_CALLING(Channel.value(validateBedFile()))
    }
    // Route to appropriate sub-workflow
    //if (params.variant_calling_outdir) {
    //    RUN_FROM_VARIANT_CALLING_OUTDIR(params.variant_calling_outdir, bed_ch)
    //} 
    //else if (params.post_samplesheet) {
    //    RUN_FROM_POST_SAMPLESHEET(params.post_samplesheet, bed_ch)
    //} 
}

// ═══════════════════════════════════════════════════════════════════════════
// WORKFLOW COMPLETION
// ═══════════════════════════════════════════════════════════════════════════

workflow.onComplete {
    log.info """
    ╔════════════════════════════════════════════════════════════╗
    ║  Pipeline completed!                                       ║
    ║  Status: ${workflow.success ? '✓ SUCCESS' : '✗ FAILED'}
    ║  Duration: ${workflow.duration}
    ║  Results: ${params.outdir}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error """
    ╔════════════════════════════════════════════════════════════╗
    ║  ✗ Pipeline failed                                         ║
    ║  Error: ${workflow.errorMessage}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
