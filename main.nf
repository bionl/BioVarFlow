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

include { CONSENSUS_CALLING } \
  from './modules/consensus.nf'

include { DB_QC_EXPORT } \
  from './modules/db_qc_export.nf'

include { MANIFEST } \
  from './modules/manifest.nf'

// ═══════════════════════════════════════════════════════════════════════════
// PARAMETERS & VALIDATION
// ═══════════════════════════════════════════════════════════════════════════

// Default parameters
params.input                  = params.input ?: params.samplesheet
params.outdir                 = params.outdir ?: params.output
params.bed                    = params.bed ?: "${workflow.projectDir}/data/annotated_merged_MANE_deduped.bed"
//params.run_variant_calling    = params.run_variant_calling instanceof Boolean ? params.run_variant_calling : true
params.create_consensus       = params.create_consensus instanceof Boolean ? params.create_consensus : true
params.run_db_qc              = params.run_db_qc instanceof Boolean ? params.run_db_qc : true
// Run identifier propagated to the run-output manifest (consumed by the
// downstream DB ingestion pipeline). Falls back to workflow.runName at
// manifest-build time if left null.
params.run_id                 = params.run_id ?: null
params.ref_fasta              = params.ref_fasta ?: params.vep_fasta
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
            ref_fasta_ch = Channel.value(file(params.ref_fasta))
            ref_fai_ch = Channel.value(file(params.ref_fasta + ".fai"))
            
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
        POST_SAREK(vcf_ch, bam_ch, bed_ch)
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
        POST_SAREK(vcf_ch, bam_ch, bed_ch)
}

workflow RUN_FULL_VARIANT_CALLING {
    take:
        bed_ch

    main:
        log.info """
        ╔════════════════════════════════════════════════════════════╗
        ║  Running full variant calling pipeline (Sarek)             ║
        ║  Consensus calling: ${params.create_consensus ? 'ENABLED' : 'DISABLED'}
        ╚════════════════════════════════════════════════════════════╝
        """.stripIndent()

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

        if (params.create_consensus) {
            ref_fasta_ch = Channel.value(file(params.ref_fasta))
            ref_fai_ch   = Channel.value(file(params.ref_fasta + ".fai"))

            CONSENSUS_CALLING(
                COLLECT_VARIANT_CALLING_OUTPUTS.out.dv_vcf,
                COLLECT_VARIANT_CALLING_OUTPUTS.out.hc_vcf,
                ref_fasta_ch,
                ref_fai_ch
            )

            final_vcf_ch = CONSENSUS_CALLING.out.consensus_vcf
        } else {
            final_vcf_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.dv_vcf
        }

        def vcf_with_meta_ch = final_vcf_ch.map { sample, vcf ->
            def meta = [ sample: sample, assay: assayMap.get(sample, 'NA') ]
            tuple(meta, vcf)
        }

        def bam_with_meta_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.bam.map { sample, bam, bai ->
            def meta = [ sample: sample, assay: assayMap.get(sample, 'NA') ]
            tuple(meta, bam, bai)
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

        POST_SAREK(vcf_with_meta_ch, bam_with_meta_ch, bed_ch)
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow {
    
    // Validate and load BED file
    def bedFile = validateBedFile()
    bed_ch = Channel.value(bedFile)
    RUN_FULL_VARIANT_CALLING(bed_ch)
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
