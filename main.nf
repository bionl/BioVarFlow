nextflow.enable.dsl=2

// ═══════════════════════════════════════════════════════════════════════════
// IMPORTS
// ═══════════════════════════════════════════════════════════════════════════

include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION } \
  from './external/sarek/subworkflows/local/utils_nfcore_sarek_pipeline'

include { NFCORE_SAREK } \
  from './external/sarek/main.nf'

include { POST_SAREK; POST_SAREK_SOMATIC } \
  from './modules/vep.nf'

include { CONSENSUS_CALLING } \
  from './modules/consensus.nf'

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
params.somatic_bed            = params.somatic_bed ?: "${workflow.projectDir}/data/Somatic_125genes_MANE_50bp.bed"
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

// Somatic reporting panel BED — only required in somatic mode.
def validateSomaticBedFile() {
    def bedFile = params.somatic_bed ? file(params.somatic_bed) : null
    if (!bedFile?.exists()) {
        error "❌ Somatic panel BED file not found: ${params.somatic_bed}"
    }
    return bedFile
}

// Reads the somatic samplesheet and injects a duplicate status=0 row (with a
// _germline patient/sample suffix) for every tumor-only sample — i.e. any
// patient that has a status=1 row but no matching status=0 row.  The two rows
// share the same BAM/BAI path so Sarek runs HaplotypeCaller on the tumor BAM
// as if it were a germline sample, without affecting Mutect2 pairing (Sarek
// pairs by patient ID, and the injected row carries a distinct patient ID).
def injectGermlineRows(String inputPath) {
    def inputFile  = file(inputPath)
    def lines      = inputFile.readLines()
    def header     = lines[0]
    def colNames   = header.split(',')*.trim()

    def rows = lines.tail().findAll { it.trim() }.collect { line ->
        def vals = line.split(',', -1)*.trim()
        [colNames, vals].transpose().collectEntries { k, v -> [(k): v] }
    }

    // Identify patients that have only tumor rows (status=1, no status=0)
    def patientStatuses = [:]
    rows.each { row ->
        def p = row['patient']
        if (!patientStatuses.containsKey(p)) patientStatuses[p] = [] as Set
        patientStatuses[p] << (row['status'] as Integer)
    }
    def tumorOnlyPatients = patientStatuses
        .findAll { p, statuses -> 1 in statuses && !(0 in statuses) }
        .keySet()

    if (!tumorOnlyPatients) {
        log.info "✓ No tumor-only samples found — samplesheet unchanged"
        return inputPath
    }

    log.info "✓ Injecting germline rows for tumor-only patient(s): ${tumorOnlyPatients.join(', ')}"

    def injected = []
    rows.each { row ->
        if (row['patient'] in tumorOnlyPatients && (row['status'] as Integer) == 1) {
            def fake = new LinkedHashMap(row)
            fake['patient'] = row['patient'] + '_germline'
            fake['sample']  = row['sample']  + '_germline'
            fake['status']  = '0'
            injected << fake
        }
    }

    def modifiedFile = file("${workflow.workDir}/samplesheet_somatic_with_germline.csv")
    modifiedFile.parent.mkdirs()
    modifiedFile.text = header + '\n' +
        (rows + injected)
            .collect { row -> colNames.collect { row[it] ?: '' }.join(',') }
            .join('\n') + '\n'

    log.info "✓ Modified samplesheet written to: ${modifiedFile}"
    return modifiedFile.toString()
}

// Reads the somatic samplesheet and builds a map from normal sample name →
// tumor sample name for every paired patient (a patient that has both a
// status=0 and a status=1 row).  Used to remap HC VCFs produced for the real
// normal into the tumor sample namespace so they can be joined with the
// somatic VEP channel (which is keyed by tumor name) for combined reporting.
//
// Example: PatientA / NormalA (status=0) + PatientA / TumorA (status=1)
//          → [ "NormalA": "TumorA" ]
def buildNormalToTumorMap(String inputPath) {
    def inputFile = file(inputPath)
    def lines     = inputFile.readLines()
    def header    = lines[0]
    def colNames  = header.split(',')*.trim()

    def rows = lines.tail().findAll { it.trim() }.collect { line ->
        def vals = line.split(',', -1)*.trim()
        [colNames, vals].transpose().collectEntries { k, v -> [(k): v] }
    }

    // Group by patient → collect normal and tumor sample names
    def byPatient = [:]
    rows.each { row ->
        def p = row['patient']
        if (!byPatient.containsKey(p)) byPatient[p] = [ normals: [], tumors: [] ]
        if ((row['status'] as Integer) == 0) byPatient[p].normals << row['sample']
        if ((row['status'] as Integer) == 1) byPatient[p].tumors  << row['sample']
    }

    // Only keep patients that have BOTH a normal and a tumor (paired case)
    def normalToTumor = [:]
    byPatient.each { patient, samples ->
        if (samples.normals && samples.tumors) {
            // For each normal, map it to its paired tumor(s).
            // In the common 1-normal : 1-tumor case this is a simple 1:1 map.
            samples.normals.each { normal ->
                // Map to the first tumor; extend if multi-tumor per normal is needed.
                normalToTumor[normal] = samples.tumors[0]
            }
            log.info "✓ Paired patient ${patient}: normal(s) ${samples.normals} → tumor(s) ${samples.tumors}"
        }
    }

    if (!normalToTumor) {
        log.info "✓ No paired tumor-normal patients found — normalToTumorMap is empty"
    }

    return normalToTumor
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

        // Run post-processing (no somatic VCF in this path)
        POST_SAREK(vcf_ch, bam_ch, bed_ch, Channel.empty())
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
        POST_SAREK(vcf_ch, bam_ch, bed_ch, Channel.empty())
}

workflow RUN_FULL_VARIANT_CALLING {
    take:
        bed_ch

    main:
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

        // In somatic mode, inject fake status=0 rows for tumor-only samples so
        // Sarek runs HaplotypeCaller on them for ACMG SF reporting.
        def samplesheetPath = (params.somatic_mode && params.input)
            ? injectGermlineRows(params.input.toString())
            : params.input

        PIPELINE_INITIALISATION(
            params.version,
            params.validate_params,
            args,
            params.outdir,
            samplesheetPath,
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

            // Build normal→tumor map early — needed to pick the right Mutect2 VCF.
            // When a patient has a paired run, Sarek emits TWO mutect2 dirs:
            //   mutect2/HCC1395T/         ← tumor-only run (wrong for paired patients)
            //   mutect2/HCC1395T_vs_HCC1395N/ ← paired run (correct for paired patients)
            // Stripping _vs_* from both yields the same key "HCC1395T", causing a
            // duplicate channel entry. We must pick only the paired dir for paired
            // patients and only the non-_vs_ dir for tumor-only patients.
            def samplesheetForMap0 = params.input?.toString() ?: ''
            def normalToTumorMap0  = samplesheetForMap0 ? buildNormalToTumorMap(samplesheetForMap0) : [:]
            // Collect all tumor sample names that have a paired normal
            def pairedTumors = normalToTumorMap0.values().toSet()

            def mutect2_vcf_ch = NFCORE_SAREK.out.multiqc_report
                .flatMap {
                    file("${params.outdir}/variant_calling/mutect2/*/*.mutect2.filtered.vcf.gz", checkIfExists: !isGCS)
                }
                .filter { vcf -> vcf.name.endsWith('.vcf.gz') && !vcf.name.endsWith('.tbi') }
                .filter { vcf ->
                    def dirName = vcf.parent.name
                    def isPaired = dirName.contains('_vs_')
                    // For paired patients: only keep the _vs_ dir (paired Mutect2 output).
                    // For tumor-only patients: only keep the non-_vs_ dir.
                    if (isPaired) {
                        def tumorName = dirName.replaceAll(/_vs_.+$/, '')
                        return pairedTumors.contains(tumorName)
                    } else {
                        return !pairedTumors.contains(dirName)
                    }
                }
                .map { vcf ->
                    // Sarek names paired dirs as TumorA_vs_NormalA — strip _vs_* to get
                    // the tumor sample name so publishDir and downstream joins use it.
                    def sampleName = vcf.parent.name.replaceAll(/_vs_.+$/, '')
                    tuple(sampleName, vcf)
                }

            def rescue_script_ch = Channel.value(file("${params.scriptdir}/mutect2_rescue.py"))
            MUTECT2_RESCUE(mutect2_vcf_ch, rescue_script_ch)
            final_vcf_ch = MUTECT2_RESCUE.out.map { sample, vcf, tbi -> tuple(sample, vcf) }
        } else if (params.create_consensus) {
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

        // POST_SAREK for germline annotation:
        // — germline mode:      consensus VCF as usual
        // — somatic tumor-only: HC VCF from injected _germline rows, suffix stripped
        // — somatic paired:     HC VCF from real normal, remapped to tumor name via normalToTumorMap
        if (params.somatic_mode) {
            def isGCS = isGcsPath(params.outdir)

            // Build normal→tumor map for paired patients (empty map if all tumor-only)
            def samplesheetForMap = params.input?.toString() ?: ''
            def normalToTumorMap  = samplesheetForMap ? buildNormalToTumorMap(samplesheetForMap) : [:]

            // ── Somatic VEP annotation (rescued Mutect2 VCF → somatic.vep.vcf) ──
            // The full in-panel set is annotated once; PASS and the tumor-only
            // post-hoc thresholds are applied in the report script, which infers
            // paired vs tumor-only from the VCF's own ##tumor_sample header.
            def somatic_bed_ch = Channel.value(validateSomaticBedFile())
            POST_SAREK_SOMATIC(vcf_with_meta_ch, somatic_bed_ch)

            def somatic_vep_ch = POST_SAREK_SOMATIC.out.somatic_vep

            // ── HC VCF collection — all filtered VCFs from haplotypecaller output ──
            def all_hc_vcf_ch = NFCORE_SAREK.out.multiqc_report
                .flatMap {
                    file("${params.outdir}/variant_calling/haplotypecaller/*/*filtered.vcf.gz", checkIfExists: !isGCS)
                }
                .filter { vcf ->
                    vcf.name.endsWith('.vcf.gz') &&
                    !vcf.name.contains('.g.vcf.gz') &&
                    !vcf.name.endsWith('.tbi')
                }

            // ── Tumor-only: dirs ending in _germline → strip suffix to get tumor name ──
            def tumor_only_vcf_ch = all_hc_vcf_ch
                .filter { vcf -> vcf.parent.name.endsWith('_germline') }
                .map { vcf ->
                    def origSample = vcf.parent.name.replaceAll(/_germline$/, '')
                    def meta = [ sample: origSample, assay: assayMap.get(origSample, 'NA') ]
                    tuple(meta, vcf)
                }

            // ── Paired: real normal dirs → remap NormalA → TumorA via normalToTumorMap ──
            def paired_vcf_ch = all_hc_vcf_ch
                .filter { vcf ->
                    !vcf.parent.name.endsWith('_germline') &&
                    normalToTumorMap.containsKey(vcf.parent.name)
                }
                .map { vcf ->
                    def tumorSample = normalToTumorMap[vcf.parent.name]
                    def meta = [ sample: tumorSample, assay: assayMap.get(tumorSample, 'NA') ]
                    tuple(meta, vcf)
                }

            // Merge both branches — each item is already keyed by tumor sample name
            def hc_germline_vcf_ch = tumor_only_vcf_ch.mix(paired_vcf_ch)

            // ── BAM collection — same two-branch logic ──
            def all_bam_ch = COLLECT_VARIANT_CALLING_OUTPUTS.out.bam

            def tumor_only_bam_ch = all_bam_ch
                .filter { sample, bam, bai -> sample.endsWith('_germline') }
                .map { sample, bam, bai ->
                    def origSample = sample.replaceAll(/_germline$/, '')
                    def meta = [ sample: origSample, assay: assayMap.get(origSample, 'NA') ]
                    tuple(meta, bam, bai)
                }

            def paired_bam_ch = all_bam_ch
                .filter { sample, bam, bai ->
                    !sample.endsWith('_germline') &&
                    normalToTumorMap.containsKey(sample)
                }
                .map { sample, bam, bai ->
                    def tumorSample = normalToTumorMap[sample]
                    def meta = [ sample: tumorSample, assay: assayMap.get(tumorSample, 'NA') ]
                    tuple(meta, bam, bai)
                }

            def hc_germline_bam_ch = tumor_only_bam_ch.mix(paired_bam_ch)

            hc_germline_vcf_ch.ifEmpty {
                log.warn "⚠️  No HC germline VCFs found — ACMG SF reporting skipped for all somatic samples."
            }

            // POST_SAREK: HC germline VCF (ACMG SF) + somatic VEP — both keyed by tumor name
            POST_SAREK(hc_germline_vcf_ch, hc_germline_bam_ch, bed_ch, somatic_vep_ch)
        } else {
            POST_SAREK(vcf_with_meta_ch, bam_with_meta_ch, bed_ch, Channel.empty())
        }
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
