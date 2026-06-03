// subworkflows/db_ingest.nf
// ───────────────────────────────────────────────────────────────────────────
// DB_INGEST — reads the run-output manifest produced by MANIFEST, filters
// samples whose qc_recommendation == "READY", then runs the full db_export
// pipeline (VEP annotation → variant extraction → mosdepth coverage →
// coverage load → variant load) for each passing sample.
//
// Inputs:
//   manifest_file   — path to run_output_manifest.tsv
//                     columns: run_id, sample_id, assay, final_vcf,
//                              coverage_bam, qc_status, qc_recommendation
//
// The BAM index is derived as <coverage_bam>.bai — Sarek always produces
// <sample>.md.cram.bai alongside the BAM, but the published BAM is named
// *.sorted.bam and its index is *.sorted.bam.bai.
// ───────────────────────────────────────────────────────────────────────────
nextflow.enable.dsl = 2

include { VEP_Annotate_DB } from '../external/db_export/modules/vep_annotate'
include { EXTRACT_VARIANTS    } from '../external/db_export/modules/extract_variants'
include { MOSDEPTH_THRESHOLDS } from '../external/db_export/modules/mosdepth_thresholds'
include { CONVERT_THRESHOLDS  } from '../external/db_export/modules/thresholds_to_coverage'
include { LOAD_VARIANTS       } from '../external/db_export/modules/load_variants'
include { LOAD_DENSE_DEPTH    } from '../external/db_export/modules/load_dense_depth'

workflow DB_INGEST {
    take:
        manifest_file   // path — run_output_manifest.tsv

    main:
        // ── Parse manifest and filter READY samples ───────────────────────
        ch_pass = manifest_file
            .splitCsv(header: true, sep: '\t')
            .filter { row ->
                def rec = (row.qc_recommendation ?: '').trim()
                if (rec != 'READY') {
                    log.info "[DB_INGEST] Skipping ${row.sample_id} (qc_recommendation=${rec ?: 'missing'})"
                }
                rec == 'READY'
            }
            .map { row ->
                def meta      = [ sample: row.sample_id, assay: (row.assay ?: 'NA') ]
                def vcf       = file(row.final_vcf,    checkIfExists: false)
                def bam       = file(row.coverage_bam, checkIfExists: false)
                def bam_index = file("${row.coverage_bam}.bai", checkIfExists: false)
                [ meta, vcf, bam, bam_index ]
            }

        ch_pass.count().subscribe { n ->
            if (n == 0) {
                log.warn "[DB_INGEST] No samples passed QC (qc_recommendation == READY). Skipping DB ingestion."
            } else {
                log.info "[DB_INGEST] ${n} sample(s) passed QC — starting DB ingestion."
            }
        }

        // ── Shared reference channels ─────────────────────────────────────
        ch_bins_bed      = Channel.value(file(params.bins_bed,    checkIfExists: true))
        ch_vep_cache     = Channel.value(file(params.vep_cache,   checkIfExists: true))
        ch_vep_fasta     = Channel.value(file(params.vep_fasta,   checkIfExists: true))
        ch_vep_fasta_fai = Channel.value(file(params.vep_fasta + ".fai", checkIfExists: true))
        ch_vep_plugins   = Channel.value(file(params.vep_plugins))
        ch_clinvar_vcf   = Channel.value(file(params.clinvar_vcf,     checkIfExists: true))
        ch_clinvar_vcf_tbi = Channel.value(file(params.clinvar_vcf + ".tbi", checkIfExists: true))

        // ── Step 1: VEP annotation ────────────────────────────────────────
        ch_vcf_only = ch_pass.map { meta, vcf, bam, bai -> [ meta, vcf ] }

        VEP_Annotate_DB(
            ch_vcf_only,
            ch_vep_cache,
            ch_vep_fasta,
            ch_vep_fasta_fai,
            ch_clinvar_vcf,
            ch_clinvar_vcf_tbi,
            ch_vep_plugins
        )

        // ── Step 2: Extract variants from VEP-annotated VCF ──────────────
        ch_variants_script = Channel.value(
            file("${workflow.projectDir}/external/db_export/bins/db_vep_vcf_to_variants_all.py")
        )
        EXTRACT_VARIANTS(VEP_Annotate_DB.out.vep_vcf, ch_variants_script)

        // ── Step 3: BAM coverage ──────────────────────────────────────────
        ch_bam_input = ch_pass.map { meta, vcf, bam, bai -> [ meta, bam, bai ] }
        MOSDEPTH_THRESHOLDS(ch_bam_input, ch_bins_bed)

        // ── Step 4: Convert mosdepth thresholds to coverage TSV ──────────
        ch_coverage_script = Channel.value(
            file("${workflow.projectDir}/external/db_export/bins/db_depth_to_coverage_new.py")
        )
        CONVERT_THRESHOLDS(MOSDEPTH_THRESHOLDS.out.per_base_bed, ch_coverage_script)

        // ── Step 5: Load dense depth ──────────────────────────────────────
        LOAD_DENSE_DEPTH(CONVERT_THRESHOLDS.out.coverage_tsv)

        // ── Step 6: Load variants (after coverage succeeds) ───────────────
        ch_variants_ready = EXTRACT_VARIANTS.out.variants_tsv
            .join(LOAD_DENSE_DEPTH.out.loaded_sample.map { meta -> [ meta, true ] })
            .map { meta, variants_tsv, _done -> [ meta, variants_tsv ] }

        LOAD_VARIANTS(ch_variants_ready)

    emit:
        loaded = LOAD_VARIANTS.out.loaded_sample
}
