// subworkflows/db_ingest.nf
// ───────────────────────────────────────────────────────────────────────────
// DB_INGEST — reads the run-output manifest produced by MANIFEST, filters
// samples whose qc_recommendation == "READY", then runs the db_export
// pipeline for each passing sample:
//   variant extraction (no VEP) → mosdepth coverage → coverage load → variant load
//
// The consensus VCF is used directly without re-annotation. Annotation
// columns in the variants TSV (gene, hgvsc, impact, clinvar) are left
// empty and populated downstream by the database annotation layer.
//
// Inputs:
//   manifest_file   — path to run_output_manifest.tsv
//                     columns: run_id, sample_id, assay, final_vcf,
//                              coverage_bam, qc_status, qc_recommendation
//
// The BAM index is derived as <coverage_bam>.bai
// ───────────────────────────────────────────────────────────────────────────
nextflow.enable.dsl = 2

include { EXTRACT_VARIANTS    } from '../external/db_export/modules/extract_variants'
include { MOSDEPTH_THRESHOLDS } from '../external/db_export/modules/mosdepth_thresholds'
include { CONVERT_THRESHOLDS  } from '../external/db_export/modules/thresholds_to_coverage'
include { LOAD_VARIANTS       } from '../modules/load_variants'
include { LOAD_DENSE_DEPTH    } from '../modules/load_dense_depth'

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
        ch_bins_bed = Channel.value(file(params.bins_bed, checkIfExists: true))

        // ── Step 1: Extract variants directly from the consensus VCF ──────
        // Uses a local script that does not require VEP/CSQ annotation.
        ch_vcf_input       = ch_pass.map { meta, vcf, bam, bai -> [ meta, vcf ] }
        ch_variants_script = Channel.value(
            file("${workflow.projectDir}/scripts/db_vcf_to_variants.py")
        )
        EXTRACT_VARIANTS(ch_vcf_input, ch_variants_script)

        // ── Step 2: BAM coverage ──────────────────────────────────────────
        ch_bam_input = ch_pass.map { meta, vcf, bam, bai -> [ meta, bam, bai ] }
        MOSDEPTH_THRESHOLDS(ch_bam_input, ch_bins_bed)

        // ── Step 3: Convert mosdepth thresholds to coverage TSV ──────────
        ch_coverage_script = Channel.value(
            file("${workflow.projectDir}/external/db_export/bins/db_depth_to_coverage_new.py")
        )
        CONVERT_THRESHOLDS(MOSDEPTH_THRESHOLDS.out.per_base_bed, ch_coverage_script)

        // ── Step 4: Load dense depth ──────────────────────────────────────
        LOAD_DENSE_DEPTH(CONVERT_THRESHOLDS.out.coverage_tsv)

        // ── Step 5: Load variants (after coverage succeeds) ───────────────
        ch_variants_ready = EXTRACT_VARIANTS.out.variants_tsv
            .join(LOAD_DENSE_DEPTH.out.loaded_sample.map { meta -> [ meta, true ] })
            .map { meta, variants_tsv, _done -> [ meta, variants_tsv ] }

        LOAD_VARIANTS(ch_variants_ready)

    emit:
        loaded = LOAD_VARIANTS.out.loaded_sample
}
