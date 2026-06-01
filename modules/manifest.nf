// modules/manifest.nf
// ───────────────────────────────────────────────────────────────────────────
// Run-Output Manifest — final aggregation step.
//
// Produces a single TSV (one row per sample) at:
//     ${params.outdir}/manifest/run_output_manifest.tsv
//
// Columns:
//     run_id, sample_id, assay, final_vcf, coverage_bam, qc_status, qc_recommendation
//
// The manifest is consumed by the DB_INGEST subworkflow (subworkflows/db_ingest.nf),
// which filters READY samples and triggers bionl_lean_db_export for each.
// Sample filtering (qc_recommendation == "READY") is done in DB_INGEST.
//
// All paths in the manifest are absolute, *published* paths (i.e. derived
// from ${params.outdir}, not from the temporary work directory). Both
// local paths and gs:// URIs are supported transparently.
// ───────────────────────────────────────────────────────────────────────────
nextflow.enable.dsl=2

params.outdir   = params.outdir   ?: "results"
params.scriptdir = params.scriptdir ?: "${workflow.projectDir}/scripts"

/********************  PROCESS  ********************/

process CREATE_OUTPUT_MANIFEST {
    tag "${run_id}"
    publishDir "${params.outdir}/manifest", mode: 'copy', overwrite: true

    input:
        val  run_id                       // String — run identifier
        val  sample_info                  // List<String> "sample\tvcf\tbam"
        path qc_jsons, stageAs: 'qc/*'    // List<Path> per-sample QC JSON files
        path script                       // build_run_output_manifest.py

    output:
        path "run_output_manifest.tsv"

    script:
    """
    set -euo pipefail

    cat > sample_info.tsv <<'__EOF__'
${sample_info.join('\n')}
__EOF__

    python3 ${script} \\
        --run-id '${run_id}' \\
        --sample-info sample_info.tsv \\
        --qc-dir qc \\
        --out run_output_manifest.tsv
    """
}

/********************  WORKFLOW  ********************/

workflow MANIFEST {
    take:
        final_vcf_ch    // tuple(meta,  vcf)        — final consensus VCF
        bam_ch          // tuple(meta,  bam, bai)   — coverage/alignment BAM
        qc_json_ch      // tuple(meta,  json, tsv)  — QC Gate per-sample JSON

    main:
        // Resolve the run identifier. Prefer params.run_id when supplied,
        // otherwise fall back to Nextflow's auto-generated runName so that
        // the manifest is always emitted.
        def resolved_run_id = (params.containsKey('run_id') && params.run_id) \
            ? params.run_id.toString() \
            : workflow.runName.toString()

        // Build canonical, *published* per-sample paths under params.outdir.
        // We deliberately construct these from params.outdir (rather than
        // using the work-dir path of the staged file) so the manifest always
        // points at the published artefact (gs:// URI when running on GCS).
        def vcf_path_ch = final_vcf_ch.map { meta, vcf ->
            def sample = meta.sample
            def assay  = meta.assay ?: 'NA'
            def published = "${params.outdir}/${sample}/consensus/${vcf.name}"
            tuple(sample, assay, published)
        }

        def bam_path_ch = bam_ch.map { meta, bam, bai ->
            def sample = meta.sample
            def published = "${params.outdir}/preprocessing/mapped/${sample}/${bam.name}"
            tuple(sample, published)
        }

        def qc_json_path_ch = qc_json_ch.map { meta, json, tsv ->
            tuple(meta.sample, json)
        }

        // Inner-join by sample so we only emit the manifest once every
        // required artefact (final VCF, BAM, QC JSON) is available for the
        // sample. Missing artefacts simply skip the sample.
        def joined_ch = qc_json_path_ch
            .join(vcf_path_ch)            // (sample, json, assay, vcf_path)
            .join(bam_path_ch)            // (sample, json, assay, vcf_path, bam_path)

        // Split into two parallel branches that we collect independently:
        //   info: a String row per sample (sample_id, assay, vcf, bam)
        //   json: the staged QC JSON file for that sample
        joined_ch.multiMap { sample, json, assay, vcf_path, bam_path ->
            info: "${sample}\t${assay}\t${vcf_path}\t${bam_path}"
            json: json
        }.set { split_ch }

        def script_ch = Channel
            .fromPath("${params.scriptdir}/build_run_output_manifest.py")
            .first()

        CREATE_OUTPUT_MANIFEST(
            Channel.value(resolved_run_id),
            split_ch.info.collect(),
            split_ch.json.collect(),
            script_ch
        )

    emit:
        manifest = CREATE_OUTPUT_MANIFEST.out
}
