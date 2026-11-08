// ═══════════════════════════════════════════════════════════════════════════
// CNV post-processing
//
// Sarek's `cnvkit batch` gives us per-sample coverage only. Everything that
// makes the calls usable happens here:
//
//   CNVKIT_POSTPROCESS  per sample: check the PoN fits, re-fix against it,
//                       re-segment with sensitive CBS settings, call with the
//                       samplesheet sex, and sanity-check chrX
//   CNV_FILTER          run level: tiered filter + annotation over every
//                       sample's calls, producing a TSV per sample and one
//                       Excel workbook
//
// Why not sarek's own .cnr/.cns (validated in docs/cnv_validation.md):
//   * with --intervals the recalibrated CRAMs keep only on-target reads, so
//     every antitarget bin has zero coverage and sarek's .cnr carries ~12k
//     fake "neutral" bins. `fix` against a real PoN drops them.
//   * `cnvkit batch` does not expose -m cbs -t 0.001 --smooth-cbs
//     --drop-low-coverage, and its defaults miss small deletions and merge
//     mosaic gains into whole arms.
// ═══════════════════════════════════════════════════════════════════════════

params.scriptdir = params.scriptdir ?: "${workflow.projectDir}/scripts"

// Sarek's samplesheet carries XX / XY / NA; CNVkit wants female / male.
// Sex is not optional here: with empty antitargets CNVkit cannot guess it
// reliably, and a wrong guess puts a normal male chrX at cn = 0.
def cnvkitSex(meta) {
    switch (meta.sex?.toString()?.toUpperCase()) {
        case 'XX': return 'female'
        case 'XY': return 'male'
        default:
            error "❌ CNV calling needs the sex of ${meta.sample} in the samplesheet " +
                  "(got '${meta.sex}'). Set it to XX or XY, or drop cnvkit from --tools."
    }
}

process CNVKIT_POSTPROCESS {
    tag "${meta.sample}"
    label 'process_medium'

    publishDir "${params.outdir}/cnv/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(target_cnn), path(antitarget_cnn)
    path pon
    path pon_check_script
    path chrx_check_script

    output:
    tuple val(meta), path("${meta.sample}.call.cns"),    emit: calls
    tuple val(meta), path("${meta.sample}.cnr"),         emit: cnr
    tuple val(meta), path("${meta.sample}.metrics.tsv"), emit: metrics
    tuple val(meta), path("${meta.sample}.cnv_qc.tsv"),  emit: qc
    path "${meta.sample}.scatter.png", optional: true
    path "versions.yml", emit: versions

    script:
    def sample   = meta.sample
    def sex      = cnvkitSex(meta)
    // Only pass -y when the PoN itself was built with it. A -y against a
    // female-normalised reference made a normal male chrX read as cn = 0.
    def male_ref = params.cnv_pon_haploid_x ? '-y' : ''
    def haploid  = params.cnv_pon_haploid_x ? '--haploid-x-reference' : ''
    """
    # Stop now if the PoN describes different bins than this sample has.
    python3 ${pon_check_script} ${target_cnn} ${pon}

    # CNVkit's own guess: recorded for the QC table, not used for calling.
    cnvkit.py sex ${target_cnn} > ${sample}.sexcheck.tsv

    cnvkit.py fix ${target_cnn} ${antitarget_cnn} ${pon} -o ${sample}.cnr

    cnvkit.py segment ${sample}.cnr \\
        -m cbs -t ${params.cnv_segment_threshold} --smooth-cbs --drop-low-coverage \\
        -o ${sample}.cns

    cnvkit.py call ${sample}.cns \\
        --sample-sex ${sex} ${male_ref} --center median \\
        -o ${sample}.call.cns

    cnvkit.py metrics ${sample}.cnr -s ${sample}.cns -o ${sample}.metrics.tsv

    # Cosmetic; a plotting failure must not cost the run its calls.
    cnvkit.py scatter ${sample}.cnr -s ${sample}.cns -o ${sample}.scatter.png || \\
        echo "WARNING: scatter plot failed for ${sample}"

    python3 ${chrx_check_script} \\
        --cnr ${sample}.cnr \\
        --sample ${sample} \\
        --sex ${sex} \\
        --sexcheck ${sample}.sexcheck.tsv \\
        --tolerance ${params.cnv_chrx_tolerance} ${haploid} \\
        --out ${sample}.cnv_qc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version)
    END_VERSIONS
    """
}

process CNV_MANIFEST_CHECK {
    label 'process_low'

    publishDir "${params.outdir}/cnv", mode: 'copy'

    input:
    path pon
    path targets_bed
    path common_bed
    path manifest
    path check_script

    output:
    path "cnv_pon_provenance.tsv", emit: provenance

    script:
    // Optional inputs arrive as the NO_FILE placeholder. A missing manifest is a
    // warning inside the script, not a reason to skip the process: the run should
    // still record in the report that the PoN's provenance was never verified.
    def manifest_arg = manifest.name   == 'NO_FILE' ? '' : "--manifest ${manifest}"
    def common_arg   = common_bed.name == 'NO_FILE' ? '' : "--common-bed ${common_bed}"
    def haploid_arg  = params.cnv_pon_haploid_x ? '--haploid-x-reference' : ''
    """
    python3 ${check_script} \\
        --pon ${pon} \\
        --targets-bed ${targets_bed} \\
        ${manifest_arg} \\
        ${common_arg} \\
        ${haploid_arg} \\
        --cnvkit-version ${params.cnv_cnvkit_version} \\
        --min-pon-size ${params.cnv_min_pon_size} \\
        --out cnv_pon_provenance.tsv
    """
}

process CNV_FILTER {
    label 'process_low'

    publishDir "${params.outdir}/cnv", mode: 'copy'

    input:
    path calls, stageAs: 'calls/*'       // every sample's *.call.cns
    path metrics, stageAs: 'metrics/*'   // every sample's *.metrics.tsv
    path cnrs, stageAs: 'cnr/*'          // every sample's *.cnr, for the hom-del scan
    path sex_tsv
    path filter_script
    path pon
    path panel_bed
    path gnomad_sv
    path clingen
    path segdup
    path gnomad_exome_cnv
    path common_bed
    path provenance      // cnv_pon_provenance.tsv, for the report's Thresholds sheet

    output:
    path "cnv_filtered_all.tsv", emit: all_calls
    path "*.cnv_filtered.tsv",   emit: per_sample
    path "cnv_report.xlsx",      emit: xlsx

    script:
    // Optional resources arrive as the NO_FILE placeholder when unset.
    def exome_cnv_arg = gnomad_exome_cnv.name == 'NO_FILE' ? '' : "--gnomad-exome-cnv ${gnomad_exome_cnv}"
    def common_arg    = common_bed.name       == 'NO_FILE' ? '' : "--common-bed ${common_bed}"
    def pon_size_arg  = params.cnv_pon_size   ? "--pon-size ${params.cnv_pon_size}" : ''
    def haploid_arg   = params.cnv_pon_haploid_x ? '--haploid-x-reference' : ''
    """
    # `cnvkit metrics` names its single row by .cnr path and the filter keys on
    # that basename, so one concatenated table covering every sample is enough.
    head -1 \$(ls metrics/*.metrics.tsv | head -1) > metrics_all.tsv
    for m in metrics/*.metrics.tsv; do tail -n +2 \$m >> metrics_all.tsv; done

    python3 ${filter_script} \\
        --calls calls/*.call.cns \\
        --cnr cnr/*.cnr \\
        --sex ${sex_tsv} \\
        --pon ${pon} ${pon_size_arg} ${haploid_arg} \\
        --panel-bed ${panel_bed} \\
        --panel-name '${params.cnv_panel_name}' \\
        --gnomad-sv ${gnomad_sv} \\
        --clingen ${clingen} \\
        --segdup ${segdup} \\
        ${exome_cnv_arg} \\
        ${common_arg} \\
        --metrics metrics_all.tsv \\
        --max-bivar ${params.cnv_max_bivar} \\
        --min-batch-refs ${params.cnv_min_pon_size} \\
        --min-probes ${params.cnv_min_exons} \\
        --del-log2 ${params.cnv_del_log2} \\
        --dup-log2 ${params.cnv_dup_log2} \\
        --mosaic-probes ${params.cnv_mosaic_exons} \\
        --mosaic-log2 ${params.cnv_mosaic_log2} \\
        --min-depth ${params.cnv_min_depth} \\
        --max-spread ${params.cnv_max_spread} \\
        --max-af ${params.cnv_max_af} \\
        --recip ${params.cnv_recip_overlap} \\
        --homdel-log2 ${params.cnv_homdel_log2} \\
        --homdel-min-bins ${params.cnv_homdel_min_bins} \\
        --homdel-cohort-frac ${params.cnv_homdel_cohort_frac} \\
        --homdel-min-cohort ${params.cnv_homdel_min_cohort} \\
        --provenance ${provenance} \\
        --small-exons ${params.cnv_small_exons} \\
        --outdir . \\
        --out-xlsx cnv_report.xlsx
    """
}

process CNV_QC_SUMMARY {
    label 'process_low'

    publishDir "${params.outdir}/cnv", mode: 'copy'

    input:
    path qc_rows, stageAs: 'qc/*'

    output:
    path "cnv_qc_summary.tsv", emit: summary

    script:
    """
    # One table so a flagged sex or chrX status is visible without opening
    # every per-sample folder.
    head -1 \$(ls qc/*.cnv_qc.tsv | head -1) > cnv_qc_summary.tsv
    for q in qc/*.cnv_qc.tsv; do tail -n +2 \$q >> cnv_qc_summary.tsv; done

    if awk -F'\\t' 'NR>1 && (\$8=="FAIL" || \$4=="no")' cnv_qc_summary.tsv | grep -q .; then
        echo "WARNING: some samples failed the chrX check or disagree with cnvkit's sex guess:"
        awk -F'\\t' 'NR>1 && (\$8=="FAIL" || \$4=="no") {print "  "\$1": sex="\$2" guess="\$3" chrX="\$5" ("\$8")"}' cnv_qc_summary.tsv
    fi
    """
}

// ═══════════════════════════════════════════════════════════════════════════
// CNV_CALLING — wires the processes together
// ═══════════════════════════════════════════════════════════════════════════

workflow CNV_CALLING {
    take:
        coverage_ch     // tuple(meta, targetcoverage.cnn, antitargetcoverage.cnn)
        sex_tsv_ch      // sample<TAB>male|female for every sample in the run
        panel_bed_ch    // exon BED of reportable genes, for the relevance level

    main:
        if (!params.cnv_pon) {
            error "❌ cnvkit is in --tools but --cnv_pon is not set. This pipeline does not " +
                  "build panels of normals; supply the customer's CNVkit reference .cnn, " +
                  "built from the same BED as --cnv_targets_bed."
        }

        // The relevance panel arrives as a channel from main.nf, where the
        // fallback to params.bed is resolved: params assigned at a script's top
        // level are not visible inside an included module.
        //
        // These have no sensible default, and a null reaches file() as an
        // "Argument of `file()` function cannot be null" that names nothing.
        [ cnv_gnomad_sv: params.cnv_gnomad_sv,
          cnv_clingen  : params.cnv_clingen,
          cnv_segdup   : params.cnv_segdup ].each { name, value ->
            if (!value) {
                error "❌ cnvkit is in --tools but --${name} is not set. " +
                      "See docs/cnv.md for the resources CNV filtering needs."
            }
        }

        pon_ch        = Channel.value(file(params.cnv_pon, checkIfExists: true))
        pon_check_ch  = Channel.value(file("${params.scriptdir}/cnv_pon_check.py", checkIfExists: true))
        chrx_check_ch = Channel.value(file("${params.scriptdir}/cnv_chrx_check.py", checkIfExists: true))
        filter_ch     = Channel.value(file("${params.scriptdir}/filter_cnv.py", checkIfExists: true))

        CNVKIT_POSTPROCESS(coverage_ch, pon_ch, pon_check_ch, chrx_check_ch)

        def noFile = file("${workflow.projectDir}/assets/NO_FILE", checkIfExists: true)

        // manifest.json ships beside the .cnn, so default to that unless told
        // otherwise. file().parent keeps this working for gs:// paths too.
        def manifestPath = params.cnv_pon_manifest ?: "${file(params.cnv_pon).parent}/manifest.json"
        def manifestFile = file(manifestPath)
        if (params.cnv_pon_manifest && !manifestFile.exists()) {
            error "❌ --cnv_pon_manifest was given as ${params.cnv_pon_manifest} but that file " +
                  "does not exist. Omit it to look for manifest.json next to --cnv_pon, or fix the path."
        }

        common_bed_ch = Channel.value(
            params.cnv_common_bed ? file(params.cnv_common_bed, checkIfExists: true) : noFile)

        CNV_MANIFEST_CHECK(
            pon_ch,
            Channel.value(file(params.cnv_targets_bed, checkIfExists: true)),
            common_bed_ch,
            Channel.value(manifestFile.exists() ? manifestFile : noFile),
            Channel.value(file("${params.scriptdir}/cnv_manifest_check.py", checkIfExists: true))
        )

        CNV_FILTER(
            CNVKIT_POSTPROCESS.out.calls.map   { _meta, cns -> cns }.collect(),
            CNVKIT_POSTPROCESS.out.metrics.map { _meta, m   -> m   }.collect(),
            CNVKIT_POSTPROCESS.out.cnr.map     { _meta, cnr -> cnr }.collect(),
            sex_tsv_ch,
            filter_ch,
            pon_ch,
            panel_bed_ch,
            Channel.value(file(params.cnv_gnomad_sv, checkIfExists: true)),
            Channel.value(file(params.cnv_clingen,   checkIfExists: true)),
            Channel.value(file(params.cnv_segdup,    checkIfExists: true)),
            Channel.value(params.cnv_gnomad_exome_cnv ? file(params.cnv_gnomad_exome_cnv, checkIfExists: true) : noFile),
            common_bed_ch,
            CNV_MANIFEST_CHECK.out.provenance
        )

        CNV_QC_SUMMARY(CNVKIT_POSTPROCESS.out.qc.map { _meta, qc -> qc }.collect())

    emit:
        calls    = CNVKIT_POSTPROCESS.out.calls
        filtered = CNV_FILTER.out.all_calls
        xlsx     = CNV_FILTER.out.xlsx
        qc       = CNV_QC_SUMMARY.out.summary
        provenance = CNV_MANIFEST_CHECK.out.provenance
}
