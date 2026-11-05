# Structural variants (Manta)

Manta is on by default: `manta` is in `params.tools`. It needs no user input.

It was added after a method comparison against a 3billion 3B-EXOME report
(`EPF26-XDMT`), which calls SVs with Manta v1.6.0 alongside its depth-based
CNV caller. Depth-based CNV calling and paired-end/split-read SV calling are
orthogonal: CNVkit finds copy-number change but not its breakpoints, Manta
finds breakpoints and also the balanced events (inversions, translocations)
that leave copy number unchanged and are therefore invisible to CNVkit.
A Manta call with an exact breakpoint is the one that can be confirmed by
Sanger.

## What runs

Sarek's `BAM_VARIANT_CALLING_GERMLINE_MANTA` subworkflow, which is
`MANTA_GERMLINE` on the recalibrated/mapped CRAM. Only the diploid SV VCF is
kept and annotated; the candidate small-indel and candidate-SV VCFs are
Manta's own intermediates.

## Output

```
${outdir}/variant_calling/manta/<sample>/<sample>.manta.diploid_sv.vcf.gz
${outdir}/variant_calling/manta/<sample>/<sample>.manta.diploid_sv.vcf.gz.tbi
```

These are **not** in `manifest/run_output_manifest.tsv`, whose columns are
fixed at sample / assay / consensus VCF / BAM. Adding an SV column means
changing `scripts/build_run_output_manifest.py` too.

## Exome wiring

Manta must know it is looking at exome data, or it reports off-target noise as
structural variants. `nextflow.config` sets one param for that, whenever
`cnvkit` or `manta` is in `params.tools`:

| param | effect |
|---|---|
| `wes = true` | sarek turns this into Manta's `--exome`, and into CNVkit's `--method hybrid` |

`params.intervals` is deliberately left unset, so Manta gets **no**
`--callRegions` and scans genome-wide. That is slower and noisier off-target,
but not wrong, and it is the price of not restricting the small-variant
callers.

### Why intervals is not set

Sarek has no per-tool intervals -- one `PREPARE_INTERVALS` feeds every caller
-- and, confirmed by grep across the whole submodule, no interval padding
anywhere. `params.intervals` is used verbatim. Setting it to the capture BED
therefore restricts DeepVariant and HaplotypeCaller to it too.

Measured, `data/ACMG_SF_plus_HemOnc_MANE_exons_50bp.refseq.bed` against
`data/xgen_v2_targets.hg38.cnvkit.bed`:

- 41% of panel bases fall inside the targets, but almost all of the gap is
  3'UTR: 95 panel regions of >=3 kb hold 506 kb of it, the terminal exons of
  `APC`, `CBL`, `SMAD4`, `BMPR1A`, `TTN`. Capture kits skip long UTRs, so
  there are no reads there either way. Not a real loss.
- Restricted to CDS-sized panel regions (<1 kb): 53% covered by the bare
  targets, 90.6% by the targets padded +-50bp. The probes sit flush with the
  coding exon.
- Of the 9,732 exon-boundary +-10bp windows, 46 (0.47%) are fully inside the
  targets, 8,710 (89.50%) only partially, and 976 (10.03%) not at all. The
  probes stop at the exon edge, so the *intronic* half of nearly every splice
  window is outside the BED -- the canonical `+1/+2/-1/-2` positions.

So a bare capture BED as `intervals` means DeepVariant and HaplotypeCaller
never call canonical splice variants across ~90% of the panel's exons.
`BRCA2` has 7 such windows, `MUTYH` 7, `SCN5A` 9, `MYH11` 13. Panel filtering
and Exomiser select *from* the callset, so neither can recover a variant that
was never called.

CNVkit does not need intervals at all: with `--cnvkit_reference` set,
`reference_exists` is true and
`modules/nf-core/cnvkit/batch/main.nf` drops `--targets`, taking bins from the
`.cnn`. The CRAM itself is never interval-restricted either --
`intervals_for_preprocessing` goes only to `CRAM_SAMPLEQC`, and BQSR is
skipped via `skip_tools = baserecalibrator` -- so CNVkit's coverage is
unaffected. (Traced through the channel wiring, not confirmed in a live run.)

### What leaving it unset costs

`PREPARE_INTERVALS` falls back to `BUILD_INTERVALS` on the `.fai`, giving
whole-chromosome intervals. Consequences:

| | effect |
|---|---|
| CNVkit | none -- `--targets` suppressed anyway, hybrid mode from `wes` |
| DeepVariant / HaplotypeCaller | unrestricted, splice sites intact; scatter over the whole genome, so slower |
| Manta | no `--callRegions`, genome-wide scan, slower and noisier off-target |
| sarek on-target QC | **lost** -- mosdepth reports genome-wide coverage, so there is no "% of targets at >=20x" |

The panel-level QC in this pipeline's own QC gate uses `params.bed` and is
unaffected.

Because sarek disables mosdepth's `--fast-mode` whenever `params.wes` is true,
whole-chromosome intervals would make it compute per-base depth across the
genome and emit a per-base BED of that size. `nextflow.config` therefore pins
`MOSDEPTH` back to `-n --fast-mode --by 500`: the on-target metric is already
gone, so the genome-wide QC runs in genome-wide mode.

The alternative, if the on-target QC matters more than caller runtime, is to
set `intervals` to a merged +-100bp padding of `cnv_targets_bed` while leaving
`cnv_targets_bed` itself bare for the PoN. Padding is safe for CNVkit (it
never sees intervals) and recovers the splice sites; panel coverage plateaus
there, 43% -> 67.6% at +-50bp -> 69.0% at +-100bp.

Manta is skipped in `somatic_mode`, which overrides `params.tools` to
`mutect2` before this wiring runs.

## Not covered

Annotation, filtering and reporting of the SV calls. The VCF is produced and
published; nothing downstream reads it yet. VEP's SV consequence annotation,
a gnomAD-SV frequency filter (the `cnv_gnomad_sv` resource is already in the
config for CNVkit and applies equally here), and a report section are the
next steps.
