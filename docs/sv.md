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
structural variants. Two params do that, and `nextflow.config` sets both from
`cnv_targets_bed`:

| param | effect |
|---|---|
| `wes = true` | sarek turns this into Manta's `--exome` |
| `intervals = cnv_targets_bed` | sarek bgzips/tabixes it into Manta's `--callRegions` |

These used to be set only inside the CNVkit branch, and only when a panel of
normals was supplied. They are now set whenever `cnvkit` **or** `manta` is in
`params.tools`, so Manta stays exome-correct on a run that has no PoN and
therefore skips CNV calling.

Consequence worth knowing: on a no-PoN run, SNV/indel calling is now also
restricted to `cnv_targets_bed` instead of running across the whole reference.
That is the correct behaviour for a capture experiment — and it is what a
with-PoN run already did — but it is a change for any no-PoN run that relied
on off-target calls.

Manta is skipped in `somatic_mode`, which overrides `params.tools` to
`mutect2` before this wiring runs.

## Not covered

Annotation, filtering and reporting of the SV calls. The VCF is produced and
published; nothing downstream reads it yet. VEP's SV consequence annotation,
a gnomAD-SV frequency filter (the `cnv_gnomad_sv` resource is already in the
config for CNVkit and applies equally here), and a report section are the
next steps.
