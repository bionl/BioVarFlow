# CNV calling

Germline CNV calling for exome and panel projects. Sarek computes CNVkit
coverage; BioVarFlow adds the post-processing that makes the calls reportable.
Validated on the KHCC project — see [cnv_validation.md](cnv_validation.md).

Building panels of normals is **out of scope**. The customer, or a separate PoN
pipeline, delivers a ready CNVkit reference `.cnn`. This pipeline checks that it
fits the data and uses it.

## Turning it on

Nothing to turn on: `cnvkit` is in the default `tools`, so CNV calling is part
of the standard germline run. The **only** input it needs is the project's panel
of normals:

```bash
nextflow run main.nf -profile bionl \
  --input samplesheet.csv --outdir results \
  --cnv_pon refs/pon_final.cnn
```

Everything else -- the capture kit BED, the gene panel, the annotation
resources and every threshold -- is fixed in `nextflow.config`, not asked of the
user. `nextflow_schema.json` therefore exposes `cnv_pon` alone.

**Without `--cnv_pon`, CNV calling is skipped** and the rest of the run is
unaffected: `cnvkit` is dropped from the tool list and a warning naming the
missing param is printed. This pipeline does not build panels of normals, so a
project with no PoN yet still gets its SNV results.

With a PoN supplied, this rewires sarek:

| sarek param | set to | why |
|---|---|---|
| `wes` | `true` | without it CNVkit bins exome data as WGS at 20 kb: bivar ~ 1.0 and false whole-chromosome calls |
| `intervals` | `cnv_targets_bed` | CNVkit hybrid mode needs the capture kit's target BED |
| `cnvkit_reference` | `cnv_pon` | stops sarek building its own flat reference from the fasta |

Somatic mode forces `tools = mutect2`, so CNV calling does not run there.

## Inputs

### From the user

| Param | Notes |
|---|---|
| `cnv_pon` | CNVkit reference `.cnn` built from this project's own normals, by the customer or a separate PoN pipeline. Must come from the same BED as `cnv_targets_bed`, which the pipeline checks before using it. Omit to skip CNV calling. |

### Fixed in nextflow.config

Per project -- change these when the project changes capture kit or panel:

| Param | Default |
|---|---|
| `cnv_targets_bed` | `data/xgen_v2_targets.hg38.cnvkit.bed` (IDT xGen Exome v2, hg38). 4 columns, gene names without exon suffixes. **The PoN must be built from this exact BED.** |
| `cnv_panel_bed` | `data/ACMG_SF_plus_HemOnc_MANE_exons_50bp.refseq.bed` -- the genes a CNV must hit to reach Tier 1. The gene panel, not the capture kit. |
| `cnv_panel_name` | `ACMG SF + HemOnc`, the label on the report |
| `cnv_common_bed` | `null`. Cohort-recurrent CNV regions, if the PoN was delivered with them; left null the recurrent-region filter is not applied. |
| `cnv_pon_haploid_x` | `false`. Set `true` only if the PoN was built with `cnvkit reference -y`. |
| `cnv_pon_size` | `null`. Normals in the PoN, for the low-confidence flag. Left null, no call is flagged on PoN size. |

Shared by every hg38 project, alongside the VEP references:

| Param | Default |
|---|---|
| `cnv_gnomad_sv` | `gs://vep-data-vaic/cnv/gnomad_v4.1_sv_DEL_DUP_CNV.tsv` (prepared once, below) |
| `cnv_gnomad_exome_cnv` | `gs://vep-data-vaic/cnv/gnomad.v4.1.cnv.all.vcf.gz` |
| `cnv_clingen` | `gs://vep-data-vaic/cnv/ClinGen_gene_curation_list_GRCh38.tsv` |
| `cnv_segdup` | `gs://vep-data-vaic/cnv/genomicSuperDups.txt.gz` |

**Sample sex is mandatory.** The samplesheet's `sex` column must be `XX` or
`XY` for every sample; `NA` stops the run. With `--intervals` the antitarget
bins are empty, so CNVkit cannot guess sex reliably, and a wrong guess makes a
normal male chrX read as cn = 0.

All thresholds are params (`cnv_max_bivar`, `cnv_min_exons`, `cnv_del_log2`, ...)
with the validated defaults. Change them in the config per project rather than
editing the filter.

## What runs

### CNVKIT_POSTPROCESS (per sample)

Container pinned to CNVkit 0.9.13 — the version the method was validated on.
Sarek ships 0.9.10.

1. **PoN compatibility check** (`scripts/cnv_pon_check.py`) — compares target-bin
   count and coordinates between the sample and the PoN, and stops with a clear
   message if they differ. `cnvkit.py fix` does not check this itself and
   quietly misaligns bins when the PoN came from another BED.
2. `cnvkit.py sex` on the target coverage — recorded for the QC table, not used
   for calling.
3. `cnvkit.py fix` from sarek's `*.targetcoverage.cnn` / `*.antitargetcoverage.cnn`
   against the PoN. **Sarek's own `.cnr` is not used:** with `--intervals` the
   recalibrated CRAMs keep only on-target reads, so every antitarget bin reads
   zero and sarek's `.cnr` carries ~12k fake "neutral" bins. `fix` against a
   real PoN drops them.
4. `cnvkit.py segment -m cbs -t 0.001 --smooth-cbs --drop-low-coverage`.
   `cnvkit batch` exposes none of these, and its defaults missed a 3-exon
   deletion and merged a 5 Mb mosaic gain into a whole arm.
5. `cnvkit.py call --sample-sex <from samplesheet> --center median`, with `-y`
   only when `cnv_pon_haploid_x` is set.
6. `cnvkit.py metrics` and `scatter`.
7. **chrX sanity check** (`scripts/cnv_chrx_check.py`) — against a
   female-normalised PoN a normal male chrX sits at median log2 ≈ −1 and a
   female ≈ 0; a `-y` PoN shifts both up by 1. A sample more than
   `cnv_chrx_tolerance` from its expectation is marked `FAIL` in the QC table.
   This is what catches a badly built PoN: one prototype PoN had females at
   +0.7 and males at −0.33.

### CNV_FILTER (run level)

`scripts/filter_cnv.py` over every sample's `.call.cns`. Four levels; nothing is
dropped silently — each call keeps the reasons it passed or failed.

| Level | Rule | Default |
|---|---|---|
| 1 Sample QC | bivar | ≤ 0.25 |
| | PoN size | < 10 normals → low-confidence flag |
| 2 Technical | Exons in call | ≥ 3 |
| | Germline signal vs expected level | DEL ≤ −0.4, DUP ≥ +0.3 |
| | Mosaic route | ≥ 50 exons and \|Δlog2\| ≥ 0.2; not on haploid X/Y |
| | Segment depth | ≥ 20× |
| | Median PoN spread over the call | < 0.3 |
| 3 Common | Overlap with `cnv_common_bed` | ≥ 50% → filtered |
| | gnomAD-SV PASS, same type, 50% reciprocal | AF > 0.01 → filtered |
| | Segmental duplication | > 50% → flag only |
| 4 Relevance | Panel gene overlap | required for Tier 1 |
| | ClinGen | DEL hitting HI=3 / DUP hitting TS=3 → annotated |

Tiers: **Tier 1** passes 1–3, hits a panel gene, not in a segdup. **Tier 2**
passes 1–3 but has no panel gene, or sits in a segdup. **Tier 3** filtered, with
the reason. Calls with ≤ 5 exons are flagged "small: confirm (MLPA/ddPCR)".

The mosaic route is kept deliberately: blood-derived "germline" samples can
carry somatic clones, and a KHCC sample had mosaic −7q / −14q / −18q / +17q25 /
+21 at 35–50% of cells.

### CNV_QC_SUMMARY

Collects the per-sample QC rows into `cnv/cnv_qc_summary.tsv` and prints a
warning naming any sample that failed the chrX check or disagrees with CNVkit's
own sex guess.

## Outputs

```
results/cnv/
├── <sample>/
│   ├── <sample>.cnr  .cns  .call.cns     CNVkit, re-fixed against the PoN
│   ├── <sample>.metrics.tsv              bivar etc.
│   ├── <sample>.cnv_qc.tsv               sex agreement + chrX status
│   ├── <sample>.sexcheck.tsv
│   └── <sample>.scatter.png
├── <sample>.cnv_filtered.tsv             one row per non-neutral segment
├── cnv_filtered_all.tsv                  all samples
├── cnv_qc_summary.tsv
└── cnv_report.xlsx                       Summary / Tier 1 / Tier 2 / Tier 3 / Thresholds
```

## One-off resource prep

The slim gnomAD SV table, from `gnomad.v4.1.sv.sites.bed.gz` (1.3 GB, gnomAD
downloads → "Genome SV BED"). Keeps DEL/DUP/CNV records ≥ 200 bp. `freq` is AF
(column 49), or CN_NONREF_FREQ (column 62) for multiallelic CNVs.

```bash
gzcat gnomad.v4.1.sv.sites.bed.gz | awk -F'\t' -v OFS='\t' \
  'NR==1{print "#chrom","start","end","svtype","freq","grpmax_af","filter","name"; next}
   ($5=="DEL"||$5=="DUP"||$5=="CNV") && ($3-$2)>=200 {f=($49!=""&&$49!="NA"&&$49!=".")?$49:$62; print $1,$2,$3,$5,f,$625,$626,$4}' \
  > gnomad_v4.1_sv_DEL_DUP_CNV.tsv
```

Check the column numbers against the file header if gnomAD publishes a new
version. Also download "Exome CNV VCF" from the same page for
`--cnv_gnomad_exome_cnv`.

## Known limits

- **One PoN per run.** Customer batches needing different PoNs would need a
  per-sample PoN column in the samplesheet.
- **Report integration not done.** Tier 1 and Tier 2 CNVs are not yet in the
  sample report; they belong in their own CNV section, not mixed with SNV ACMG
  classes.
- **Manta SV filtering not built.** Manta runs with `--exome` under `--wes`, but
  filtering SVs against the panel needs whole gene spans, not the exon BED.
- **TIDDIT** does not fit exome data; keep it out of `--tools`.
- **Off-target reads:** computing coverage from the markduplicates CRAM
  (pre-BQSR) might bring back off-target reads and give the antitarget bins real
  signal. Untested.
