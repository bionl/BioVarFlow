# CNV calling

Germline CNV calling for exome and panel projects. Sarek computes CNVkit
coverage; BioVarFlow adds the post-processing that makes the calls reportable.
Validated on the KHCC project — see [cnv_validation.md](cnv_validation.md).

BioVarFlow does not build panels of normals. A separate pipeline,
[`cnv_normal`](#building-a-pon), does, and delivers `pon.<batch>.cnn` +
`common_cnv.bed` + `manifest.json`. BioVarFlow checks that delivery against the
data and uses it; see [Building a PoN](#building-a-pon) for the three client
situations and which one applies.

## Turning it on

Nothing to turn on: `cnvkit` is in the default `tools`, so CNV calling is part
of the standard germline run. The **only** input it needs is the project's panel
of normals:

```bash
nextflow run main.nf -profile bionl \
  --input samplesheet.csv --outdir results \
  --cnv_pon refs/pon.A.cnn \
  --cnv_common_bed refs/common_cnv.bed
```

Everything else -- the capture kit BED, the gene panel, the annotation
resources and every threshold -- is fixed in `nextflow.config`, not asked of the
user. `nextflow_schema.json` therefore exposes `cnv_pon` and `cnv_common_bed`
and nothing else.

Those two are **one delivery, not two inputs.** The `common_cnv.bed` is built
from the same normals as the `.cnn` by the `cnv_normal` pipeline, and a BED from
a different PoN is worse than no BED at all. Supply the pair or neither.

**Without `--cnv_pon`, CNV calling is skipped** and the rest of the run is
unaffected: `cnvkit` is dropped from the tool list and a warning naming the
missing param is printed. This pipeline does not build panels of normals, so a
project with no PoN yet still gets its SNV results.

With a PoN supplied, this rewires sarek:

| sarek param | set to | why |
|---|---|---|
| `wes` | `true` | this is what selects CNVkit's `--method hybrid`. Without it CNVkit bins exome data as WGS at 20 kb: bivar ~ 1.0 and false whole-chromosome calls |
| `cnvkit_reference` | `cnv_pon` | stops sarek building its own flat reference from the fasta |

`intervals` is deliberately **not** set. An earlier version of this doc said
CNVkit hybrid mode needs the capture BED as `--intervals`; that is wrong.
`--method hybrid` comes from `wes`, and when `--cnvkit_reference` is given
sarek's cnvkit module sets `reference_exists` and drops `--targets`
altogether, taking its bins from the `.cnn`. Sarek has no per-tool intervals
and no interval padding, so setting `intervals` would restrict every caller
to the BED verbatim and silently drop canonical splice positions from the
SNV/indel callset. See `docs/sv.md` for the measurements.

Somatic mode forces `tools = mutect2`, so CNV calling does not run there.

## Inputs

### From the user

| Param | Notes |
|---|---|
| `cnv_pon` | CNVkit reference `.cnn` built from this project's own normals, by the customer or the `cnv_normal` pipeline. Must come from the same BED as `cnv_targets_bed`. Omit to skip CNV calling. |
| `cnv_common_bed` | The `common_cnv.bed` delivered with that `.cnn`: cohort-recurrent CNV regions, columns chrom, start, end, type (`DEL`, `DUP`, `DEL/DUP`, `HOMDEL`), `n_samples=N`. Supplying `cnv_pon` without it prints a warning, because the recurrent artifacts it catches then reach Tier 1 and Tier 2. |
| `cnv_pon_manifest` | The `manifest.json` from the same build. Leave empty and it is looked for next to `cnv_pon`, which is where a delivery puts it, so it normally needs no flag. Setting it to a path that does not exist is an error, rather than silently falling back. |

### Fixed in nextflow.config

Per project -- change these when the project changes capture kit or panel:

| Param | Default |
|---|---|
| `cnv_targets_bed` | `data/xgen_v2_targets.hg38.cnvkit.bed` (IDT xGen Exome v2, hg38). 4 columns, gene names without exon suffixes. **The PoN must be built from this exact BED.** Not passed to sarek -- it records the PoN's provenance and is the BED to build a new PoN against. |
| `cnv_panel_bed` | `data/ACMG_SF_plus_HemOnc_MANE_exons_50bp.refseq.bed` -- the genes a CNV must hit to reach Tier 1. The gene panel, not the capture kit. |
| `cnv_panel_name` | `ACMG SF + HemOnc`, the label on the report |
| `cnv_pon_haploid_x` | `false`. Set `true` only if the PoN was built with `cnvkit reference -y`. |
| `cnv_pon_size` | `null`. Normals in the PoN, for the low-confidence flag. Left null, no call is flagged on PoN size. |

`cnv_common_bed` moved out of this table: it is a per-PoN input, not a per-project
constant. Measured on KHCC, filtering XZJR against `pon.A.cnn` without it moved
Tier 1 from 0 to 2 (a chr2 `TTN` DUP at +0.36 and a chr2 `COL3A1` DUP at +0.53,
both cohort-recurrent artifacts) and Tier 2 from 7 to 24. EYTA went 0 to 1
Tier 1 and 11 to 14 Tier 2. Its `HOMDEL` rows matter most in single-sample runs,
where the in-run recurrence test (>= 3 samples homozygous) cannot fire.

Shared by every hg38 project, alongside the VEP references:

| Param | Default |
|---|---|
| `cnv_gnomad_sv` | `gs://vep-data-vaic/cnv/gnomad_v4.1_sv_DEL_DUP_CNV.tsv` (prepared once, below) |
| `cnv_gnomad_exome_cnv` | `gs://vep-data-vaic/cnv/gnomad.v4.1.cnv.all.vcf.gz` |
| `cnv_clingen` | `gs://vep-data-vaic/cnv/ClinGen_gene_curation_list_GRCh38.tsv` |
| `cnv_segdup` | `gs://vep-data-vaic/cnv/genomicSuperDups.txt.gz` |

**Sample sex is mandatory.** The samplesheet's `sex` column must be `XX` or
`XY` for every sample; `NA` stops the run. CNVkit's own guess leans on the
antitarget bins, and whether those carry any signal depends on the capture and
on how the CRAM was made -- it is not guaranteed. On KHCC's coverage it is
exactly zero (`antitarget_frac_with_reads` is 0.000 for all 32 samples in
`cnv_normal`'s `qc/coverage_check.tsv`), so there is nothing to guess from. A
wrong sex makes a normal male chrX read as cn = 0, so it is taken from the
samplesheet rather than inferred. `CNVKIT_POSTPROCESS` still records CNVkit's
guess in the QC row, and `CNV_QC_SUMMARY` warns when the two disagree.

All thresholds are params (`cnv_max_bivar`, `cnv_min_exons`, `cnv_del_log2`, ...)
with the validated defaults. Change them in the config per project rather than
editing the filter.

## What runs

### CNV_MANIFEST_CHECK (once per run)

Checks the PoN against its `manifest.json` before anything is called, and writes
`cnv_pon_provenance.tsv`, which `CNV_FILTER` appends to the report's Thresholds
sheet so a delivered workbook records which normals produced it.

Two findings **stop the run**, both because they corrupt results without looking
like errors:

| Check | Why it stops |
|---|---|
| `targets_bed.md5` differs from `cnv_targets_bed` | different bins, so every log2 ratio is wrong. Stronger than `cnv_pon_check.py`, which compares bin counts and so only catches a size difference |
| `reference_sex_model` disagrees with `cnv_pon_haploid_x` | `-y` against a female-normalised PoN makes a normal male chrX read as `cn = 0`, which looks like a reportable homozygous deletion |

These warn and continue, because a usable third-party or older PoN may
legitimately not satisfy them:

- no `manifest.json` at all — the Thresholds sheet then says the provenance is
  unverified, rather than leaving the row out
- `cnvkit_version` differs from `cnv_cnvkit_version`
- no `pons` entry matches the `.cnn` given, so `n_normals` cannot be checked
- the `.cnn`'s own md5 differs from the manifest's, i.e. the file changed since
- `n_normals` below `cnv_min_pon_size` — the low-confidence flag already marks
  the calls themselves
- `cnv_common_bed`'s md5 differs from the manifest's: almost certainly a BED
  from a different PoN, which is worse than supplying none, since its regions
  describe another cohort's artifacts

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
   against the PoN. **Sarek's own `.cnr` is not used,** because its `.cns` and
   calls come from `cnvkit batch`'s defaults rather than the segmentation in
   step 4, and re-running `fix` against the project PoN is what makes the
   `.cnr` comparable to that PoN. A side effect on PoNs whose antitargets are
   empty (KHCC's): `fix` drops every antitarget bin, rather than carrying ~14k
   zero-coverage bins through as "neutral". See Known limits.
4. `cnvkit.py segment -m cbs -t 0.001 --smooth-cbs --drop-low-coverage`.
   `cnvkit batch` exposes none of these, and its defaults missed a 3-exon
   deletion and merged a 5 Mb mosaic gain into a whole arm.
   Note what `--drop-low-coverage` costs: it removes 0-read bins *before*
   segmentation, so a homozygous deletion never reaches a segment at all. That
   is why `CNV_FILTER` scans the `.cnr` separately for them (see its rules
   table).
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
| | *Homozygous deletions are exempt from level 2* | 0 reads is its own evidence |
| | Germline signal vs expected level | DEL ≤ −0.4, DUP ≥ +0.3 |
| | Mosaic route | ≥ 50 exons and \|Δlog2\| ≥ 0.2; not on haploid X/Y |
| | Segment depth | ≥ 20× |
| | Median PoN spread over the call | < 0.3 |
| 3 Common | Overlap with `cnv_common_bed` | ≥ 50% → filtered |
| | gnomAD-SV PASS, same type, 50% reciprocal | AF > 0.01 → filtered |
| | Segmental duplication | > 50% → flag only |
| 4 Relevance | Panel gene overlap | required for Tier 1 |
| | ClinGen | DEL hitting HI=3 / DUP hitting TS=3 → annotated |

Plus a separate scan for events level 2 can never see:

| Homozygous deletion | Rule | Default |
|---|---|---|
| Detection | consecutive `.cnr` target bins at log2 below the expected level | ≥ 2 bins, log2 < −3 |
| Confirmation | other samples with normal coverage (log2 > −1) at those bins | ≥ 50% of them |
| Small runs | below this many other samples, the cohort check is skipped and the PoN alone is trusted | < 5 |
| Frequency | gnomAD-SV **AF²** (the homozygous genotype frequency), matched by containment rather than reciprocal overlap | AF² > 0.01 → common |
| Recurrence | samples in the run homozygous at the same place | ≥ 3 → recurrent |

Containment, not reciprocal overlap, because an exon-bounded call is shorter
than the real event. The `HOMDEL` rows in `cnv_common_bed` matter most in
single-sample runs, where the in-run recurrence test cannot fire: measured on
XZJR alone, dropping them let a `GPAT2` hom-del escape to Tier 2, because its
gnomAD-SV AF of 0.0052 is under the `cnv_max_af` ceiling. `CFHR3`/`CFHR1` stay
Tier 3 either way at AF 0.2293.

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

## Coverage must match the PoN's

A PoN is only valid for coverage computed the same way it was. BioVarFlow takes
its coverage from sarek's `cnvkit batch` on the preprocessed CRAM; `cnv_normal`
computes it with `cnvkit.py coverage` directly. They agree because both leave
these at CNVkit's defaults, and `cnv_normal` pins them explicitly so a future
CNVkit version cannot silently move one:

| | value | where |
|---|---|---|
| antitarget average size | 150000 | `cnv_normal`: `antitarget_avg_size`; sarek: CNVkit default |
| minimum MAPQ | 0 | `cnv_normal`: `min_mapq = 0`, so no `--min-mapq` is passed |
| metric | read depth, not read counts | `cnv_normal`: `count_reads = false`, so no `--count` |
| target BED | `cnv_targets_bed` | md5-checked against the manifest |

Two of these only matter at PoN build time. Once a PoN exists, sarek's
`cnvkit batch` runs with `--reference`, which suppresses `--targets` and takes
every bin from the `.cnn` — so the patient run inherits the PoN's binning
whatever `antitarget_avg_size` would have been. MAPQ and the depth-vs-counts
choice do apply to the patient's own coverage, and a mismatch there shifts every
log2 ratio.

Any PoN from outside `cnv_normal` has to match this table. `CNV_MANIFEST_CHECK`
can only verify what the manifest records, so for a third-party PoN the coverage
settings are an assumption, not a check.

## Building a PoN

Use the `cnv_normal` pipeline, branch `pon-build-loo` (a separate repository,
not a submodule of this one). Do not build a reference by hand: the `.cnn` alone is not a delivery, and `common_cnv.bed` and
`manifest.json` are what make it checkable and the tiers trustworthy.

Which flow applies depends on what the client brings:

| Client brings | Mode | Flow |
|---|---|---|
| A ready PoN (`.cnn` + `common_cnv.bed` + `manifest.json`) | A | BioVarFlow with `--cnv_pon --cnv_common_bed` |
| Normals (FASTQ) + patients | B | normals through sarek preprocessing, then `cnv_normal`, then BioVarFlow on the patients with that PoN |
| Only patients (>= 20, same kit and lab) | C | patients through BioVarFlow **without** a PoN (CNV skipped), then `cnv_normal` on their CRAMs; the leave-one-out calls are their CNV results |

### Mode C: the cohort is its own PoN

Mode C needs no change in BioVarFlow — with `cnv_pon` unset, `cnvkit` is dropped
from the tool list, CNV calling is skipped and the SNV side is unaffected.

**Do not then re-run those same samples through BioVarFlow against the PoN built
from them.** A sample inside its own PoN partly normalises its own CNVs away, so
real events come back at reduced amplitude or not at all. Their results are
`cnv_normal`'s leave-one-out calls, in `loo/<sample>/<sample>.call.cns`, where
each sample was held out of the reference used to call it.

To tier those calls, run the filter directly on them:

```bash
python3 scripts/filter_cnv.py \
  --calls 'loo/*/*.call.cns' --cnr 'loo/*/*.cnr' \
  --sex sex.tsv --batch batch.tsv \
  --pon A=pon.A.cnn --pon-size A=27 \
  --panel-bed data/ACMG_SF_plus_HemOnc_MANE_exons_50bp.refseq.bed \
  --panel-name 'ACMG SF + HemOnc' \
  --common-bed common_cnv.bed \
  --gnomad-sv ... --clingen ... --segdup ... \
  --outdir out --out-xlsx out/cnv_report.xlsx
```

A `CNV_FILTER`-only entry point (`--cnv_loo_dir`) would make this a pipeline step
rather than a manual command. Not built.

This caveat is not only about mode C. It is why the KHCC numbers in
[cnv_validation.md](cnv_validation.md) need reading carefully: 27 of the 32
samples are inside `pon.A.cnn`, so a pipeline run on them is partly
self-normalised. IRIC is the exception — `manifest.json` records
`"excluded": ["IRIC"]`, so its result is clean.

## Known limits

- **One PoN per run,** and KHCC already needs two. Batch B (FBWA, SZRM, TIAM,
  XDMT) is technically distinct: FBWA's bivar is 0.28 against the full cohort
  and 0.17 against its own batch. `cnv_normal` therefore builds `pon.A.cnn` (27
  normals) and `pon.B.cnn` (4). Until there is a per-sample PoN column in the
  samplesheet, run batch B separately with `--cnv_pon pon.B.cnn`; its 4 normals
  are below `cnv_min_pon_size`, so its calls carry the low-confidence flag.
- **No clean end-to-end validation on KHCC.** 27 of the 32 samples are inside
  `pon.A.cnn`, so a pipeline run on them is partly self-normalised. What is left
  is IRIC (excluded from the PoN, so clean) and the four batch-B samples (whose
  PoN is below the size threshold). The tier counts in
  [cnv_validation.md](cnv_validation.md) are a method check, not a clinical
  validation; that needs a cohort that is not also the PoN.
- **CBS misses some short events that have a clear signal.** About 25 against
  3billion on KHCC, mostly 3-6 exon duplications; XDMT's `PADI3` DUP, 10 exons
  at +0.57, was left inside one 121 Mb neutral segment. Lowering `-t` to 0.01 or
  0.0001, and dropping `--smooth-cbs`, do not help. `hmm-germline` was tested on
  all 32 samples and **rejected**: it found 35 of the 67 clear-signal events
  against CBS's 42, lost the validated IQMM `RLBP1` deletion, found 44 of IRIC's
  52 against CBS's 50, and doubled Tier 2 (224 calls against 117). Running both
  methods adds 3 events. This is a known limit, not a TODO.
- **The homozygous-deletion scan depends on the PoN.** A gene with no coverage
  in the PoN cannot be called, because `cnvkit fix` drops bins the PoN does not
  cover before the scan ever sees them.
- **Report integration not done.** Tier 1 and Tier 2 CNVs are not yet in the
  sample report; they belong in their own CNV section, not mixed with SNV ACMG
  classes.
- **Manta SV filtering not built.** Manta runs with `--exome` under `--wes`, but
  filtering SVs against the panel needs whole gene spans, not the exon BED.
- **TIDDIT** does not fit exome data; keep it out of `--tools`.
- **Off-target reads are not used, and that is not caused by `--intervals`.**
  Whether the antitarget bins carry signal depends on the capture and on how the
  CRAM was made. `params.intervals` never restricted this pipeline's CRAMs: the
  only sarek step that subsets reads by intervals is ApplyBQSR
  (`external/sarek/modules/nf-core/gatk4/applybqsr/main.nf`), and
  `skip_tools = "baserecalibrator"` means it never runs. So removing
  `params.intervals` changed nothing here.

  On KHCC's coverage the antitargets are simply empty --
  `antitarget_frac_with_reads` is 0.000 for all 32 samples in `cnv_normal`'s
  `qc/coverage_check.tsv` -- so `fix` drops all ~14k antitarget bins and calling
  is target-only. That is harmless, not a defect. The off-target signal is only
  usable when both the PoN's normals and the patient have it; check
  `coverage_check.tsv` for a new capture before assuming either way.
