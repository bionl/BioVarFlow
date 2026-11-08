# CNV post-processing — validation on the KHCC test data

Validates the CNVkit post-processing method on the KHCC test data. The pipeline
wiring that followed is described in [cnv.md](cnv.md); this file is the evidence
behind its defaults.

**Read the caveat in "What this does and does not validate" before quoting any
number here.**

## What was run

CNVkit 0.9.13 (`etal/cnvkit`; the binary is at `/opt/conda/bin/cnvkit.py`, **not** on the
login PATH — a `bash -lc` wrapper in a Nextflow process will not find it).

Per sample, against `new_run/pon_final.A.cnn` (27 normals):

```
cnvkit.py sex     S.targetcoverage.cnn
cnvkit.py fix     S.targetcoverage.cnn S.antitargetcoverage.cnn pon_final.A.cnn -o S.cnr
cnvkit.py segment S.cnr -m cbs -t 0.001 --smooth-cbs --drop-low-coverage -o S.cns
cnvkit.py call    S.cns --sample-sex <sex> --center median -o S.call.cns
cnvkit.py metrics S.cnr -s S.cns -o S.metrics.tsv
```

Then `filter_cnv.py` over all three `.call.cns` with the ACMG SF + HemOnc panel BED,
gnomAD-SV v4.1, gnomAD exome CNV, ClinGen, UCSC segdups and `new_run/common_cnv.bed`.

`fix` keeps 189,674 of 197,769 target bins and **0 of 14,421 antitarget bins** for all three
samples. `cnv_normal`'s `qc/coverage_check.tsv` puts this on the whole cohort:
`antitarget_frac_with_reads` is 0.000 for all 32 samples, so the antitargets of
this capture carry no signal at all and `fix` correctly drops them. The "most
bins have no or very low coverage" warning comes from that and is expected.

This is a property of the capture and of how the CRAMs were made, **not** of
`--intervals`: the only sarek step that subsets reads by intervals is ApplyBQSR,
which never runs in BioVarFlow (`skip_tools = "baserecalibrator"`). Earlier
versions of these docs blamed `--intervals`; that was wrong. The reason to
re-run `fix` and `segment` instead of using sarek's `.cnr`/`.cns` is the
segmentation settings, not the antitargets.

## Level 1 — sample QC

| Sample | bivar | ≤ 0.25 |
|---|---|---|
| IQMM | 0.148 | pass |
| IRIC | 0.244 | pass (mosaic sample, close to the threshold) |
| EYTA | 0.123 | pass |

## Rule 6 — chrX sanity check (median target log2 on the `.cnr`)

| Sample | Sex (samplesheet) | `cnvkit.py sex` | chrX | chrY | autosomes |
|---|---|---|---|---|---|
| IQMM | male | — (chrX confirms male) | −0.995 | +0.171 | −0.001 |
| IRIC | female | Female (+0.262) | −0.009 | −25.592 | −0.001 |
| EYTA | male | Male (−0.752) | −0.979 | +0.178 | −0.000 |

Males ≈ −1, female ≈ 0, autosomes centred: the PoN is correctly female-normalised, so
`-y` must stay off (`cnv_pon_haploid_x = false` for this PoN).

## Expected calls

3 samples, 161 calls: **Tier 1 6, Tier 2 24, Tier 3 131**.

| Sample | Expected | Observed | Verdict |
|---|---|---|---|
| IQMM | chr15:89,210,284–89,211,901 DEL ≈ −0.75, 3 exons (RLBP1), Tier 2 + "small: confirm" | identical coordinates, log2 −0.76, 3 exons, `all_genes=RLBP1`, Tier 2, flag "small: confirm (MLPA/ddPCR)" | match |
| IQMM | chr11:544,719–561,892 DUP ≈ +0.58 (LRRC56/LMNTD2/RASSF7), Tier 2 | identical coordinates, log2 +0.58, 19 exons, `all_genes=LMNTD2,LRRC56,RASSF7`, Tier 2 | match |
| IQMM | chr1:196,825,476–196,918,406 DEL (CFHR1–CFHR4), Tier 3 common in gnomAD-SV | identical coordinates, log2 −1.04, `all_genes=CFHR1,CFHR4`, Tier 3, reason "common in gnomAD-SV (AF 0.015)" | match |
| EYTA | chr2:85,539,287–85,867,617 DUP ≈ +0.58, 83 exons, Tier 2 | chr2:85,541,082–85,867,617, log2 +0.58, 82 exons, Tier 2 | match (start 1.8 kb in, 82 vs 83 exons) |
| IRIC | mosaic 7q DEL 95.5–159.1 Mb | four contiguous segments: 95.5–102.5, 102.7–122.7, 122.7–152.0, 152.0–159.1 Mb, log2 −0.35 to −0.44 | match (same territory, more breakpoints) |
| IRIC | 14q DEL 39.0–45.2 Mb | chr14:39,045,162–45,247,286 DEL, −0.36, 159 exons | match |
| IRIC | 17q25 DUP 78.0–83.1 Mb | chr17:78,049,062–83,094,575 DUP, +0.23, 1051 exons | match |
| IRIC | 18q DEL 51.1–80.2 Mb | chr18:51,078,255–80,247,348 DEL, −0.33, 822 exons | match |
| IRIC | +21 | chr21 DUP in three segments, 13.6–21.5, 25.6–41.3, 41.3–46.7 Mb, +0.28 to +0.35 | match |
| IRIC | Tier 1 where POT1, KCNH2, FLNC, PRKAG2, XRCC2, FANCM, FASN, GAA, SMAD4, RUNX1 are hit | all ten appear in Tier 1 `panel_genes` (FAAP100 and IL10RB additionally) | match |

The mosaic route fired as intended on IRIC: estimated cell fractions 31–46%, in the
35–50% range the prototype reported. The 7q and 21 arms are split across more segments than the
leave-one-out prototype, which is expected: **IQMM and EYTA are inside**
`pon_final.A.cnn`, so their own signal partly normalises itself away and shifts
breakpoints. IRIC is **not** in the PoN — `manifest.json` records
`"excluded": ["IRIC"]` — so its breakpoint spread is segmentation, not
self-normalisation. Every panel gene and every tier assignment still lands where
the prototype put it.

## Change to the ported script

`scripts/filter_cnv.py` is the prototype with two BioVarFlow-local changes, both
tagged `# BioVarFlow-local` so a re-port stays mechanical:

1. `--haploid-x-reference`: with a PoN built with `-y`, `expected()` puts a
   normal male chrX at log2 0 and a female at +1, instead of −1 and 0.
2. `--provenance`: appends `CNV_MANIFEST_CHECK`'s `cnv_pon_provenance.tsv` to the
   workbook's Thresholds sheet, so a delivered report records which normals
   produced it.

Everything else is the prototype verbatim. Re-running the prototype's own
command through the ported script reproduces its numbers exactly: 32 samples,
1517 calls, Tier 1 6, Tier 2 113, Tier 3 1398.

## Comparison with 3billion's reports, all samples

Our CNVkit calls were compared with 3billion's CNV reports for 31 KHCC samples
(`compare_3billion.py` in the prototype; reports in `3billion/cnv_reports/`).

| 3billion depth calls (their caller, 3bCNV) | n | found by us, any tier |
|---|---|---|
| IRIC (mosaic sample), all calls | 52 | 50 |
| Other 29: clear signal in our data (≥ 3 exons, \|log2\| ≥ 0.3 DUP / 0.4 DEL) | 67 | 42 |
| Other 29: homozygous CFHR3/CFHR1 deletion | 22 | **0 before the hom-del scan, 22 after** |
| Other 29: weak or no signal in our data | 40 | 2 |
| Other 29: < 3 exons, or not in the capture | 14 | 1 |
| 3billion **Manta** calls (sub-exon SVs, mostly intronic) | 42 | 1 — out of scope for CNVkit |

Per reported event: IQMM 2 of 4 (`RLBP1` and `LRRC56` found; `PRSS2` and
`CNNM2`/`NT5C2` have weak signal in our data), IRIC 5 of 5 mosaic events, and
EYTA's chr2:85.54–85.87 Mb DUP is in their report too.

Going the other way, 42 of our 118 Tier 1/Tier 2 calls in samples that have a
3billion report also appear in it. Many of the rest sit in segmental
duplications or are small events already flagged "confirm". One recurring set of
DUPs (chr2:95.5 Mb, chr3:130.6 Mb, chr4:69.9–70.6 Mb) appears only in XZJR,
NIOZ, YIOO and IRIC — the samples sequenced between the two technical batches —
so it is probably technical.

## Homozygous deletions

`segment --drop-low-coverage` removes 0-read bins before segmentation, so a
2-copy loss never reaches a segment. CFHR3/CFHR1 sit at log2 ≈ −27 and were
missed in all 7 samples carrying them: 0 of 3billion's 22 calls. The `.cnr` scan
in `CNV_FILTER` recovers them.

- 90 homozygous deletions across the 32 samples, 86 in Tier 3 and 4 in Tier 2.
- CFHR3/CFHR1 found in all 7 samples 3billion reported it in (DVFZ, EXGT, KAXU,
  NYOM, WZIU, XZJR, YIOO) — 22 of 22 of their calls.
- Almost all 90 are known common polymorphisms: `LCE3B`/`LCE3C`, `UGT2B17`,
  `ANKRD36`/`ANKRD36B`, `CFHR3`/`CFHR1`, `GBP3`, `USP17L`, `BTNL3`/`BTNL8`.
- Tier counts over 32 samples went from T1 6 / T2 109 / T3 1312 to
  **T1 6 / T2 113 / T3 1398**. All non-hom-del calls are unchanged.
- Single-sample mode: XZJR alone still reports the `CFHR3`, `CFHR1`,
  `LCE3B`/`LCE3C` and `GBP3` hom-dels in Tier 3, so the fallback to the PoN
  works when the in-run cohort check cannot.

## Segmentation methods tested

CBS at `-t 0.001 --smooth-cbs` is the choice. It still misses about 25 short
events that have a clear signal, mostly 3–6 exon duplications — XDMT's `PADI3`
DUP, 10 exons at +0.57, was left inside one 121 Mb neutral segment. Lowering
`-t` to 0.01 or 0.0001, and dropping `--smooth-cbs`, do not help.

`hmm-germline` was tested on all 32 samples and **rejected**:

| | CBS | hmm-germline |
|---|---|---|
| of the 67 clear-signal 3billion events | 42 | 35 |
| of IRIC's 52 events | 50 | 44 |
| IQMM `RLBP1` deletion (validated) | found | **lost** |
| Tier 2 calls | 117 | 224 |

Running both methods and merging adds 3 events. Not worth doubling Tier 2.

## What this does and does not validate

This is a **method check, not a clinical validation.**

- **27 of the 32 samples are inside `pon.A.cnn`.** A sample in its own PoN
  partly normalises its own CNVs away, so its numbers here are optimistic about
  breakpoints and pessimistic about amplitude. The prototype's leave-one-out
  numbers are the clean ones.
- The samples that are genuinely independent of PoN A are **IRIC** (excluded
  from it) and the four batch-B samples (FBWA, SZRM, TIAM, XDMT) — whose own
  `pon.B.cnn` has 4 normals, below `cnv_min_pon_size`, so their calls carry the
  low-confidence flag.
- So there is no clean end-to-end run on this cohort at all. A real validation
  needs a cohort that is not also the PoN.
- The three-sample run at the top of this file was done before the Nextflow
  wiring and before the hom-del scan; its tier counts (161 calls, T1 6 / T2 24 /
  T3 131) describe that run, not the current pipeline.
