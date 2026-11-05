# CNV post-processing — validation on the KHCC test data

Validates steps 2–3 of `CNV_postprocessing_agent_prompt.md` after porting
`filter_cnv.py` into `scripts/`. No Nextflow wiring yet.

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
samples — rule 4's prediction, and the reason sarek's own `.cnr` must not be used. The
"most bins have no or very low coverage" warning comes from the empty antitargets and is expected.

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
leave-one-out prototype, which is expected — IQMM, IRIC and EYTA are all *inside* `pon_final.A.cnn`,
so their own signal partly normalises itself away and shifts breakpoints. Every panel gene and
every tier assignment still lands where the prototype put it.

## Change to the ported script

`scripts/filter_cnv.py` is the prototype unchanged except for `--haploid-x-reference`
(rule 5 / step 3 of the prompt): with a PoN built with `-y`, `expected()` now puts a normal
male chrX at log2 0 and a female at +1, instead of −1 and 0.

## Not yet done

Nextflow wiring: params, sarek `--wes --intervals --cnvkit_reference`, the
`CNVKIT_POSTPROCESS` and `CNV_FILTER` processes, the PoN-compatibility fail-fast check,
and report integration (step 4, "ask first").
