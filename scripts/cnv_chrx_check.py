#!/usr/bin/env python3
"""chrX sanity check on a CNVkit .cnr, plus a per-sample CNV QC row.

Against a female-normalised PoN a normal male chrX sits at median target log2
~ -1 and a normal female at ~ 0; a PoN built with -y is already haploid on X, so
both shift up by 1. A sample far from its expectation means the samplesheet sex
is wrong, or the PoN's own X normalisation is not what --cnv_pon_haploid_x says
it is -- the condition that made a whole prototype PoN unusable. Either way the
sample's calls should not be read before someone looks, so this writes a status
column rather than failing the run.

CNVkit's own `cnvkit.py sex` guess is recorded for comparison but not trusted:
with --intervals the antitarget bins are empty and that guess gets weak.

Usage:
  cnv_chrx_check.py --cnr S.cnr --sample S --sex male|female
                    [--sexcheck S.sexcheck.tsv] [--haploid-x-reference]
                    [--tolerance 0.35] --out S.cnv_qc.tsv
"""
import argparse
import csv
import statistics as st

ap = argparse.ArgumentParser()
ap.add_argument("--cnr", required=True)
ap.add_argument("--sample", required=True)
ap.add_argument("--sex", required=True, choices=["male", "female"])
ap.add_argument("--sexcheck", help="output of `cnvkit.py sex`, for the agreement column")
ap.add_argument("--haploid-x-reference", action="store_true",
                help="the PoN was built with -y, so expected chrX log2 shifts up by 1")
ap.add_argument("--tolerance", type=float, default=0.35)
ap.add_argument("--out", required=True)
a = ap.parse_args()

COLS = ["sample", "sex_samplesheet", "sex_cnvkit_guess", "sex_agree", "chrx_median_log2",
        "chrx_expected_log2", "chrx_delta", "chrx_status", "chry_median_log2",
        "autosome_median_log2", "chrx_bins"]

vals = {"chrX": [], "chrY": [], "auto": []}
with open(a.cnr) as fh:
    rows = csv.reader(fh, delimiter="\t")
    head = next(rows)
    ic, ig, il = head.index("chromosome"), head.index("gene"), head.index("log2")
    for r in rows:
        if r[ig] == "Antitarget":      # always empty under --intervals; never informative
            continue
        vals[r[ic] if r[ic] in ("chrX", "chrY") else "auto"].append(float(r[il]))

med = {k: (st.median(v) if v else float("nan")) for k, v in vals.items()}

if a.haploid_x_reference:
    expected = 0.0 if a.sex == "male" else 1.0
else:
    expected = -1.0 if a.sex == "male" else 0.0
delta = med["chrX"] - expected
status = "PASS" if abs(delta) <= a.tolerance else "FAIL"

guessed = "unknown"
if a.sexcheck:
    try:
        lines = [l for l in open(a.sexcheck).read().splitlines() if l.strip()]
        guessed = lines[-1].split("\t")[1].strip().lower()
    except (IndexError, OSError):
        pass

with open(a.out, "w") as out:
    out.write("\t".join(COLS) + "\n")
    # "unknown", not "no", when cnvkit gave no usable guess -- the QC summary
    # flags disagreement, and a missing guess is not a disagreement.
    agree = "unknown" if guessed == "unknown" else ("yes" if guessed == a.sex else "no")
    out.write("\t".join(str(x) for x in [
        a.sample, a.sex, guessed, agree,
        f"{med['chrX']:.3f}", f"{expected:.3f}", f"{delta:+.3f}", status,
        f"{med['chrY']:.3f}", f"{med['auto']:.3f}", len(vals["chrX"]),
    ]) + "\n")

print(f"chrX check: median {med['chrX']:+.3f}, expected {expected:+.3f} -> {status}")
if status == "FAIL":
    print(f"WARNING: {a.sample} chrX is {delta:+.3f} off its expectation. Either the "
          f"samplesheet sex is wrong, or the PoN's X normalisation does not match "
          f"--cnv_pon_haploid_x={a.haploid_x_reference}.")
if guessed not in (a.sex, "unknown"):
    print(f"WARNING: {a.sample} samplesheet says {a.sex}, cnvkit guessed {guessed}.")
