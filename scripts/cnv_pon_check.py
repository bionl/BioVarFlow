#!/usr/bin/env python3
"""Fail fast when the CNVkit PoN was not built from the run's target BED.

`cnvkit.py fix` does not check that its reference describes the same bins as the
sample. Given a PoN built from a different capture BED it either errors in a way
that points nowhere useful, or produces a .cnr whose bins are quietly misaligned.
Comparing target-bin count and coordinates up front turns that into one clear
message naming the two parameters involved.

Usage: cnv_pon_check.py <sample.targetcoverage.cnn> <reference.cnn>
"""
import csv
import sys


def target_bins(path):
    """Target-bin coordinates of a .cnn, in file order, antitargets excluded."""
    with open(path) as fh:
        rows = csv.reader(fh, delimiter="\t")
        head = next(rows)
        ig = head.index("gene") if "gene" in head else 3
        return [(r[0], r[1], r[2]) for r in rows if r[ig] != "Antitarget"]


def main(sample_cnn, pon_cnn):
    s, p = target_bins(sample_cnn), target_bins(pon_cnn)
    hint = ("The PoN must be built from the same BED passed as --cnv_targets_bed.")

    if len(s) != len(p):
        sys.exit(f"ERROR: PoN/target-BED mismatch: {pon_cnn} has {len(p)} target bins, "
                 f"{sample_cnn} has {len(s)}. {hint}")

    bad = [(a, b) for a, b in zip(s, p) if a != b]
    if bad:
        a, b = bad[0]
        sys.exit(f"ERROR: PoN/target-BED mismatch: {len(bad)} of {len(s)} bins differ in "
                 f"coordinates (first: sample {a[0]}:{a[1]}-{a[2]} vs PoN {b[0]}:{b[1]}-{b[2]}). "
                 f"{hint}")

    print(f"PoN check: {len(s)} target bins match")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
