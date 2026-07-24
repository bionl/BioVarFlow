#!/usr/bin/env python3
"""
Post-hoc rescue of contamination-flagged somatic variants from FilterMutectCalls.

A variant flagged as `contamination` is rescued (emitted as PASS) when ALL of:
  TLOD  >= 10    strong somatic log-odds
  POPAF >= 5.0   gnomAD AF < 1e-5
  NLOD  >= 3.0   normal confidently ref
  GERMQ >= 20    not germline (Phred)
  normal AD[1] == 0  zero alt reads in matched normal

FORMAT columns are assumed ordered [tumor, normal].
Multi-allelic TLOD/POPAF: only the first value is evaluated.
"""

import sys
import gzip

TLOD_MIN  = 10.0
POPAF_MIN = 5.0
NLOD_MIN  = 3.0
GERMQ_MIN = 20.0


def _first_float(value: str) -> float:
    return float(value.split(",")[0])


def _info_dict(info_field: str) -> dict:
    d = {}
    for token in info_field.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        else:
            d[token] = True
    return d


def _rescue(filters: str, info: dict, fmt_keys: list, normal_fmt: str) -> bool:
    """Return True if this contamination-flagged variant should be rescued."""
    active = set(filters.split(";")) - {"PASS", "."}
    if "contamination" not in active:
        return False

    try:
        if _first_float(info.get("TLOD",  "0")) < TLOD_MIN:  return False
        if _first_float(info.get("POPAF", "0")) < POPAF_MIN: return False
        if _first_float(info.get("NLOD",  "0")) < NLOD_MIN:  return False
        if _first_float(info.get("GERMQ", "0")) < GERMQ_MIN: return False

        # Normal sample AD[1] must be 0
        fmt_vals = normal_fmt.split(":")
        ad_idx   = fmt_keys.index("AD")
        ad_parts = fmt_vals[ad_idx].split(",")
        if int(ad_parts[1]) != 0:
            return False

    except (ValueError, IndexError, KeyError):
        return False

    return True


def process(in_path: str, out_path: str) -> None:
    opener = gzip.open if in_path.endswith(".gz") else open

    rescued = 0
    total   = 0

    with opener(in_path, "rt") as fin, open(out_path, "w") as fout:
        fmt_keys = []
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT tumor normal
            if len(cols) < 11:
                fout.write(line)
                continue

            total += 1
            filter_col = cols[6]
            info_col   = cols[7]
            fmt_col    = cols[8]
            normal_col = cols[10]   # normal is the second sample

            fmt_keys = fmt_col.split(":")
            info     = _info_dict(info_col)

            if _rescue(filter_col, info, fmt_keys, normal_col):
                cols[6] = "PASS"
                rescued += 1

            fout.write("\t".join(cols) + "\n")

    print(f"[mutect2_rescue] {rescued}/{total} contamination-flagged variants rescued",
          file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.vcf[.gz] output.vcf")
    process(sys.argv[1], sys.argv[2])
