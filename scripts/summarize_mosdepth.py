#!/usr/bin/env python3
"""
Summarize mosdepth outputs into one QC line for WES/WGS reports.

- Mean depth is taken from total_region (if present) else total in *.mosdepth.summary.txt
- %>=10/20/30/50/100x are computed by aggregating counts in *.thresholds.bed(.gz)

Usage:
  summarize_mosdepth.py --prefix SAMPLE \
                        --summary SAMPLE.mosdepth.summary.txt \
                        --thresholds SAMPLE.thresholds.bed.gz \
                        --out SAMPLE_coverage_summary.overall.txt
"""
import argparse
import gzip
import os
import sys
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--prefix", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--thresholds", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()

def read_mean_depth(summary_path: str) -> float:
    df = pd.read_csv(summary_path, sep="\t")
    if "total_region" in df["chrom"].values:
        mean_depth = df.loc[df["chrom"]=="total_region", "mean"].iloc[0]
    elif "total" in df["chrom"].values:
        mean_depth = df.loc[df["chrom"]=="total", "mean"].iloc[0]
    else:
        sys.stderr.write("ERROR: Neither 'total_region' nor 'total' found in summary.\n")
        sys.exit(2)
    return float(mean_depth)

def read_thresholds(th_path: str) -> pd.DataFrame:
    opener = gzip.open if th_path.endswith(".gz") else open
    # autodetect header
    with opener(th_path, "rt") as fh:
        first = fh.readline()
    has_header = first.lower().startswith(("chrom", "#chrom"))
    df = pd.read_csv(th_path, sep="\t", header=0 if has_header else None)
    if not has_header:
        df.columns = ["chrom","start","end","region","10X","20X","30X","50X","100X"]
    return df

def main():
    args = parse_args()

    mean_depth = read_mean_depth(args.summary)
    df = read_thresholds(args.thresholds)

    df["length"] = df["end"] - df["start"]
    total_len = df["length"].sum()

    percentages = {}
    for t in ["10X","20X","30X","50X","100X"]:
        percentages[t] = df[t].sum() / total_len * 100.0

    out_df = pd.DataFrame([{
        "SAMPLE": args.prefix,
        "MeanDepth": round(mean_depth, 2),
        "Pct>=10x": round(percentages["10X"], 2),
        "Pct>=20x": round(percentages["20X"], 2),
        "Pct>=30x": round(percentages["30X"], 2),
        "Pct>=50x": round(percentages["50X"], 2),
        "Pct>=100x": round(percentages["100X"], 2),
    }])

    out_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
