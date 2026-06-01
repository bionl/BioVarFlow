#!/usr/bin/env python3
"""
Build the run-level output manifest TSV consumed by the downstream
DB ingestion pipeline (bionl_lean_db_export).

Inputs
------
--run-id        Run identifier, propagated as the first column of every row.
--sample-info   TSV file (no header) with one row per sample:
                    sample_id<TAB>assay<TAB>final_vcf<TAB>coverage_bam
                The paths must already be the absolute, *published* paths
                under params.outdir (e.g. gs://bucket/results/<run>/...).
--qc-dir        Directory containing one <sample>.qc.json file per sample,
                produced by the QC Gate (modules/db_qc_export.nf).
--out           Output TSV path.

Output
------
A single TSV with the columns:
    run_id, sample_id, assay, final_vcf, coverage_bam, qc_status, qc_recommendation

This script never filters samples — every sample in --sample-info appears
in the manifest. Filtering on qc_recommendation == "READY" is the
responsibility of the downstream DB ingestion pipeline (DB_INGEST subworkflow).
REVIEW-tier samples require manual triage; FAIL-tier samples are not ingested.

QC field resolution (lenient on purpose):
    qc_status         <- json["qc_status"]   or json["status"]
                          (PASS | REVIEW | FAIL)
    qc_recommendation <- json["db_ingestion_recommendation"]
                         or json["recommendation"]
                         (READY | REVIEW | FAIL)

Missing JSONs or missing fields produce empty strings rather than errors,
so the manifest is always emitted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


COLUMNS = [
    "run_id",
    "sample_id",
    "assay",
    "final_vcf",
    "coverage_bam",
    "qc_status",
    "qc_recommendation",
]


def load_sample_info(path: Path) -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                print(
                    f"WARNING: skipping malformed sample-info row: {line!r}",
                    file=sys.stderr,
                )
                continue
            sample, assay, vcf, bam = parts[0], parts[1], parts[2], parts[3]
            rows.append((sample, assay, vcf, bam))
    return rows


def load_qc_for_sample(qc_dir: Path, sample: str) -> tuple[str, str]:
    """Return (qc_status, qc_recommendation) for sample; empty strings if missing."""
    candidate = qc_dir / f"{sample}.qc.json"
    if not candidate.exists():
        print(f"WARNING: no QC JSON for sample {sample!r} at {candidate}", file=sys.stderr)
        return "", ""
    try:
        with candidate.open() as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"WARNING: could not read {candidate}: {exc}", file=sys.stderr)
        return "", ""

    qc_status = data.get("qc_status") or data.get("status") or ""
    qc_recommendation = (
        data.get("db_ingestion_recommendation")
        or data.get("recommendation")
        or ""
    )
    return str(qc_status), str(qc_recommendation)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--sample-info", required=True, type=Path)
    ap.add_argument("--qc-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    sample_rows = load_sample_info(args.sample_info)
    if not sample_rows:
        print("ERROR: no samples found in sample-info file", file=sys.stderr)
        return 2

    sample_rows.sort(key=lambda r: r[0])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for sample, assay, vcf_path, bam_path in sample_rows:
            qc_status, qc_recommendation = load_qc_for_sample(args.qc_dir, sample)
            row = [
                args.run_id,
                sample,
                assay,
                vcf_path,
                bam_path,
                qc_status,
                qc_recommendation,
            ]
            out.write("\t".join(row) + "\n")

    print(
        f"Wrote manifest with {len(sample_rows)} sample(s) to {args.out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
