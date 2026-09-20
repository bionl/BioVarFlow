#!/usr/bin/env python3
"""
Merge per-sample variant sheets into one cross-case table.

Source sheet per case
---------------------
Exomiser ranks a case's variants against its HPO terms, so where phenotype data
exists the 'Prioritised' sheet is the right spine: ~12 rows instead of ~29, and
it includes findings outside the 285-gene panel (NUP214 on IQMM, for instance).

Exomiser will not run at all without at least one HPO term -- verified against
15.1.0: a bare --vcf is refused with "No sample specified!", and a phenopacket
with no phenotypicFeatures fails to parse its own subject id. A generic root
term (HP:0000001) does run, but every phenotype score collapses to 0.000 and the
top combined score drops from 0.70 to 0.005, so nothing clears the threshold and
the prioritised list comes back empty. There is no useful middle ground.

So cases without HPO terms fall back to the panel Reportable sheet, and the
Ranked_By column records which happened:

  exomiser        Prioritised sheet, phenotype-ranked
  exomiser_none   Exomiser ran, nothing met the criteria (case listed, 0 rows)
  panel_only      no HPO terms; rows come from the Reportable sheet

The column matters because the two sheets do not mean the same thing. Without
it, a reader comparing 12 rows for one case against 29 for another would read
the second as having more findings, when it actually had less filtering.

Every case appears in the Cases sheet whatever happened, so a case can never
drop out of the deliverable unremarked.
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

PRIORITISED = "Prioritised"
REPORTABLE = "HemOnc (Reportable)"
ACMG = "ACMG SF (Reportable)"
# ENST -> RefSeq, derived from the NCBI MANE summary (see --mane).
DEFAULT_MANE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "data", "mane_enst_to_refseq.tsv")

# Emitted for every row regardless of which sheet it came from.
# The first six reproduce the agreed deliverable format exactly; the rest are
# carried because they answer the questions a reviewer asks next -- what kind of
# change, how common, what ClinVar says, whether QC flagged it, and where
# Exomiser ranked it.
OUT_COLS = [
    "Case #", "Case ID", "Gene", "Variant", "HGVSc", "MANE_Select",
    "Panel",
    "HGVSp", "Consequence", "Zygosity", "ClinVar", "ClinVar_Stars", "gnomAD_AF",
    "StrandBias", "Exo_Rank", "Exo_Score", "Exo_Disease",
    "Ranked_By", "Scope",
]

# Kept by --filter coding. VEP emits several per record, '&'-joined, so these
# are matched as substrings against the whole field.
CODING_CONSEQUENCES = (
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "frameshift_truncation", "frameshift_elongation",
    "inframe_insertion", "inframe_deletion",
    "splice_acceptor_variant", "splice_donor_variant",
    "protein_altering_variant", "transcript_ablation",
    "transcript_amplification",
)


def case_id_for(path, xl):
    """Sample_ID from the Sample Summary sheet, else the filename stem."""
    ss = xl.get("Sample Summary")
    if ss is not None and len(ss) and "Sample_ID" in ss.columns:
        val = ss["Sample_ID"].iloc[0]
        if pd.notna(val) and str(val).strip():
            return str(val).strip()
    return re.sub(r"_variants$", "", os.path.splitext(os.path.basename(path))[0])


def split_hgvs(value):
    """Exomiser's HGVS field -> (transcript, c., p.).

    Formatted GENE:TRANSCRIPT:c.NNN:p.(Xxx) -- e.g.
    MSH2:ENST00000233146.7:c.942+2T>G:p.?
    Note the transcript is an Ensembl ENST, where the Reportable sheet's
    MANE_ID is a RefSeq NM_. Same transcript, different accession namespace;
    they are deliberately not reconciled here.
    """
    parts = str(value or "").split(":")
    tx = c = p = None
    for part in parts:
        if part.startswith("ENST") or part.startswith("NM_"):
            tx = part
        elif part.startswith("c."):
            c = part
        elif part.startswith("p."):
            p = part
    return tx, c, p


def panel_lookup(xl):
    """{Variant: row} from 'PASS variants' -- the per-sample annotation source.

    Exomiser reports transcripts as Ensembl ENST; the deliverable uses RefSeq
    MANE_Select (NM_). Rather than translate between accession namespaces, take
    HGVSc/HGVSp/MANE_ID straight from the workbook's own VEP annotation, which
    is where the Reportable sheet gets them too. Off-panel Exomiser hits are
    absent here by construction and keep Exomiser's own values.
    """
    pv = xl.get("PASS variants")
    if pv is None or not len(pv) or "Variant" not in pv.columns:
        return {}
    return {str(r["Variant"]): r for _, r in pv.iterrows()}


def load_mane(path):
    """{ENST without version: (RefSeq_nuc, symbol)} from the MANE summary.

    Off-panel Exomiser hits never reach VEP -- BedFilterVCF removes them before
    annotation -- so the only transcript available for them is Exomiser's own
    Ensembl accession. Clinical interpretation is done against RefSeq, so map
    them here. Keyed without the version suffix because Exomiser's ENST version
    and MANE's need not agree; the accession itself is stable.
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0].startswith("ENST"):
                out[f[0].split(".")[0]] = (f[1], f[2])
    return out


def to_refseq(tx, gene, mane):
    """ENST -> MANE Select RefSeq, when we can do it safely.

    The gene symbol must agree, otherwise the mapping is left alone: a silent
    swap to the wrong transcript is far worse than an ENST a reviewer can look
    up. On the current data all 124 off-panel rows mapped with zero symbol
    mismatches.
    """
    if not tx or not str(tx).startswith("ENST"):
        return tx
    hit = mane.get(str(tx).split(".")[0])
    if hit and (not gene or str(gene) == hit[1]):
        return hit[0]
    return tx


def exomiser_lookup(xl):
    """{Variant: row} from the Exomiser tab, for attaching rank/score."""
    ex = xl.get("Exomiser")
    if ex is None or not len(ex) or "Variant" not in ex.columns:
        return {}
    return {str(r["Variant"]): r for _, r in ex.iterrows()}


def rows_superset(xl, mane, reportable_sheet):
    """Every Reportable variant, plus Exomiser's off-panel hits.

    Prioritised alone is NOT safe as the spine. It is built from Exomiser's
    output, and Exomiser filters before it ranks -- ~570 of ~240,000 variants
    survive for a typical case, with intronic, synonymous and common ones
    discarded. On IQMM that left 27 of 29 HemOnc and 17 of 22 ACMG SF
    Reportable variants absent from Prioritised, including a Pathogenic GATA2
    in another case. Panel membership already justifies a row; Exomiser is
    supporting evidence, not a gate.

    So: take the Reportable sheets whole, attach Exomiser rank/score where it
    happened to rank them, and add the off-panel hits Exomiser found on top.
    """
    exo = exomiser_lookup(xl)
    out, seen = [], {}

    for sheet, tag in ((ACMG, "ACMG SF"), (reportable_sheet, "HemOnc")):
        df = xl.get(sheet)
        if df is None or not len(df):
            continue
        for _, r in df.iterrows():
            key = str(r.get("Variant"))
            if key in seen:            # the 8 genes listed on both panels
                if tag not in seen[key]["Panel"]:
                    seen[key]["Panel"] += "+" + tag
                continue
            e = exo.get(key)
            row = {
                "Scope": "PANEL_REPORTABLE",
                "Panel": tag,
                "Gene": r.get("Gene"),
                "Variant": r.get("Variant"),
                "HGVSc": r.get("HGVSc"),
                "HGVSp": r.get("HGVSp"),
                "MANE_Select": r.get("MANE_ID"),
                "Consequence": r.get("Consequence"),
                "Zygosity": r.get("Zygosity"),
                "ClinVar": r.get("ClinVar"),
                "ClinVar_Stars": r.get("ClinVar_Stars"),
                "gnomAD_AF": r.get("gnomAD_AF"),
                "StrandBias": r.get("StrandBias"),
                # Blank where Exomiser never ranked it -- a fact about
                # Exomiser's filters, not a reason to drop a panel finding.
                "Exo_Rank": e.get("Rank") if e is not None else None,
                "Exo_Score": e.get("Combined_Score") if e is not None else None,
                "Exo_Disease": e.get("Exomiser_Disease") if e is not None else None,
            }
            seen[key] = row
            out.append(row)

    pri = xl.get(PRIORITISED)
    if pri is not None and len(pri) and "Source" in pri.columns:
        for _, r in pri[pri["Source"] == "EXOMISER_ONLY"].iterrows():
            key = str(r.get("Variant"))
            if key in seen:
                continue
            tx, c, pp = split_hgvs(r.get("HGVS"))
            row = {
                "Scope": "OFF_PANEL",
                "Panel": "",
                "Gene": r.get("Gene"),
                "Variant": r.get("Variant"),
                "HGVSc": c,
                "HGVSp": pp,
                "MANE_Select": to_refseq(tx, r.get("Gene"), mane),
                "Consequence": r.get("Consequence"),
                "Zygosity": r.get("Genotype"),
                "ClinVar": r.get("ClinVar"),
                "ClinVar_Stars": r.get("ClinVar_Stars"),
                "gnomAD_AF": None,
                "StrandBias": r.get("StrandBias"),
                "Exo_Rank": r.get("Rank"),
                "Exo_Score": r.get("Combined_Score"),
                "Exo_Disease": r.get("Exomiser_Disease"),
            }
            seen[key] = row
            out.append(row)
    return out


def rows_from_prioritised(df, panel, mane):
    out = []
    for _, r in df.iterrows():
        tx, c, p = split_hgvs(r.get("HGVS"))
        key = str(r.get("Variant"))
        ann = panel.get(key)
        out.append({
            "Scope": r.get("Source"),
            "Gene": r.get("Gene"),
            "Variant": r.get("Variant"),
            # Prefer the workbook's RefSeq annotation; fall back to Exomiser's.
            "HGVSc": (ann.get("HGVSc") if ann is not None else None) or c,
            "HGVSp": (ann.get("HGVSp") if ann is not None else None) or p,
            "MANE_Select": ((ann.get("MANE_ID") if ann is not None else None)
                            or to_refseq(tx, r.get("Gene"), mane)),
            "Consequence": (ann.get("Consequence") if ann is not None else None)
                           or r.get("Consequence"),
            "Zygosity": (ann.get("Zygosity") if ann is not None else None)
                        or r.get("Genotype"),
            "ClinVar": (ann.get("ClinVar") if ann is not None else None) or r.get("ClinVar"),
            "ClinVar_Stars": ann.get("ClinVar_Stars") if ann is not None else r.get("ClinVar_Stars"),
            "gnomAD_AF": ann.get("gnomAD_AF") if ann is not None else None,
            "StrandBias": ann.get("StrandBias") if ann is not None else None,
            "Exo_Rank": r.get("Rank"),
            "Exo_Score": r.get("Combined_Score"),
            "Exo_Disease": r.get("Exomiser_Disease"),
        })
    return out


def rows_from_reportable(df):
    out = []
    for _, r in df.iterrows():
        out.append({
            # Everything on this sheet is on the panel by construction.
            "Scope": "PANEL",
            "Gene": r.get("Gene"),
            "Variant": r.get("Variant"),
            "HGVSc": r.get("HGVSc"),
            "HGVSp": r.get("HGVSp"),
            "MANE_Select": r.get("MANE_ID"),
            "Consequence": r.get("Consequence"),
            "Zygosity": r.get("Zygosity"),
            "ClinVar": r.get("ClinVar"),
            "ClinVar_Stars": r.get("ClinVar_Stars"),
            "gnomAD_AF": r.get("gnomAD_AF"),
            "StrandBias": r.get("StrandBias"),
            "Exo_Rank": None,
            "Exo_Score": None,
            "Exo_Disease": None,
        })
    return out


def keep_coding(rows, spliceai_min, df=None):
    """Protein-altering or splice-affecting rows only."""
    kept = []
    for i, row in enumerate(rows):
        csq = str(row.get("Consequence") or "").lower()
        if any(c in csq for c in CODING_CONSEQUENCES):
            kept.append(row)
            continue
        # Deep-intronic splice effects have a benign-looking consequence, so a
        # high SpliceAI score rescues them where the column is available.
        if df is not None and "SpliceAI_DS_max" in df.columns:
            v = pd.to_numeric(df["SpliceAI_DS_max"].iloc[i], errors="coerce")
            if pd.notna(v) and v >= spliceai_min:
                kept.append(row)
    return kept


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="Directories to scan for *_variants.xlsx, and/or explicit .xlsx paths.")
    ap.add_argument("-o", "--out", required=True, help="Output .xlsx path")
    ap.add_argument("--filter", choices=["all", "coding"], default="all",
                    help="Rows to keep per case (default: all)")
    ap.add_argument("--spliceai-min", type=float, default=0.2,
                    help="With --filter coding, keep anything scoring at least this on "
                         "SpliceAI_DS_max regardless of consequence (default: 0.2)")
    ap.add_argument("--source", choices=["auto", "prioritised", "reportable"], default="auto",
                    help="auto (default) emits every Reportable variant plus Exomiser's "
                         "off-panel hits -- a superset, so no panel finding can be lost. "
                         "'prioritised' uses only the Prioritised sheet (Exomiser-filtered, "
                         "drops most panel variants); 'reportable' uses only the panel sheet.")
    ap.add_argument("--mane", default=DEFAULT_MANE,
                    help="TSV mapping Ensembl transcripts to RefSeq: "
                         "ENST<tab>NM_<tab>symbol, from the NCBI MANE summary. Used to convert "
                         "off-panel Exomiser transcripts, which are Ensembl, to the RefSeq "
                         "accessions used for interpretation. Pass '' to disable.")
    ap.add_argument("--reportable-sheet", default=REPORTABLE,
                    help=f"Fallback sheet name (default: '{REPORTABLE}')")
    args = ap.parse_args()

    paths = []
    for item in args.inputs:
        if os.path.isdir(item):
            # *variants*.xlsx, not *_variants.xlsx: files that have been through
            # a browser download arrive as "AHQX_variants (1).xlsx", and silently
            # skipping those loses whole cases from the deliverable.
            paths.extend(p for p in glob.glob(os.path.join(item, "*variants*.xlsx"))
                         if not os.path.basename(p).startswith("~$"))
        else:
            paths.append(item)
    paths = sorted(set(paths))
    if not paths:
        sys.exit("❌ No *_variants.xlsx found in the given inputs.")

    mane = load_mane(args.mane)
    if args.mane and not mane:
        print(f"⚠️  no MANE mapping loaded from {args.mane} — off-panel rows keep Ensembl IDs")

    variant_rows, case_rows, skipped = [], [], []

    for path in paths:
        try:
            xl = pd.read_excel(path, sheet_name=None)
        except Exception as e:
            skipped.append((path, f"unreadable: {e}"))
            continue

        case = case_id_for(path, xl)
        pri = xl.get(PRIORITISED)
        rep = xl.get(args.reportable_sheet)

        if args.source == "auto":
            # Superset: every Reportable variant plus Exomiser's off-panel hits.
            rows = rows_superset(xl, mane, args.reportable_sheet)
            ranked_by = ("exomiser" if (pri is not None and len(pri)) else "panel_only")
            src_df = rep
            if not rows and rep is None:
                skipped.append((path, f"no '{args.reportable_sheet}' or '{ACMG}' sheet"))
                continue
        elif args.source == "prioritised" and pri is not None and len(pri):
            ranked_by, src_df, rows = "exomiser", pri, rows_from_prioritised(pri, panel_lookup(xl), mane)
        elif args.source == "prioritised" and pri is not None:
            # Exomiser ran and produced nothing -- a real result, distinct from
            # having no phenotype data at all.
            ranked_by, src_df, rows = "exomiser_none", pri, []
        elif rep is not None:
            ranked_by, src_df, rows = "panel_only", rep, rows_from_reportable(rep)
        else:
            skipped.append((path, f"neither '{PRIORITISED}' nor '{args.reportable_sheet}'"))
            continue

        if rows and args.filter == "coding":
            rows = keep_coding(rows, args.spliceai_min, src_df)

        for r in rows:
            r["Case ID"] = case
            r["Ranked_By"] = ranked_by
            variant_rows.append(r)

        case_rows.append({"Case ID": case, "Ranked_By": ranked_by,
                          "Variants": len(rows), "Workbook": os.path.basename(path)})
        print(f"  {case:12s} {ranked_by:14s} {len(rows):4d} rows   ({os.path.basename(path)})")

    if not case_rows:
        sys.exit("❌ No cases collected — check the input paths.")

    cases = pd.DataFrame(case_rows).sort_values("Case ID", kind="stable").reset_index(drop=True)
    # Case # is a 1-based index over cases in sorted Case ID order, so it is
    # stable across reruns for an unchanged input set.
    order = {c: i for i, c in enumerate(cases["Case ID"], start=1)}
    cases.insert(0, "Case #", cases["Case ID"].map(order))

    if variant_rows:
        out = pd.DataFrame(variant_rows)
        out["Case #"] = out["Case ID"].map(order)
        out = (out.reindex(columns=OUT_COLS)
                  .sort_values(["Case #", "Exo_Rank", "Gene", "Variant"], kind="stable")
                  .reset_index(drop=True))
    else:
        out = pd.DataFrame(columns=OUT_COLS)

    with pd.ExcelWriter(args.out) as xw:
        out.to_excel(xw, index=False, sheet_name="Genetic Variants")
        cases.to_excel(xw, index=False, sheet_name="Cases")

    n_pan = int((cases["Ranked_By"] == "panel_only").sum())
    print(f"\n✓ {len(out)} rows from {len(cases)} cases (filter={args.filter}) → {args.out}")
    print(f"  phenotype-ranked: {int((cases['Ranked_By'] == 'exomiser').sum())}   "
          f"no Exomiser result: {int((cases['Ranked_By'] == 'exomiser_none').sum())}   "
          f"panel-only fallback: {n_pan}")
    if n_pan:
        print("  ⚠️  panel_only cases carry more rows per case and are NOT phenotype-ranked; "
              "see the Ranked_By column.")
    if skipped:
        print("\n⚠️  skipped:")
        for path, why in skipped:
            print(f"    {os.path.basename(path)}: {why}")


if __name__ == "__main__":
    main()
