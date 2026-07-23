#!/usr/bin/env python3
import sys, os, re, argparse
import pandas as pd
from cyvcf2 import VCF
from urllib.parse import quote_plus

# -------------------------
# CLI
# -------------------------
p = argparse.ArgumentParser(description="Build LEAN Excel deliverable")
p.add_argument("vcf", help="VCF (bgzip+tabix) with VEP+ClinVar if available")
p.add_argument("coverage_summary", help="exon_coverage_summary.sorted.txt")
p.add_argument("r1r2_summary", help="r1r2_per_exon.tsv")
p.add_argument("fr_summary", help="frstrand_per_exon.tsv")
p.add_argument("xlsx_out", help="output Excel path, e.g. results/sample.xlsx")

# Optional QC inputs for Sample Summary (all optional)
p.add_argument("--sample-id", default=None, help="Sample ID override (else inferred from VCF header or filename)")
p.add_argument("--assay", default="", help="WES/WGS (free text)")
p.add_argument("--build", default="GRCh38", help="Reference build")
p.add_argument("--flagstat", default=None, help="samtools flagstat output")
p.add_argument("--stats", default=None, help="samtools stats output")
p.add_argument("--picard-align", default=None, help="Picard CollectAlignmentSummaryMetrics.txt")
p.add_argument("--picard-insert", default=None, help="Picard CollectInsertSizeMetrics.txt")
p.add_argument("--mosdepth-summary", default=None, help="mosdepth *.summary.txt (plain, not gz)")
p.add_argument("--verifybamid", default=None, help="verifyBamID2 *.selfSM or .summary output (FREEMIX)")
p.add_argument("--sexcheck", default=None, help="file with inferred sex (one line like: Sex: Male)")
p.add_argument("--sf-genes", default=None, help="ACMG SF gene list (one gene SYMBOL per line)")
p.add_argument("--hemonc-genes", default=None,
               help="HemOnc gene list (one gene SYMBOL per line). "
                    "When provided, three extra sheets are written: 'HemOnc (P-LP)', "
                    "'HemOnc Coverage gaps', 'HemOnc Genes Coverage'.")
p.add_argument("--acmg-thresholds", default=None,
               help="mosdepth thresholds.bed(.gz) run with ACMG BED (cols: chrom start end [region] 20X 30X)")
p.add_argument("--gaps20", default=None, help="annotated gaps <20x BED (chrom start end RegionLabel)")
p.add_argument("--gaps30", default=None, help="annotated gaps <30x BED (chrom start end RegionLabel)")



args = p.parse_args()

# -------------------------
# Helpers: safe parsers
# -------------------------
def safe_float(x):
    try: return float(x)
    except: return None

def sample_from_vcf(vcf_path):
    try:
        v = VCF(vcf_path)
        if v.samples: return v.samples[0]
    except: pass
    # fallback to filename stem
    return os.path.basename(vcf_path).split(".")[0]

def parse_flagstat(path):
    """Return dict with Total_Reads, Mapped_Reads, Mapped_Percent, Duplicate_Reads, Duplicate_Percent"""
    out = {}
    if not path or not os.path.exists(path): return out
    with open(path) as f:
        for line in f:
            s = line.strip()
            # total line is often like: "123456 + 0 in total ..."
            if " in total " in s and "QC-passed" in s:
                n = s.split()[0]
                if n.isdigit(): out["Total_Reads"] = int(n)
            # mapped line: "123456 + 0 mapped (99.97% : N/A)"
            elif re.search(r"\bmapped\s*\(", s):
                n = s.split()[0]
                out["Mapped_Reads"] = int(n) if n.isdigit() else None
                m = re.search(r"\(([\d\.]+)%", s)
                if m: out["Mapped_Percent"] = safe_float(m.group(1))
            # duplicates line (may be absent if not marked)
            elif re.search(r"\bduplicates\b", s) and "(" in s:
                n = s.split()[0]
                out["Duplicate_Reads"] = int(n) if n.isdigit() else None
                m = re.search(r"\(([\d\.]+)%", s)
                if m: out["Duplicate_Percent"] = safe_float(m.group(1))
    return out

def parse_samtools_stats(path):
    """Return Ti/Tv, Het/Hom (approx) from samtools stats if present"""
    out = {}
    if not path or not os.path.exists(path): return out
    ti = tv = het = homalt = mean_insert_size = 0
    with open(path) as f:
        for line in f:
            # SN  number of first fragments:   123
            if line.startswith("SN"):
                if "number of heterozygous sites" in line:
                    het = int(line.strip().split()[-1])
                if "number of homozygous-ALT sites" in line:
                    homalt = int(line.strip().split()[-1])
                if "TSTV" in line:
                    # some builds print: "SN\tts/tv\t2.02"
                    parts = line.strip().split()
                    try: out["TiTv"] = float(parts[-1])
                    except: pass
                if "insert size average" in line:
                    mean_insert_size = float(line.strip().split()[-1])  
            # For bcftools stats you’d parse different keys; kept minimal here.
    if het or homalt:
        try: out["Het_Hom_Ratio"] = round(het / max(homalt,1), 3)
        except: pass
    if mean_insert_size:
        out["Mean_Insert_Size"] = mean_insert_size
    return out

def parse_picard_alignment(path):
    """Picard CollectAlignmentSummaryMetrics: %Q30, %mapped, mean read length, etc."""
    out = {}
    if not path or not os.path.exists(path): return out
    # Picard has a header section; data lines often start with CATEGORY
    # We'll pick the FIRST row for 'PAIR' or 'UNPAIRED' if present
    df = None
    try:
        with open(path) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
        # Find header line (tab-delimited)
        hdr_idx = next(i for i,l in enumerate(lines) if l.startswith("CATEGORY"))
        header = lines[hdr_idx].split("\t")
        row = lines[hdr_idx+1].split("\t")
        d = dict(zip(header,row))
        # Common keys (presence depends on Picard version)
        out["Pct_Q30"] = safe_float(d.get("PF_READS_ALIGNED_Q20_BASES", ""))  # fallback if Q30 missing
        if "PCT_PF_READS_ALIGNED" in d: out["Pct_Mapped"] = safe_float(d["PCT_PF_READS_ALIGNED"])*100 if d["PCT_PF_READS_ALIGNED"] not in (None,"") else None
        if "MEAN_READ_LENGTH" in d: out["Mean_Read_Length"] = safe_float(d["MEAN_READ_LENGTH"])
        if "MEAN_INSERT_SIZE" in d: out["Mean_Insert_Size"] = safe_float(d["MEAN_INSERT_SIZE"])
        if "PCT_ADAPTER" in d: out["Pct_Adapter"] = safe_float(d["PCT_ADAPTER"])*100 if d["PCT_ADAPTER"] not in (None,"") else None
    except Exception:
        pass
    return out

def parse_picard_insert(path):
    """Picard InsertSizeMetrics: mean insert size if not captured above."""
    out = {}
    if not path or not os.path.exists(path): return out
    try:
        with open(path) as f:
            lines = [l.strip() for l in f if l.strip()]
        hdr_idx = next(i for i,l in enumerate(lines) if l.startswith("MEDIAN_INSERT_SIZE"))
        header = lines[hdr_idx].split("\t")
        row = lines[hdr_idx+1].split("\t")
        d = dict(zip(header,row))
        if "MEAN_INSERT_SIZE" in d:
            out["Mean_Insert_Size"] = safe_float(d["MEAN_INSERT_SIZE"])
    except Exception:
        pass
    return out

def parse_mosdepth_summary(path):
    """Return Mean_Coverage and %>=10x/20x/30x if mosdepth summary available."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return out
        row = df.iloc[0]  # assume one sample per file

        out["Mean_Coverage"] = float(row.get("MeanDepth", 0))

        for col in ["Pct>=10x","Pct>=20x","Pct>=30x","Pct>=50x","Pct>=100x"]:
            if col in df.columns:
                try:
                    out[col] = float(row[col])
                except Exception:
                    pass
    except Exception:
        # silently ignore parse errors
        return out
    return out
def parse_mosdepth_summary_regions(path):
    """Parse mosdepth summary.txt to get mean coverage and min/max across regions."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return out

        # Focus on on-target regions (rows ending with _region)
        region_df = df[df["chrom"].astype(str).str.endswith("_region")].copy()
        if region_df.empty:
            region_df = df.copy()

        # Compute weighted mean coverage across regions
        cov_mean = (region_df["mean"] * region_df["length"]).sum() / region_df["length"].sum()
        out["Mean_Coverage"] = round(float(cov_mean), 2)

    except Exception as e:
        print(f"[WARN] Failed to parse mosdepth summary {path}: {e}")
        return out
    return out
def _num(x):
    # handles "78,05" and "100%" etc.
    return pd.to_numeric(str(x).replace("%","").replace(",", "."), errors="coerce")

import pandas as pd

def build_genes_coverage(df_in, sf_genes=None):
    """
    Expect columns EXACTLY: Gene, MANE_ID, ExonLen, Pct>=20x
    Return per-gene table: Gene, MANE_ID, Pct20  (length-weighted mean).
    """
    required = ["Gene","MANE_ID","ExonLen","Pct>=20x"]
    missing = [c for c in required if c not in df_in.columns]
    if missing:
        raise ValueError(f"Genes Coverage: missing columns: {', '.join(missing)}")

    df = df_in[required].copy()

    # Optional: restrict to ACMG set
    if sf_genes:
        df = df[df["Gene"].isin(sf_genes)]

    # Clean & coerce
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df = df[df["Gene"].notna() & (~df["Gene"].isin(["","NA","."]))]

    df["MANE_ID"] = df["MANE_ID"].astype(str).str.strip()

    df["ExonLen"] = pd.to_numeric(df["ExonLen"], errors="coerce")

    # Handle values like "78,05" or "100%"
    p20 = (df["Pct>=20x"].astype(str)
                     .str.replace("%","", regex=False)
                     .str.replace(",",".", regex=False))
    df["Pct>=20x"] = pd.to_numeric(p20, errors="coerce")

    df = df.dropna(subset=["ExonLen","Pct>=20x"])
    df = df[df["ExonLen"] > 0]

    if df.empty:
        return pd.DataFrame(columns=["Gene","MANE_ID","Pct20"])

    # Choose transcript per gene = MANE_ID with largest total exon length (ignore NA/.)
    tx_pref = (df[~df["MANE_ID"].isin(["","NA","."])]
                 .groupby(["Gene","MANE_ID"], as_index=False)["ExonLen"].sum()
                 .sort_values(["Gene","ExonLen"], ascending=[True, False]))
    best_tx = (tx_pref.groupby("Gene", as_index=False)
                      .first()[["Gene","MANE_ID"]])

    # Fallback to "NA" if a gene only had NA transcripts
    all_genes = pd.DataFrame({"Gene": sorted(df["Gene"].unique())})
    best_tx = (all_genes.merge(best_tx, on="Gene", how="left")
                        .fillna({"MANE_ID":"NA"}))

    # Length-weighted mean %≥20x across all regions of each gene
    df["_covered20_bp"] = df["ExonLen"] * (df["Pct>=20x"] / 100.0)
    agg = (df.groupby("Gene", as_index=False)
             .agg(total_len=("ExonLen","sum"),
                  covered20_bp=("_covered20_bp","sum")))
    agg["Pct20"] = (agg["covered20_bp"] / agg["total_len"] * 100.0).round(2)

    out = (agg.merge(best_tx, on="Gene", how="left")
             [["Gene","MANE_ID","Pct20"]]
             .sort_values("Gene")
             .reset_index(drop=True))

    return out


def read_acmg_thresholds(th_path, sf_genes=None):
    """
    Read mosdepth thresholds generated WITH ACMG BED (thresholds 20,30).
    If 'region' column is missing, synthesize RegionLabel as 'chrom:start-end'.
    Returns per-exon rows with RegionLabel, Chrom, ExonStart, ExonEnd, ExonLen, Pct>=20x, Pct>=30x.
    """
    import pandas as pd, os
    if not th_path or not os.path.exists(th_path):
        return pd.DataFrame(columns=["RegionLabel","Chrom","ExonStart","ExonEnd","ExonLen","Pct>=20x","Pct>=30x","Gene","MANE_ID","Exon"])

    df = pd.read_csv(th_path, sep="\t", comment="#", low_memory=False)
    # Normalize column names
    cols = {c.lower(): c for c in df.columns}
    def col(name, default=None): return cols.get(name.lower(), default)

    # If header is missing / different, coerce
    if col("chrom") is None or col("start") is None or col("end") is None:
        # assume fixed layout: chrom start end region 20X 30X ...
        df.columns = ["chrom","start","end","region","20X","30X"] + list(df.columns[6:])
        cols = {c.lower(): c for c in df.columns}

    chrom = col("chrom") or "chrom"
    start = col("start") or "start"
    end   = col("end")   or "end"
    region= col("region")  # may be None if no label in BED

    # Synthesize RegionLabel if missing/blank
    if (region is None) or (region == 'NA') or (region not in df.columns):
        df["RegionLabel"] = df[chrom].astype(str) + ":" + df[start].astype(int).astype(str) + "-" + df[end].astype(int).astype(str)
    else:
        # If region exists but is empty/unknown, still fall back to coords
        df["RegionLabel"] = df[region].astype(str)
        mask_empty = df["RegionLabel"].isin(["", ".", "unknown", "UNKNOWN", "None"])
        df.loc[mask_empty, "RegionLabel"] = df.loc[mask_empty, chrom].astype(str) + ":" + df.loc[mask_empty, start].astype(int).astype(str) + "-" + df.loc[mask_empty, end].astype(int).astype(str)

    # Optional gene parsing (will be '.' if no names)
    parts = df["RegionLabel"].astype(str).str.split("|", n=2, expand=True)
    df["Gene"] = parts[0] if parts.shape[1] > 0 else "."
    df["MANE_ID"] = parts[1] if parts.shape[1] > 1 else "."
    df["Exon"] = parts[2] if parts.shape[1] > 2 else "."

    if sf_genes:
        df = df[df["Gene"].isin(sf_genes)]

    df["ExonLen"] = df[end] - df[start]
    if "20X" not in df.columns or "30X" not in df.columns:
        # if thresholds file has more cols, pick those with digits '20' / '30'
        for t in ("20","30"):
            cand = [c for c in df.columns if "".join(ch for ch in str(c) if ch.isdigit()) == t]
            if cand:
                df[f"{t}X"] = df[cand[0]]

    df["Pct>=20x"] = (df["20X"] / df["ExonLen"] * 100).round(2)
    df["Pct>=30x"] = (df["30X"] / df["ExonLen"] * 100).round(2)

    out = df.rename(columns={chrom:"Chrom", start:"ExonStart", end:"ExonEnd"})[
        ["RegionLabel","Gene","MANE_ID","Exon","Chrom","ExonStart","ExonEnd","ExonLen","Pct>=20x","Pct>=30x"]
    ].drop_duplicates(subset=["RegionLabel","Chrom","ExonStart","ExonEnd"])
    return out


import os
import pandas as pd

def acmg_pct_regions_covered(th_path, min_depth=20, gene_set=None, key_prefix="ACMG"):
    """
    Compute % of regions fully covered at >= min_depth using mosdepth thresholds output.
    Expects columns like: #chrom, start, end, region, 10X, 20X, 30X, ...

    Parameters
    ----------
    gene_set : iterable[str] or None
        If provided, restrict the calculation to rows whose `region` column
        starts with a GENE symbol in this set (assumes `GENE|TRANSCRIPT` labels
        in the source BED). If None, all regions are used (legacy behaviour).
    key_prefix : str
        Prefix for the returned metric keys. Default "ACMG" preserves the
        legacy keys ("ACMG_PctRegionsCovered", ...). Pass e.g. "HemOnc" to
        produce parallel HemOnc-scoped metrics.

    Returns
    -------
    dict with keys: {prefix}_PctRegionsCovered, {prefix}_RegionsCovered, {prefix}_RegionsTotal
    """
    out = {}
    if not th_path or not os.path.exists(th_path):
        return out

    # 1) Read WITH header; infer gzip; do NOT treat '#' as comment so we keep the header row
    df = pd.read_csv(th_path, sep="\t", compression="infer", header=0)

    # 2) Normalize column names (strip leading '#', lower for lookups)
    df.columns = [str(c).lstrip("#") for c in df.columns]
    cols = list(df.columns)

    # 3) Identify coordinate and region columns
    def pick_exact(names, default=None):
        for n in cols:
            if str(n).lower() in names:
                return n
        return default

    chrom = pick_exact({"chrom", "chr"}, cols[0] if cols else None)
    start = pick_exact({"start"}, cols[1] if len(cols) > 1 else None)
    end   = pick_exact({"end"},   cols[2] if len(cols) > 2 else None)
    region_col = pick_exact({"region", "target", "name"}, None)  # mosdepth uses 'region' when --by had a 4th column

    if chrom is None or start is None or end is None:
        return out

    # 3b) Optional gene-set scoping: keep only rows whose region label's gene
    #     portion (before the first '|') is in gene_set.
    if gene_set:
        if region_col and region_col in df.columns:
            gene_series = df[region_col].astype(str).str.split("|", n=1).str[0]
            df = df[gene_series.isin(set(gene_set))].copy()
        # If no region column exists we can't scope by gene → keep behaviour
        # equivalent to "no rows match" so metrics reflect the empty scope.
        else:
            df = df.iloc[0:0].copy()

    # 4) Find the threshold column (prefer exact "20X", else any col whose digits == 20)
    preferred = f"{int(min_depth)}X"
    if preferred in df.columns:
        thr_col = preferred
    else:
        thr_col = None
        for c in cols:
            digits = "".join(ch for ch in str(c) if ch.isdigit())
            if digits and int(digits) == int(min_depth):
                thr_col = c
                break
    if thr_col is None:
        return out  # threshold column not present

    # 5) Coerce numeric and compute interval length
    s = pd.to_numeric(df[start], errors="coerce")
    e = pd.to_numeric(df[end],   errors="coerce")
    thr = pd.to_numeric(df[thr_col], errors="coerce")
    df = df.assign(_len=(e - s), _thr=thr).dropna(subset=["_len"]).copy()

    # 6) Define "region id" to aggregate:
    #    - If 'region' exists, treat each *label* as one region (matches mosdepth behavior when --by had names)
    #    - Otherwise, use exact coordinates as the region id
    if region_col and region_col in df.columns:
        grp = df.groupby(region_col, as_index=False).agg(len_sum=("_len","sum"), thr_sum=("_thr","sum"))
    else:
        grp = (df.assign(_id = df[chrom].astype(str)+":"+
                               s.astype("Int64").astype(str)+"-"+
                               e.astype("Int64").astype(str))
                 .groupby("_id", as_index=False).agg(len_sum=("_len","sum"), thr_sum=("_thr","sum")))

    total_regions = int(len(grp))
    if total_regions == 0:
        return out

    # 7) A region is "fully covered" if all bases in that region have depth >= min_depth:
    covered = int((grp["thr_sum"] >= grp["len_sum"]).sum())
    pct = round(covered / total_regions * 100.0, 2)

    out[f"{key_prefix}_PctRegionsCovered"] = pct
    out[f"{key_prefix}_RegionsCovered"] = covered
    out[f"{key_prefix}_RegionsTotal"] = total_regions
    return out



def load_gap_bed_for_agg(path, threshold_tag, sf_genes=None):
    """
    Read gaps BED (chrom start end [RegionLabel]).
    If RegionLabel is missing, synthesize it from chrom:start-end of the *exon* is not known —
    but since we intersected with the exon BED when annotating, col4 should already be RegionLabel=exon coords.
    Aggregates per RegionLabel: counts, total bp, and intervals list.
    """
    import pandas as pd, os
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["RegionLabel", f"Gaps{threshold_tag}_n", f"Gaps{threshold_tag}_bp", f"Gaps{threshold_tag}_intervals"])

    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0:str})
    # If only 3 cols, synthesize RegionLabel from the *gap* coords (less ideal). We prefer 4th col from annotated step.
    if df.shape[1] >= 4:
        df.columns = ["Chrom","Start","End","RegionLabel"] + [f"extra{i}" for i in range(df.shape[1]-4)]
    else:
        df.columns = ["Chrom","Start","End"]
        df["RegionLabel"] = df["Chrom"].astype(str) + ":" + df["Start"].astype(int).astype(str) + "-" + df["End"].astype(int).astype(str)

    # Optional filter by SF genes if RegionLabel encodes gene (it won't here, so skip)
    df["GapLen"] = pd.to_numeric(df["End"], errors="coerce") - pd.to_numeric(df["Start"], errors="coerce")
    df = df.dropna(subset=["Start","End","GapLen"])
    df = df.drop_duplicates(subset=["Chrom","Start","End","RegionLabel"])
    def _intervals(g):
        return ";".join(f"{r.Chrom}:{int(r.Start)}-{int(r.End)}" for _, r in g.sort_values(["Chrom","Start","End"]).iterrows())

    agg = (df.groupby("RegionLabel", as_index=False)
             .agg(**{
                 f"Gaps{threshold_tag}_n": ("GapLen","size"),
                 f"Gaps{threshold_tag}_bp": ("GapLen","sum"),
             }))
    ints = (df.groupby("RegionLabel").apply(_intervals)
             .reset_index(name=f"Gaps{threshold_tag}_intervals"))
    out = agg.merge(ints, on="RegionLabel", how="left")
    out[f"Gaps{threshold_tag}_n"] = out[f"Gaps{threshold_tag}_n"].astype("Int64")
    out[f"Gaps{threshold_tag}_bp"] = out[f"Gaps{threshold_tag}_bp"].astype("Int64")
    return out


def parse_verifybamid(path):
    """Return FREEMIX contamination if VerifyBamID2 present."""
    out = {}
    if not path or not os.path.exists(path): return out
    try:
        with open(path) as f:
            for l in f:
                if "FREEMIX" in l:
                    m = re.search(r"FREEMIX[:=\s]+([\d\.eE+-]+)", l)
                    if m: out["FREEMIX"] = safe_float(m.group(1))
                # selfSM format: first line often: SAMPLE\tFREEMIX\t... ; handle tabular too
                if "\t" in l and l.upper().startswith("SAMPLE"):
                    # next line contains number
                    row = next(f).strip().split("\t")
                    for i,h in enumerate(l.strip().split("\t")):
                        if h.upper()=="FREEMIX":
                            out["FREEMIX"] = safe_float(row[i])
                            break
                    break
    except Exception:
        pass
    return out

def parse_sexcheck(path):
    if not path or not os.path.exists(path): return {}
    txt = open(path).read().strip()
    m = re.search(r"sex\s*[:=]\s*(\w+)", txt, re.IGNORECASE)
    return {"Sex_Check": m.group(1)} if m else {"Sex_Check": txt}
import re


def _norm_revstat(s: str) -> str:
    """
    Decode VEP CSQ encoding of ClinVar review status text.

    VEP writes CSQ subfields with `_` in place of spaces and `&` in place
    of commas (because `,` is CSQ's transcript separator). So a raw
    ClinVar status like "criteria provided, conflicting classifications"
    arrives here as "criteria_provided&_conflicting_classifications".
    """
    if not s: return ""
    s = s.lower()
    # Reverse VEP encoding: & -> "," and _ -> " "
    s = s.replace("&", ",").replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s

def clinvar_stars_from_revstat(revstat: str) -> int:
    """
    Map ClinVar review status text to 0–4 stars (germline).

    Handles multiple statuses separated by | or ; — returns the maximum.
    Commas are part of ClinVar status text (e.g. "criteria provided,
    conflicting classifications") and MUST NOT be treated as a
    separator; we normalize (decode VEP encoding) first so that any `&`
    from VEP becomes a real `,` inside the status string, then split
    only on genuine multi-status separators.
    """
    if not revstat: return None
    normalized = _norm_revstat(revstat)
    toks = re.split(r"[|;]+", normalized)
    def map_one(t: str) -> int:
        t = t.strip()
        if not t: return -1
        if t.startswith("practice guideline"): return 4
        if t.startswith("reviewed by expert panel"): return 3
        if t.startswith("criteria provided, multiple submitters, no conflicts"): return 2
        if t.startswith("criteria provided, conflicting"): return 1
        if t.startswith("criteria provided, single submitter"): return 1
        if t.startswith("no assertion criteria provided"): return 0
        if t.startswith("no classification provided"): return 0
        if t.startswith("no classification for the individual variant"): return 0
        # default conservative
        return 0
    return max(map_one(t) for t in toks if t is not None)

def clinvar_star_glyph(n):
    return "★"*int(n) if isinstance(n,(int,float)) and n>=0 else ""

def clinvar_url_from_row(r):
    """
    Build a ClinVar URL for a variant row.
    Prefers VCV ID, then ALLELEID; falls back to a search on HGVS/Variant.
    """
    def first_token(x):
        if x is None: return ""
        # take first token if values are delimited by | , ;
        return re.split(r"[|,;]", str(x))[0].strip()

    # Try VCV-style id if present (e.g., VCV000123456 or 123456 -> padded)
    #vcv = first_token(r.get("ClinVar_VCV") or r.get("VCV") or r.get("CLNVID"))
    #print(vcv)
    #if vcv:
    #    vcv = vcv.strip()
    #    if vcv.startswith("VCV"):
    #        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vcv}/"
    #    if vcv.isdigit():
    #        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV{int(vcv):09d}/"

    # Else try ALLELEID
    alleleid = first_token(r.get("ALLELEID"))
    if alleleid:
        # keep only digits
        m = re.search(r"\d+", alleleid)
        if m:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{int(m.group(0))}/"

    # Fallback: search by HGVS (or Variant string)
    #query = (r.get("HGVSc") or r.get("Variant") or "").strip()
    #print(query)
    #if query:
    #    return "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + quote_plus(query.split(":")[1])

    return ""
def select_csq_entry(var, csq_format, prefer_clinvar=True):
    """
    Return a dict for the CSQ row that best matches the current alt(s).
    - Prefer entries whose Allele matches the first ALT.
    - Among those, prefer one that has ClinVar_ALLELEID (if prefer_clinvar=True).
    Fallback: any CSQ row.
    """
    csq = var.INFO.get("CSQ")
    if not csq:
        return None

    alts = [str(a) for a in (var.ALT or [])]
    alt0 = alts[0] if alts else None

    rows = []
    for entry in csq.split(","):
        cols = entry.split("|")
        d = dict(zip(csq_format, cols + [""]*(len(csq_format)-len(cols))))
        rows.append(d)

    # filter by allele match if possible
    if alt0:
        rows_alt = [d for d in rows if d.get("Allele") == alt0]
    else:
        rows_alt = rows[:]

    if prefer_clinvar:
        rows_cv = [d for d in rows_alt if (d.get("ClinVar"))]
        if rows_cv:
            return rows_cv[0]

    return (rows_alt[0] if rows_alt else rows[0]) if rows else None

def _num_or_none(x):
    try:
        if x is None or x == ".":
            return None
        return float(x)
    except Exception:
        return None

def parse_spliceai(val):
    """Return the max SpliceAI DS value and its event name."""
    if not val:
        return None, None
    best_ds, best_event = None, None
    for entry in str(val).split(","):
        parts = entry.split("|")
        if len(parts) < 10:
            continue
        try:
            ds_ag, ds_al, ds_dg, ds_dl = map(float, parts[2:6])
        except ValueError:
            continue
        ds_map = {
            "DS_AG": ds_ag,
            "DS_AL": ds_al,
            "DS_DG": ds_dg,
            "DS_DL": ds_dl,
        }
        event, ds = max(ds_map.items(), key=lambda kv: kv[1])
        if best_ds is None or ds > best_ds:
            best_ds, best_event = ds, event
    return best_ds, best_event

# -------------------------
# Coverage summary loader (your format)
# -------------------------
coverage_data = {}
with open(args.coverage_summary) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 5: continue
        region = parts[0]             # e.g. "1:17018000-17019000"
        cov = {}
        for p in parts[1:]:
            key, val = p.split(":")
            cov[key.strip()] = val.strip().replace("%","")
        coverage_data[region] = cov

def get_exon_cov(chrom, pos):
    for region, cov in coverage_data.items():
        rchrom, coords = region.split(":")
        start, end = map(int, coords.split("-"))
        if str(chrom) == rchrom and start <= pos <= end:
            return cov
    return {">=20x":"NA", ">=30x":"NA", ">=50x":"NA", ">=100x":"NA"}

# R1/R2 and FWD/REV
r1r2_data = {}
with open(args.r1r2_summary) as f:
    for line in f:
        chrom, start, end, gene, r1, r2, ratio = line.strip().split("\t")
        r1r2_data[(chrom, int(start), int(end))] = (r1, r2, gene, ratio)

def get_r1r2(chrom, pos):
    for (c,s,e), vals in r1r2_data.items():
        if str(chrom)==c and s<=pos<=e: return vals
    return ("NA","NA","NA")

fr_data = {}
with open(args.fr_summary) as f:
    for line in f:
        chrom, start, end, gene, fwd, rev, ffrac, rfrac = line.strip().split("\t")
        fr_data[(chrom, int(start), int(end))] = (gene, fwd, rev, ffrac, rfrac)

def get_fr(chrom, pos):
    for (c,s,e), vals in fr_data.items():
        if str(chrom)==c and s<=pos<=e: return vals
    return ("NA","NA","NA","NA")

# Optional SF gene list
sf_genes = set()
if args.sf_genes and os.path.exists(args.sf_genes):
    with open(args.sf_genes) as f:
        for l in f:
            g = l.strip().split()[0]
            if g: sf_genes.add(g)

# Optional HemOnc gene list (BioVarFlow_HemOnc branch)
hemonc_genes = set()
if args.hemonc_genes and os.path.exists(args.hemonc_genes):
    with open(args.hemonc_genes) as f:
        for l in f:
            g = l.strip().split()[0] if l.strip() else ""
            if g and not g.startswith("#"): hemonc_genes.add(g)

# -------------------------
# Parse VCF to dataframe
# -------------------------
vcf = VCF(args.vcf)
sample_id = args.sample_id or sample_from_vcf(args.vcf)

# CSQ fields (if present)
csq_format = []
for h in vcf.raw_header.split("\n"):
    if h.startswith("##INFO=<ID=CSQ"):
        try:
            csq_format = h.split("Format: ")[1].rstrip('">').split("|")
        except: pass
        break

records = []
for var in vcf:
    chrom, pos, ref = var.CHROM, var.POS, var.REF
    alt = ",".join(var.ALT) if var.ALT else ""
    filt = var.FILTER or "PASS"
    qual = _num_or_none(getattr(var, "QUAL", None))
    dp = var.format("DP")[0][0] if "DP" in var.FORMAT else var.INFO.get("DP")
    vaf = var.format("VAF")[0][0] if "VAF" in var.FORMAT else None
    gq  = var.format("GQ")[0][0] if "GQ" in var.FORMAT else None

    gt = var.genotypes[0][:2]
    zyg = "het" if gt[0] != gt[1] else ("hom" if gt[0] == 1 else "ref")

    ad_ref = ad_alt = None
    if "AD" in var.FORMAT:
        try: ad_ref, ad_alt = var.format("AD")[0]
        except: pass

    gene = transcript = hgvsc = hgvsp = consequence = exon = intron = impact = None
    clinvar = alleleid = gnomad_af = revel = spliceai_fields = spliceai_event = spliceai_ds = cadd = mane_id = bayesdel_score = clinvar_review_status = stars = am_pathogenicity = am_class = None
    ann = {}
    if csq_format:
        ann = select_csq_entry(var, csq_format) or {}
        gene        = ann.get("SYMBOL")
        transcript  = ann.get("Feature")
        hgvsc       = ann.get("HGVSc")
        hgvsp       = ann.get("HGVSp")
        consequence = ann.get("Consequence")
        mane_id = ann.get("MANE_SELECT")
        exon        = ann.get("EXON")
        intron      = ann.get("INTRON")
        impact      = ann.get("IMPACT")
        am_pathogenicity = ann.get("am_pathogenicity")
        am_class    = ann.get("am_class")
        gnomad_af   = ann.get("gnomADg_AF") or ann.get("gnomADe_AF") or ann.get("AF")
        revel       = ann.get("REVEL")
        spliceai_fields = {
        "DS_AG": ann.get("SpliceAI_pred_DS_AG"),
        "DS_AL": ann.get("SpliceAI_pred_DS_AL"),
        "DS_DG": ann.get("SpliceAI_pred_DS_DG"),
        "DS_DL": ann.get("SpliceAI_pred_DS_DL"),
        }
        bayesdel_score = ann.get("BayesDel")
        clinvar     = ann.get("ClinVar_CLNSIG") or ann.get("CLIN_SIG")
        alleleid    = ann.get("ClinVar")
        clinvar_review_status = ann.get("ClinVar_CLNREVSTAT") or ann.get("CLNREVSTAT")
        cadd        = ann.get("CADD_PHRED") or ann.get("CADDraw")
    
    cov = get_exon_cov(chrom, pos)
    r1r2 = get_r1r2(chrom, pos)
    fr = get_fr(chrom, pos)
    revstat = clinvar_review_status if clinvar_review_status is not None else ""
    stars = clinvar_stars_from_revstat(revstat)

    spliceai_scores = {k: float(v) for k, v in spliceai_fields.items() 
                      if v not in [None, ""]}
    spliceai_event, spliceai_ds = (
        max(spliceai_scores.items(), key=lambda kv: kv[1]) 
        if spliceai_scores else (None, None)
    )

    try:
        r1 = int(r1r2[0]); r2 = int(r1r2[1])
        r1r2_abs = round(r1 / r2, 2) if r2 > 0 else None
    except: r1r2_abs = None
    try:
        fwd = int(fr[0]); rev = int(fr[1])
        fr_abs = round(fwd / rev, 2) if rev > 0 else None
    except: fr_abs = None

    rec = {
        "Sample": sample_id,
        "Chrom": chrom, "Pos": pos, "Ref": ref, "Alt": alt,
        "Variant": f"{chrom}:{pos}:{ref}:{alt}",
        "FILTER": filt,
        "QUAL": qual,
        "Gene": gene, "Transcript": transcript,
        "Consequence": consequence, "Impact": impact,
        "Exon": exon, "Intron": intron,
        "HGVSc": hgvsc, "HGVSp": hgvsp,
        "MANE_ID": mane_id,
        "GT": f"{gt[0]}/{gt[1]}",
        "AD_Ref": ad_ref, "AD_Alt": ad_alt,
        "DP": dp, "GQ": gq, "VAF": vaf, "Zygosity": zyg,
        "ExonCov20": cov.get(">=20x","NA"),
        "ExonCov30": cov.get(">=30x","NA"),
        "ExonCov50": cov.get(">=50x","NA"),
        "ExonCov100": cov.get(">=100x","NA"),
        "R1": r1r2[0], "R2": r1r2[1], "R1R2_frac": r1r2[2], "R1R2_abs": r1r2_abs,
        "FWD": fr[0], "REV": fr[1], "FWD_frac": fr[2], "REV_frac": fr[3], "FR_abs": fr_abs,
        "gnomAD_AF": gnomad_af, "REVEL": revel, "SpliceAI_DS_max": spliceai_ds, "SpliceAI_Event": spliceai_event, "BayesDel_score": bayesdel_score, "CADD": cadd, "AM_Pathogenicity": am_pathogenicity, "AM_Class": am_class,
        "ClinVar": clinvar, "ClinVar_ReviewStatus": revstat, "ClinVar_Stars": stars, "ALLELEID": alleleid, "ClinVar_ALLELEID": alleleid,  
        "ClinVar_StarsGlyph" : clinvar_star_glyph(stars),
        "HGVS_full": None
    }

    # Combined HGVS (Gene c. p.)
    if hgvsc:
        rec["HGVSc"] = hgvsc.split(":")[1] if ":" in hgvsc else hgvsc
    if hgvsp:
        rec["HGVSp"] = hgvsp.split(":")[1] if ":" in hgvsp else hgvsp
    if gene or hgvsc or hgvsp:
        rec["HGVS_full"] = f"{gene or ''} {hgvsc or ''} {hgvsp or ''}".strip()
    records.append(rec)

df = pd.DataFrame(records)

# -------------------------
# Build tabs
# -------------------------
with pd.ExcelWriter(args.xlsx_out) as xw:

    # 1) Sample Summary (single row)
    ss = {
        "Sample_ID": sample_id,
        "Assay": args.assay,
        "Build": args.build
    }
    ss.update(parse_flagstat(args.flagstat))
    ss.update(parse_samtools_stats(args.stats))
    # Insert size & mapping/Q30 from Picard if present
    tmp = parse_picard_alignment(args.picard_align);  ss.update(tmp)
    tmp = parse_picard_insert(args.picard_insert);    ss.update({k:v for k,v in tmp.items() if v is not None})
    # Mosdepth means
    ss.update(parse_mosdepth_summary_regions(args.mosdepth_summary))
    # Scope the ACMG region-coverage metric to the ACMG SF gene set (so it isn't
    # diluted by HemOnc regions when using a merged ACMG+HemOnc BED). Fall back
    # to the legacy "all regions" behaviour when no SF list is provided.
    ss.update(acmg_pct_regions_covered(args.acmg_thresholds,
                                       gene_set=sf_genes if sf_genes else None,
                                       key_prefix="ACMG"))
    # Parallel HemOnc-scoped metric (only emitted when a HemOnc list is provided).
    if hemonc_genes:
        ss.update(acmg_pct_regions_covered(args.acmg_thresholds,
                                           gene_set=hemonc_genes,
                                           key_prefix="HemOnc"))
    # Contamination & sex
    ss.update(parse_verifybamid(args.verifybamid))
    ss.update(parse_sexcheck(args.sexcheck))

    # Derive Ti/Tv if not provided (rough, SNP only)
    try:
        snps = df[(df["Ref"].str.len()==1) & (df["Alt"].str.len()==1)]
        transitions = {("A","G"),("G","A"),("C","T"),("T","C")}
        ti = sum((r.Ref,r.Alt) in transitions for _,r in snps.iterrows())
        tv = len(snps)-ti
        if tv>0 and "TiTv" not in ss: ss["TiTv"] = round(ti/tv,3)
    except: pass

    # Het/Hom summary from VCF if not in stats
    try:
        het = (df["Zygosity"]=="het").sum()
        hom = (df["Zygosity"]=="hom").sum()
        if hom>0 and "Het_Hom_Ratio" not in ss: ss["Het_Hom_Ratio"] = round(het/hom,3)
    except: pass

    pd.DataFrame([ss]).to_excel(xw, index=False, sheet_name="Sample Summary")

    # ---------------------------------------------------------------
    # Helper: write the three per-gene-set sheets (variants, coverage
    # gaps, per-gene rollup). Shared by the ACMG SF and HemOnc blocks.
    # ---------------------------------------------------------------
    def make_hyperlink(u):
        return f'=HYPERLINK("{u}","{u}")' if u else ""

    variant_cols = [
        "Gene","Variant","HGVSc","HGVSp","MANE_ID","Zygosity","GT","AD_Ref","AD_Alt","DP","GQ","QUAL",
        "Consequence","Exon","Intron","ClinVar","ClinVar_ReviewStatus","ClinVar_Stars","ClinVar_StarsGlyph","ClinVar_Link","gnomAD_AF","REVEL","SpliceAI_DS_max","SpliceAI_Event","BayesDel_score","AM_Pathogenicity","AM_Class","HGVS_full"
    ]
    empty_gaps_cols = [
        "Gene","MANE_ID","Exon","RegionLabel","Chrom","ExonStart","ExonEnd","ExonLen",
        "Pct>=20x","Pct>=30x",
        "Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
        "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"
    ]

    def write_gene_set_sheets(gene_set, variants_sheet, gaps_sheet, genes_cov_sheet):
        """
        Emit three sheets to the open ExcelWriter `xw`, scoped to `gene_set`:
          1. `variants_sheet`   — ClinVar-pathogenic variants in genes of `gene_set`
          2. `gaps_sheet`       — per-exon coverage + gap intervals for genes of `gene_set`
          3. `genes_cov_sheet`  — per-gene length-weighted %≥20x rollup

        An empty `gene_set` produces the legacy "no gene filter" behaviour for
        the variants sheet (kept for backwards compatibility with existing
        ACMG runs that omit --sf-genes). Coverage-gap sheets are always
        gene-scoped when a set is provided (they can't sensibly be "all
        genes" because thresholds/gaps files come from the merged BED).
        """
        # ---- variants sheet ----
        variants = df.copy()
        variants["ClinVar_URL"] = variants.apply(clinvar_url_from_row, axis=1)
        variants["ClinVar_Link"] = variants["ClinVar_URL"].map(make_hyperlink)
        if "ClinVar" in variants.columns:
            variants = variants[variants["ClinVar"].astype(str).str.contains("pathogenic", case=False, na=False)]
            #variants = variants[variants["ClinVar_Stars"].fillna(0).astype(int) >= 2]
        else:
            variants = variants.iloc[0:0]
        if gene_set:
            variants = variants[variants["Gene"].isin(gene_set)]
        variants[variant_cols].to_excel(xw, index=False, sheet_name=variants_sheet)

        # ---- coverage gaps sheet ----
        combined = pd.DataFrame(columns=empty_gaps_cols)
        try:
            th_exon = read_acmg_thresholds(args.acmg_thresholds, sf_genes=gene_set)  # gene-filtered
            g20_agg = load_gap_bed_for_agg(args.gaps20, "<20x", sf_genes=gene_set)
            g30_agg = load_gap_bed_for_agg(args.gaps30, "<30x", sf_genes=gene_set)

            combined = (th_exon
                        .merge(g20_agg, on="RegionLabel", how="left")
                        .merge(g30_agg, on="RegionLabel", how="left"))

            for c in ["Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
                      "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"]:
                if c not in combined.columns:
                    combined[c] = pd.NA

            combined = (combined[
                ["Gene","MANE_ID","Exon","Chrom","ExonStart","ExonEnd","ExonLen",
                 "Pct>=20x","Pct>=30x",
                 "Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
                 "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"]]
                .sort_values(["Chrom","ExonStart","ExonEnd"])
            )
        except Exception:
            pass

        combined.to_excel(xw, index=False, sheet_name=gaps_sheet)

        # ---- per-gene rollup sheet ----
        try:
            genes_cov = build_genes_coverage(combined, sf_genes=gene_set if gene_set else None)
        except Exception:
            genes_cov = pd.DataFrame(columns=["Gene","MANE_ID","Pct20"])

        # Ensure every gene in `gene_set` appears in the rollup — genes with
        # no coverage data (e.g. because the mosdepth BED used at run time
        # didn't include them) show up with Pct20 = NaN. This keeps the sheet
        # honest about scope so downstream consumers (HTML report grid) can
        # render "NA" for unmeasured genes rather than silently omitting them.
        if gene_set:
            all_genes = pd.DataFrame({"Gene": sorted(gene_set)})
            genes_cov = (all_genes.merge(genes_cov, on="Gene", how="left")
                                  .fillna({"MANE_ID": "NA"}))
            # Pct20 stays NaN for genes with no data (excluded from fillna on
            # purpose) so downstream code can distinguish "unmeasured" from
            # "measured and 0%".

        genes_cov.to_excel(xw, index=False, sheet_name=genes_cov_sheet)

    # 2-4) ACMG SF sheets (existing behaviour, sheet names preserved for
    # backwards compatibility with the HTML report reader).
    write_gene_set_sheets(
        gene_set=sf_genes,
        variants_sheet="ACMG SF (P-LP)",
        gaps_sheet="Coverage gaps",
        genes_cov_sheet="ACMG SF Genes Coverage",
    )

    # 4b) HemOnc sheets (BioVarFlow_HemOnc branch). Only written when a
    # HemOnc gene list was supplied via --hemonc-genes.
    if hemonc_genes:
        write_gene_set_sheets(
            gene_set=hemonc_genes,
            variants_sheet="HemOnc (P-LP)",
            gaps_sheet="HemOnc Coverage gaps",
            genes_cov_sheet="HemOnc Genes Coverage",
        )

    # 5) PASS variant table
    pass_cols = [
        "Variant","Gene","HGVSc","HGVSp","MANE_ID","Transcript","Consequence","GT","Zygosity","AD_Ref","AD_Alt","DP","GQ","QUAL",
        "FILTER","VAF","gnomAD_AF","ClinVar","ClinVar_ReviewStatus","ClinVar_Stars","ClinVar_StarsGlyph","REVEL","SpliceAI_DS_max","SpliceAI_Event","BayesDel_score","AM_Pathogenicity","AM_Class","HGVS_full"
    ]
    df_pass = df[df["FILTER"]=="PASS"].copy()
    df_pass[pass_cols].to_excel(xw, index=False, sheet_name="PASS variants")

    # 6) PASS variants filtered to gnomAD AF < 1% (or absent from gnomAD)
    #    and/or ClinVar Pathogenic/Likely pathogenic
    af_numeric = df_pass["gnomAD_AF"].map(_num_or_none)
    is_rare = af_numeric.isna() | (af_numeric < 0.01)
    is_plp = df_pass["ClinVar"].astype(str).str.contains("pathogenic", case=False, na=False)
    df_pass[is_rare | is_plp][pass_cols].to_excel(xw, index=False, sheet_name="PASS Rare or P-LP")

print(f"Wrote Excel report → {args.xlsx_out}")
