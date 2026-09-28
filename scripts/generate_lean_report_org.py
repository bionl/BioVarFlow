#!/usr/bin/env python3
import sys, os, re, argparse
import pandas as pd
from cyvcf2 import VCF
from urllib.parse import quote_plus, unquote

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
p.add_argument("--acmg-thresholds", default=None,
               help="mosdepth thresholds.bed(.gz) run with ACMG BED (cols: chrom start end [region] 20X 30X)")
p.add_argument("--gaps20", default=None, help="annotated gaps <20x BED (chrom start end RegionLabel)")
p.add_argument("--gaps30", default=None, help="annotated gaps <30x BED (chrom start end RegionLabel)")
p.add_argument("--mane", default=None,
               help="TSV mapping Ensembl transcripts to RefSeq (ENST<tab>NM_<tab>symbol), from "
                    "the NCBI MANE summary. Off-panel variants never reach VEP, so Exomiser's "
                    "Ensembl accession is all we have for them; interpretation is done against "
                    "RefSeq.")
p.add_argument("--sites-vcf", default=None,
               help="StrandBiasPileup's <sample>_sb_sites.vcf.gz -- the unfiltered consensus "
                    "subset to panel + Exomiser positions. The only source of FORMAT/AD for "
                    "OFF-panel variants, which never reach PASS variants because BedFilterVCF "
                    "removes them before annotation. Used to give the Exomiser and Prioritised "
                    "tabs a VAF.")
p.add_argument("--exomiser", default=None,
               help="Exomiser <sample>_exomiser.variants.tsv. Adds an 'Exomiser' tab and a "
                    "'Prioritised' tab. A zero-byte NO_FILE placeholder means Exomiser did not "
                    "run for this sample and both tabs are skipped.")
p.add_argument("--exomiser-top", type=int, default=0,
               help="Rows kept on the Exomiser tab, by rank. 0 (the default) keeps every variant "
                    "Exomiser emitted -- 745 for IQMM. Set a positive number to cap the tab; the "
                    "full TSV and HTML are published under <outdir>/exomiser either way.")
p.add_argument("--panel-genes", default=None,
               help="Gene symbols on the capture panel, one per line. Used ONLY to split the "
                    "Prioritised tab into PANEL_VUS and OFF_PANEL -- unlike --sf-genes it does "
                    "not filter any sheet. Without it the panel is inferred from the genes "
                    "present in this sample's report, which misses panel genes the sample has "
                    "no calls in.")
p.add_argument("--prioritised-min-score", type=float, default=0.05,
               help="Prioritised tab: minimum Exomiser combined score for an OFF-panel variant "
                    "to be listed (default 0.05). Panel variants ranked by Exomiser are always "
                    "listed regardless of score. Set very high to make the tab a pure "
                    "panel/Exomiser intersection.")
p.add_argument("--strand-bias", default=None,
               help="strand_bias.py TSV: per-variant ALT_FWD/ALT_REV and Fisher p from the raw alignment")



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

    # RegionLabel is the JOIN KEY against the gaps BED, so it must be unique per
    # exon. The BED's 4th column is GENE|TRANSCRIPT and repeats across every exon
    # of a gene -- keying on it collapsed gaps to gene level and stamped one
    # gene-wide total onto every exon row. Key on exon coordinates instead;
    # CoverageGapsAnnotation labels gaps the same way.
    df["RegionLabel"] = (df[chrom].astype(str) + ":" +
                         df[start].astype(int).astype(str) + "-" +
                         df[end].astype(int).astype(str))

    # Gene / transcript still come from the BED's label column when present.
    if (region is not None) and (region != 'NA') and (region in df.columns):
        lbl = df[region].astype(str)
        lbl = lbl.where(~lbl.isin(["", ".", "unknown", "UNKNOWN", "None"]), other=".")
        parts = lbl.str.split("|", n=2, expand=True)
        df["Gene"]    = parts[0] if parts.shape[1] > 0 else "."
        df["MANE_ID"] = parts[1] if parts.shape[1] > 1 else "."
    else:
        df["Gene"] = "."
        df["MANE_ID"] = "."

    # Exon number: the BED carries no exon index, so derive one by ordering each
    # gene's intervals along the genome. Previously this column was always "."
    # because splitting GENE|TRANSCRIPT never yielded a third field.
    _ord = df.assign(_s=pd.to_numeric(df[start], errors="coerce"))
    df["Exon"] = (_ord.groupby("Gene")["_s"].rank(method="dense").astype("Int64").astype(str))
    df.loc[df["Gene"] == ".", "Exon"] = "."

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

def acmg_pct_regions_covered(th_path, min_depth=20):
    """
    Compute % of ACMG regions fully covered at >= min_depth using mosdepth thresholds output.
    Expects columns like: #chrom, start, end, region, 10X, 20X, 30X, ...
    Returns: {"ACMG_PctRegionsCovered": 97.35, "ACMG_RegionsCovered": N, "ACMG_RegionsTotal": M}
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

    if chrom is None or start is None or end is None:
        return out

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

    # 6) Each row is one exon interval — use coordinates as the region id so we
    #    evaluate per-exon coverage instead of collapsing all exons of a gene.
    #    (Grouping by the BED's 4th-column label would collapse them: the labels
    #    are GENE|TRANSCRIPT and repeat across every exon of a gene, so a gene
    #    would only count as covered if 100% of ALL its bases were >= min_depth,
    #    which almost nothing passes even at high mean coverage.)
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

    out["ACMG_PctRegionsCovered"] = pct
    out["ACMG_RegionsCovered"] = covered
    out["ACMG_RegionsTotal"] = total_regions
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
    if not s: return ""
    # Lowercase, replace underscores with spaces, collapse spaces/commas/semicolons
    s = s.lower().replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s

EXOMISER_COLS = [
    ("#RANK",                          "Rank"),
    ("GENE_SYMBOL",                    "Gene"),
    ("HGVS",                           "HGVS"),
    ("GENOTYPE",                       "Genotype"),
    ("FUNCTIONAL_CLASS",               "Consequence"),
    ("MOI",                            "MOI"),
    ("EXOMISER_GENE_COMBINED_SCORE",   "Combined_Score"),
    ("EXOMISER_GENE_PHENO_SCORE",      "Pheno_Score"),
    ("EXOMISER_VARIANT_SCORE",         "Variant_Score"),
    ("CONTRIBUTING_VARIANT",           "Contributing"),
    ("EXOMISER_ACMG_CLASSIFICATION",   "Exomiser_ACMG"),
    ("EXOMISER_ACMG_DISEASE_NAME",     "Exomiser_Disease"),
    ("CLINVAR_PRIMARY_INTERPRETATION", "ClinVar"),
    ("CLINVAR_STAR_RATING",            "ClinVar_Stars"),
    ("MAX_FREQ",                       "Max_Freq_pct"),
    ("GENE_CONSTRAINT_LOEUF",          "LOEUF"),
]


def load_mane(path):
    """{ENST without version: (RefSeq_nuc, symbol)} from the MANE summary."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0].startswith("ENST"):
                out[f[0].split(".")[0]] = (f[1], f[2])
    return out


def split_exomiser_hgvs(value):
    """Exomiser's GENE:TRANSCRIPT:c.:p. field -> (transcript, c., p.)."""
    tx = c = pr = None
    for part in str(value or "").split(":"):
        if part.startswith("ENST") or part.startswith("NM_"):
            tx = part
        elif part.startswith("c."):
            c = part
        elif part.startswith("p."):
            pr = part
    return tx, c, pr


def to_refseq(tx, gene, mane):
    """ENST -> MANE Select RefSeq, only when the gene symbol agrees.

    A silent swap to the wrong transcript is worse than an ENST a reviewer can
    look up, so a symbol mismatch leaves the accession untouched.
    """
    if not tx or not str(tx).startswith("ENST"):
        return tx
    hit = mane.get(str(tx).split(".")[0])
    return hit[0] if hit and (not gene or str(gene) == hit[1]) else tx


def load_sites_vaf(path):
    """{(chrom,pos,ref,alt): VAF} from the strand-bias sites VCF.

    Computed against SITE-level depth, not per record. The consensus is already
    split by `norm -m -any`, so each record carries AD=[ref, this_alt] and
    alt/(ref+alt) would drop the sibling allele at a 1/2 site -- the same error
    the multiallelic reconciliation pass exists to correct, and the one that let
    18 sub-20% rows through the VAF floor. Records are therefore grouped by
    position and the denominator is AD_ref + every ALT's AD.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    recs = collections.defaultdict(list)
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                keys = f[8].split(":")
                if "AD" not in keys:
                    continue
                try:
                    ad = [int(x) for x in f[9].split(":")[keys.index("AD")].split(",")]
                except ValueError:
                    continue
                if len(ad) < 2:
                    continue
                recs[(f[0], f[1])].append((f[3].upper(), f[4].upper(), ad))
    except OSError:
        return {}

    out = {}
    for (chrom, pos), rows in recs.items():
        site = rows[0][2][0] + sum(r[2][1] for r in rows)   # AD_ref + all ALTs
        if site <= 0:
            continue
        for ref, alt, ad in rows:
            out[(chrom, pos, ref, alt)] = round(ad[1] / site, 4)
    return out


def load_strand_bias(path):
    """{(chrom, pos, ref, alt): row} from strand_bias.py's TSV."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    try:
        sb = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return {}
    if "StrandBias_Flag" not in sb.columns:
        return {}
    sb = sb.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
    return {(r["CHROM"], r["POS"], r["REF"], r["ALT"]): r for _, r in sb.iterrows()}


def attach_strand_bias(view, sb, vaf=None):
    """Add the strand-bias columns to an Exomiser view, keyed on Variant.

    The check now runs on the RAW BAM over the union of panel and
    Exomiser-ranked sites, so off-panel rows get a real result instead of a
    blank. Before that they were untestable: the panel-filtered BAM does not
    contain their reads, and an empty cell reads as "clean" to a reviewer --
    which is exactly the failure the status labels exist to prevent.
    """
    vaf = vaf or {}
    cols = {"VAF": [], "ALT_FWD": [], "ALT_REV": [], "StrandBias_P": [], "StrandBias": []}
    for v in view["Variant"].astype(str):
        parts = v.split(":")
        key = tuple(parts) if len(parts) == 4 else None
        # VAF for these tabs comes from the sites VCF, which covers off-panel
        # variants too; PASS variants holds only the panel ones.
        cols["VAF"].append(vaf.get(key, pd.NA))
        row = sb.get(key) if key else None
        if row is None:
            cols["ALT_FWD"].append(pd.NA); cols["ALT_REV"].append(pd.NA)
            cols["StrandBias_P"].append(pd.NA)
            cols["StrandBias"].append("NOT_TESTED_NO_PILEUP" if sb else "NOT_RUN")
        else:
            cols["ALT_FWD"].append(row["ALT_FWD"]); cols["ALT_REV"].append(row["ALT_REV"])
            cols["StrandBias_P"].append(row["StrandBias_P"])
            cols["StrandBias"].append(row["StrandBias_Flag"])
    out = view.copy()
    for k, v in cols.items():
        out[k] = v
    return out


def load_exomiser(path):
    """Exomiser variants TSV -> DataFrame, or None when it did not run.

    LeanReport is handed a zero-byte assets/NO_FILE placeholder for samples
    Exomiser skipped, so that the workbook is still produced. Treat any empty
    or header-less file the same way.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return None
    if "#RANK" not in df.columns or df.empty:
        return None
    df["Variant"] = ("chr" + df["CONTIG"].astype(str).str.removeprefix("chr")
                     + ":" + df["START"].astype(str)
                     + ":" + df["REF"].astype(str) + ":" + df["ALT"].astype(str))
    df["_rank"] = pd.to_numeric(df["#RANK"], errors="coerce")
    return df.sort_values("_rank", kind="stable")


def exomiser_view(df):
    """Rename to the reviewer-facing columns, keeping only those present."""
    cols = [(src, dst) for src, dst in EXOMISER_COLS if src in df.columns]
    out = df[["Variant"] + [c for c, _ in cols]].copy()
    out = out.rename(columns=dict(cols))
    return out[["Rank", "Gene", "Variant"] +
               [d for _, d in cols if d not in ("Rank", "Gene")]]


def is_pathogenic_clinvar(clinvar_val):
    """
    True only for an exact Pathogenic / Likely_pathogenic classification
    (including combined "Pathogenic/Likely_pathogenic" or "Pathogenic&Likely_pathogenic").
    Rejects substring-adjacent but distinct terms like
    "Conflicting_classifications_of_pathogenicity" or "risk_factor".
    """
    if not clinvar_val:
        return False
    exact = {"pathogenic", "likely pathogenic"}
    for tok in re.split(r"[&/|,]", str(clinvar_val)):
        if tok.strip().lower().replace("_", " ") in exact:
            return True
    return False


def clinvar_stars_from_revstat(revstat: str) -> int:
    """
    Map ClinVar review status text to 0–4 stars (germline).
    Handles multiple statuses separated by | , ; — returns the maximum.
    """
    if not revstat: return None
    toks = re.split(r"[|,;]+", revstat)
    def map_one(t: str) -> int:
        t = _norm_revstat(t).strip()
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
    print(alleleid)
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
        # VEP percent-encodes reserved characters inside CSQ subfields, so a
        # synonymous change arrives as "p.Gln342%3D" instead of "p.Gln342=".
        # Decode here so both the Excel and the HTML report show valid HGVS.
        hgvsc       = unquote(ann["HGVSc"]) if ann.get("HGVSc") else ann.get("HGVSc")
        hgvsp       = unquote(ann["HGVSp"]) if ann.get("HGVSp") else ann.get("HGVSp")
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


# ---------------------------------------------------------------------------
# Multiallelic site reconciliation
# ---------------------------------------------------------------------------
# `bcftools norm -m -any` splits a multiallelic site into one biallelic record
# per ALT. That is required for annotation, but it loses site-level context in
# ways that misrepresent the variant:
#
#   * GT becomes 1/0 and 0/1 (a biallelic record cannot express 1/2), so
#     Zygosity reads "het" for a site carrying no reference allele at all.
#   * FORMAT/AD is subset to [ref, this_alt], so the OTHER alt's reads vanish
#     from the VAF denominator. bcftools computes VAF as alt/(ref+alt), which
#     for a 1/2 site inflates every allele -- a site with AD 0,24,17 and DP 41
#     reported VAF 1.000 for BOTH alleles instead of 0.585 / 0.415.
#
# Rebuild the site by grouping the split records back together on coordinates,
# then recompute VAF against total allele-informative depth at the site
# (AD_Ref + sum of every ALT's AD) and restate the genotype.
if len(df) and {"Chrom", "Pos", "AD_Ref", "AD_Alt"} <= set(df.columns):
    _key = ["Chrom", "Pos"]
    _ar = pd.to_numeric(df["AD_Ref"], errors="coerce")
    _aa = pd.to_numeric(df["AD_Alt"], errors="coerce")

    _n_alt   = df.groupby(_key)["AD_Alt"].transform("size")
    _alt_sum = _aa.groupby([df[c] for c in _key]).transform("sum")
    # AD_Ref is the site's reference depth, repeated on each split row -- take
    # it once, not once per ALT.
    _site_ad = _ar + _alt_sum

    df["Multiallelic"] = (_n_alt > 1).map({True: "Y", False: "N"})
    df["Site_AD_Total"] = _site_ad.astype("Int64")

    # VAF against site-level informative depth. Single-ALT rows are unchanged
    # (ref + alt is already the site total).
    df["VAF"] = (_aa / _site_ad).where(_site_ad > 0).round(4)

    # Restate the genotype for split rows. A 1/2 site is heterozygous for two
    # different ALT alleles; "het" alone hides that there is no ref allele.
    _multi = _n_alt > 1
    df["GT_site"] = df["GT"]
    df.loc[_multi, "GT_site"] = "1/2"
    df.loc[_multi, "Zygosity"] = "het (1/2)"

    # ---- QC flags -------------------------------------------------------
    # AD counts allele-assigned reads; DP counts all reads at the locus, so a
    # gap between them is normal (uninformative/filtered reads). Surface it
    # rather than silently absorbing it into whichever denominator we chose.
    _dp = pd.to_numeric(df["DP"], errors="coerce")
    df["AD_Sum_vs_DP"] = (_site_ad - _dp).astype("Int64")

    # A hom-ALT call retaining substantial reference support is worth a look:
    # it can indicate a mis-genotype, CNV, or contamination.
    _ref_frac = (_ar / _site_ad).where(_site_ad > 0)
    _hom_with_ref = df["Zygosity"].astype(str).eq("hom") & (_ref_frac > 0.10)

    # ---- Strand bias (raw alignment) -----------------------------------
    # Joined on coordinates + alleles from strand_bias.py, which recomputes the
    # per-allele forward/reverse counts with bcftools mpileup. The caller's own
    # FS is computed after local reassembly and can badly understate the skew:
    # MSH2 chr2:47414420 T>G had FS=43.4 (passing FS>60) but 0/14 forward alt
    # reads in the raw pileup, p=4.9e-08. Advisory only -- flagged, not dropped.
    #
    # StrandBias carries the OUTCOME in words, because three very different
    # situations otherwise render identically in Excel as an empty cell or a
    # high p-value:
    #   PASS                    tested with enough reads to have found skew
    #   UNDERPOWERED            too few alt reads for any split to reach the
    #                           threshold -- absence of evidence, not evidence
    #                           of absence (SEC23B chr20:18524919, 9 alt reads,
    #                           best achievable p = 0.01)
    #   NOT_TESTED_*            indel, or no usable pileup record
    # A reviewer reading a blank cell as "no problem" is exactly the failure
    # this check exists to prevent, so the status is never left implicit.
    df["ALT_FWD"] = pd.NA
    df["ALT_REV"] = pd.NA
    df["StrandBias_P"] = pd.NA
    # NOT_RUN distinguishes "the check was never invoked" from "invoked and
    # found nothing for this variant".
    df["StrandBias"] = "NOT_RUN"
    _sb_flag = pd.Series(False, index=df.index)
    if args.strand_bias and os.path.exists(args.strand_bias):
        df["StrandBias"] = "NOT_TESTED_NO_PILEUP"
        _sb = pd.read_csv(args.strand_bias, sep="\t", dtype=str)
        _sb = _sb.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
        _sb_idx = _sb.set_index(["CHROM", "POS", "REF", "ALT"])
        _keys = list(zip(df["Chrom"].astype(str), df["Pos"].astype(str),
                         df["Ref"].astype(str), df["Alt"].astype(str)))
        for _i, _k in zip(df.index, _keys):
            if _k not in _sb_idx.index:
                continue
            _row = _sb_idx.loc[_k]
            df.at[_i, "ALT_FWD"] = _row["ALT_FWD"]
            df.at[_i, "ALT_REV"] = _row["ALT_REV"]
            df.at[_i, "StrandBias_P"] = _row["StrandBias_P"]
            df.at[_i, "StrandBias"] = _row["StrandBias_Flag"]
            _sb_flag.at[_i] = (_row["StrandBias_Flag"] == "STRAND_BIAS")

    _flags = []
    for _i in df.index:
        _f = []
        if _hom_with_ref.get(_i, False):
            _f.append(f"HOM_REF_SUPPORT={_ref_frac[_i]:.1%}")
        if pd.notna(df.at[_i, "AD_Sum_vs_DP"]) and df.at[_i, "AD_Sum_vs_DP"] != 0:
            _f.append(f"AD_DP_DIFF={df.at[_i, 'AD_Sum_vs_DP']}")
        if df.at[_i, "Multiallelic"] == "Y":
            _f.append("MULTIALLELIC")
        if _sb_flag.get(_i, False):
            _f.append(f"STRAND_BIAS(alt {df.at[_i, 'ALT_FWD']}+/{df.at[_i, 'ALT_REV']}-, "
                      f"p={df.at[_i, 'StrandBias_P']})")
        _flags.append(";".join(_f))
    df["QC_Flags"] = _flags

for _c in ("QC_Flags", "ALT_FWD", "ALT_REV", "StrandBias_P", "StrandBias", "Multiallelic"):
    if _c not in df.columns:
        df[_c] = pd.NA

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
    ss.update(acmg_pct_regions_covered(args.acmg_thresholds))
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

    # 2) ACMG SF (P/LP) — filter by ClinVar and (optionally) SF gene list
    
    acmg = df.copy()
    acmg["ClinVar_URL"] = acmg.apply(clinvar_url_from_row, axis=1)

    def make_hyperlink(u):
        return f'=HYPERLINK("{u}","{u}")' if u else ""

    acmg["ClinVar_Link"] = acmg["ClinVar_URL"].map(make_hyperlink)
    if "ClinVar" in acmg.columns:
        acmg = acmg[acmg["ClinVar"].astype(str).str.contains("pathogenic", case=False, na=False)]
    else:
        acmg = acmg.iloc[0:0]
    if sf_genes:
        acmg = acmg[acmg["Gene"].isin(sf_genes)]
    acmg_cols = [
        "Gene","Variant","HGVSc","HGVSp","MANE_ID","Zygosity","GT","AD_Ref","AD_Alt","DP","GQ","QUAL",
        "Consequence","Exon","Intron","ClinVar","ClinVar_ReviewStatus","ClinVar_Stars","ClinVar_StarsGlyph","ClinVar_Link","gnomAD_AF","REVEL","SpliceAI_DS_max","SpliceAI_Event","BayesDel_score","AM_Pathogenicity","AM_Class","HGVS_full",
        # Raw-alignment strand bias: advisory, never a filter. See strand_bias.py.
        "ALT_FWD","ALT_REV","StrandBias_P","StrandBias","QC_Flags"
    ]
    acmg[acmg_cols].to_excel(xw, index=False, sheet_name="ACMG SF (P-LP)")
    _reportable = acmg          # kept for the Prioritised sheet below

    # 3) Coverage gaps in ACMG SF genes (or all genes if no list)
    combined_df = pd.DataFrame(columns=[
        "Gene","MANE_ID","Exon","RegionLabel","Chrom","ExonStart","ExonEnd","ExonLen",
        "Pct>=20x","Pct>=30x",
        "Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
        "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"
    ])

    try:
        th_exon = read_acmg_thresholds(args.acmg_thresholds, sf_genes=sf_genes)  # RegionLabel = chrom:start-end if no names
        g20_agg = load_gap_bed_for_agg(args.gaps20, "<20x", sf_genes=sf_genes)
        g30_agg = load_gap_bed_for_agg(args.gaps30, "<30x", sf_genes=sf_genes)

        combined_df = (th_exon
                       .merge(g20_agg, on="RegionLabel", how="left")
                       .merge(g30_agg, on="RegionLabel", how="left"))

        # ensure all expected columns exist
        for c in ["Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals","Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"]:
            if c not in combined_df.columns: combined_df[c] = pd.NA

        combined_df = (combined_df[
            ["Gene","MANE_ID","Exon","Chrom","ExonStart","ExonEnd","ExonLen",
             "Pct>=20x","Pct>=30x",
             "Gaps<20x_n","Gaps<20x_bp","Gaps<20x_intervals",
             "Gaps<30x_n","Gaps<30x_bp","Gaps<30x_intervals"]]
            .sort_values(["Chrom","ExonStart","ExonEnd"])
        )
    except Exception:
        pass

    combined_df.to_excel(xw, index=False, sheet_name="Coverage gaps")

    # 4) Build from the thresholds-per-exon table if you have it
    try:
        # If you still have `th_exon` in scope (output of read_acmg_thresholds), this is best:
        genes_cov = build_genes_coverage(combined_df, sf_genes=sf_genes if len(sf_genes)>0 else None)
    except Exception:
        # Otherwise fall back to the combined_df you already created for Coverage gaps
        genes_cov = build_genes_coverage(combined_df, sf_genes=sf_genes if len(sf_genes)>0 else None)

    genes_cov.to_excel(xw, index=False, sheet_name="ACMG SF Genes Coverage")

    # 5) PASS variant table
    pass_cols = [
        "Variant","Gene","HGVSc","HGVSp","MANE_ID","Transcript","Consequence","GT","Zygosity","AD_Ref","AD_Alt","DP","GQ","QUAL",
        "FILTER","VAF","gnomAD_AF","ClinVar","ClinVar_ReviewStatus","ClinVar_Stars","ClinVar_StarsGlyph","REVEL","SpliceAI_DS_max","SpliceAI_Event","BayesDel_score","AM_Pathogenicity","AM_Class","HGVS_full",
        # Raw-alignment strand bias: advisory, never a filter. See strand_bias.py.
        "ALT_FWD","ALT_REV","StrandBias_P","StrandBias","QC_Flags"
    ]
    df[df["FILTER"]=="PASS"][pass_cols].to_excel(xw, index=False, sheet_name="PASS variants")

    # ---- Exomiser + Prioritised -------------------------------------------
    # Exomiser ranks the WHOLE callset against the patient's HPO terms,
    # including genes outside ACMG SF -- it runs on the unfiltered consensus,
    # not the panel-filtered VCF this report is built from.
    _exo = load_exomiser(args.exomiser)
    if _exo is not None:
        _sb_map  = load_strand_bias(args.strand_bias)
        _vaf_map = load_sites_vaf(args.sites_vcf)
        _mane    = load_mane(args.mane)

        _view = attach_strand_bias(exomiser_view(_exo), _sb_map, _vaf_map)
        if args.exomiser_top and args.exomiser_top > 0:
            _view = _view.head(args.exomiser_top)
        _view.to_excel(xw, index=False, sheet_name="Exomiser")

        # Prioritised: a SUPERSET, not an intersection.
        #   PANEL_REPORTABLE  every row of the ACMG SF sheet, with the Exomiser
        #                     rank/score attached where it ranked them and blank
        #                     where it did not. Exomiser filters before it ranks
        #                     -- intronic, synonymous and common variants are
        #                     dropped -- so building this from its output alone
        #                     would silently lose most reported variants.
        #   OFF_PANEL         Exomiser hits outside the ACMG SF genes, scoring
        #                     at least --prioritised-min-score.
        # What counts as "on panel". Prefer an explicit list: --panel-genes,
        # else --sf-genes if one was supplied. Failing both, fall back to the
        # genes present in this report -- the BED this pipeline filters on IS
        # the ACMG SF region set, so those genes are the panel, though a panel
        # gene with no calls in this sample will be missing from the fallback.
        # Exomiser sees the unfiltered consensus, so anything it names that is
        # absent from the panel fell outside those regions.
        _panel_label = "ACMG SF"
        _panel_genes = set()
        if args.panel_genes and os.path.exists(args.panel_genes):
            with open(args.panel_genes) as _f:
                _panel_genes = {l.split()[0] for l in _f if l.strip()}
        _panel_genes = _panel_genes or sf_genes \
                       or (set(df["Gene"].dropna().astype(str)) - {"", "NA"})

        _exo_by_variant = {str(r["Variant"]): r for _, r in exomiser_view(_exo).iterrows()}
        _sb_by_variant  = {str(r["Variant"]): r for _, r in _view.iterrows()}

        _rows, _seen = [], set()
        for _, r in _reportable.iterrows():
            _key = str(r.get("Variant"))
            if _key in _seen:
                continue
            _seen.add(_key)
            _e = _exo_by_variant.get(_key)
            _rows.append({
                "Scope": "PANEL_REPORTABLE", "Panel": _panel_label,
                "Gene": r.get("Gene"), "Variant": r.get("Variant"),
                "HGVSc": r.get("HGVSc"), "HGVSp": r.get("HGVSp"),
                "MANE_Select": r.get("MANE_ID"),
                "Consequence": r.get("Consequence"), "Zygosity": r.get("Zygosity"),
                "VAF": r.get("VAF"),
                "ClinVar": r.get("ClinVar"), "ClinVar_Stars": r.get("ClinVar_Stars"),
                "gnomAD_AF": r.get("gnomAD_AF"),
                "ALT_FWD": r.get("ALT_FWD"), "ALT_REV": r.get("ALT_REV"),
                "StrandBias_P": r.get("StrandBias_P"), "StrandBias": r.get("StrandBias"),
                "MOI": _e.get("MOI") if _e is not None else None,
                "Exo_Rank": _e.get("Rank") if _e is not None else None,
                "Exo_Score": _e.get("Combined_Score") if _e is not None else None,
                "Exo_Disease": _e.get("Exomiser_Disease") if _e is not None else None,
            })

        _ev = exomiser_view(_exo)
        _comb = pd.to_numeric(_ev["Combined_Score"], errors="coerce")
        for _i, _e in _ev.iterrows():
            _key = str(_e["Variant"])
            if _key in _seen:
                continue
            if not (pd.notna(_comb[_i]) and _comb[_i] >= args.prioritised_min_score):
                continue
            _tx, _c, _pr = split_exomiser_hgvs(_e.get("HGVS"))
            _sb = _sb_by_variant.get(_key)
            # Panel says where the GENE sits; Scope says why the row is here.
            # An Exomiser hit in a panel gene that did not qualify for the P-LP
            # sheet is PANEL_VUS, not OFF_PANEL -- the gene is on panel, only
            # the variant fell short of the reportable classification.
            _on_panel = str(_e.get("Gene")) in _panel_genes
            _rows.append({
                "Scope": "PANEL_VUS" if _on_panel else "OFF_PANEL",
                "Panel": _panel_label if _on_panel else "",
                "Gene": _e.get("Gene"), "Variant": _e.get("Variant"),
                "HGVSc": _c, "HGVSp": _pr,
                "MANE_Select": to_refseq(_tx, _e.get("Gene"), _mane),
                "Consequence": _e.get("Consequence"), "Zygosity": _e.get("Genotype"),
                "VAF": (_sb.get("VAF") if _sb is not None else None),
                "ClinVar": _e.get("ClinVar"), "ClinVar_Stars": _e.get("ClinVar_Stars"),
                # Exomiser reports MAX_FREQ as a PERCENTAGE, across all its
                # frequency sources rather than gnomAD alone.
                "gnomAD_AF": (float(_e["Max_Freq_pct"]) / 100.0
                              if pd.notna(_e.get("Max_Freq_pct"))
                              and str(_e.get("Max_Freq_pct")).strip() != "" else None),
                "ALT_FWD": (_sb.get("ALT_FWD") if _sb is not None else None),
                "ALT_REV": (_sb.get("ALT_REV") if _sb is not None else None),
                "StrandBias_P": (_sb.get("StrandBias_P") if _sb is not None else None),
                "StrandBias": (_sb.get("StrandBias") if _sb is not None else None),
                "MOI": _e.get("MOI"),
                "Exo_Rank": _e.get("Rank"), "Exo_Score": _e.get("Combined_Score"),
                "Exo_Disease": _e.get("Exomiser_Disease"),
            })

        _prio = pd.DataFrame(_rows, columns=[
            "Scope", "Panel", "Gene", "Variant", "HGVSc", "HGVSp", "MANE_Select",
            "Consequence", "Zygosity", "VAF", "ClinVar", "ClinVar_Stars", "gnomAD_AF",
            "ALT_FWD", "ALT_REV", "StrandBias_P", "StrandBias",
            "MOI", "Exo_Rank", "Exo_Score", "Exo_Disease"])
        # Reportable first, then on-panel VUS, then everything Exomiser pulled
        # in from outside the panel -- reviewer reading order, not alphabetical.
        _scope_order = {"PANEL_REPORTABLE": 0, "PANEL_VUS": 1, "OFF_PANEL": 2}
        _prio["_o"] = _prio["Scope"].map(_scope_order).fillna(9)
        _prio.sort_values(["_o", "Exo_Rank", "Gene"], kind="stable") \
             .drop(columns=["_o"]) \
             .to_excel(xw, index=False, sheet_name="Prioritised")

print(f"Wrote Excel report → {args.xlsx_out}")
