#!/usr/bin/env python3
"""Tiered filtering and annotation of CNVkit calls (KHCC germline exomes).

Usage (project-agnostic; every resource is a parameter):
  python3 filter_cnv.py --calls 'results/*/*.call.cns' --sex sex.tsv \
      --pon A=pon_A.cnn --pon B=pon_B.cnn [--batch batch.tsv] [--pon-size A=27 --pon-size B=4] \
      --panel-bed panel.bed --gnomad-sv gnomad_sv_DEL_DUP_CNV.tsv [--gnomad-exome-cnv gnomad.cnv.vcf.gz] \
      --clingen ClinGen_gene_curation_list_GRCh38.tsv --segdup genomicSuperDups.txt.gz \
      [--common-bed common_cnv.bed] [--metrics metrics.tsv] --outdir out --out-xlsx report.xlsx
  Legacy KHCC layout: python3 filter_cnv.py --run-dir new_run --annot annotation

Per sample it writes <outdir>/<S>.cnv_filtered.tsv (legacy: <run_dir>/<S>/), plus cnv_filtered_all.tsv and an
Excel workbook. Samples missing from --batch use the batch "default" (a single --pon PATH is batch "default").

Levels (every call keeps the reasons it passed or failed; nothing is dropped silently):
  1 sample QC   : bivar <= --max-bivar; batch reference size (batch B flagged low confidence)
  2 technical   : >= --min-probes exons; germline signal (DEL <= --del-log2, DUP >= --dup-log2)
                  or mosaic route (>= --mosaic-probes exons, |log2| >= --mosaic-log2; not on haploid male X/Y);
                  median depth >= --min-depth; median reference spread of its exons < --max-spread
  3 common      : >= 50% inside a cohort-recurrent region (common_cnv.bed) or gnomAD-SV PASS
                  same-type overlap (50% reciprocal) with frequency > --max-af;
                  segmental duplication > 50% is flagged, not dropped
  4 relevance   : hemonc/ACMG SF panel genes; ClinGen HI=3 (DEL) / TS=3 (DUP)
Tiers:
  Tier 1  passes 1-3 and hits a panel gene (segdup-overlapping -> Tier 2)
  Tier 2  passes 1-3, no panel gene, or panel gene inside a segdup
  Tier 3  filtered (technical or common), reason given
"""
import argparse, bisect, csv, glob, gzip, os, re, statistics as st
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--run-dir", default=None, help="legacy layout: <run_dir>/<S>/<S>.call.cns, pon_final.<batch>.cnn, sex.tsv ...")
ap.add_argument("--annot", default=None, help="legacy: folder holding the default-named annotation files")
ap.add_argument("--calls", nargs="+", default=None, help="CNVkit .call.cns files (globs allowed)")
ap.add_argument("--sex", default=None, help="sample<TAB>male|female")
ap.add_argument("--batch", default=None, help="sample<TAB>batch (optional)")
ap.add_argument("--pon", action="append", default=[], help="BATCH=reference.cnn, repeatable; or a single path")
ap.add_argument("--pon-size", action="append", default=[], help="BATCH=N normals in that PoN (for the low-confidence flag)")
ap.add_argument("--panel-bed", default=None, help="exon BED of reportable genes; name column GENE or GENE|TRANSCRIPT")
ap.add_argument("--panel-name", default="panel")
ap.add_argument("--gnomad-sv", default=None, help="slim gnomAD SV table (chrom start end svtype freq grpmax filter name)")
ap.add_argument("--gnomad-exome-cnv", default=None, help="gnomAD exome CNV VCF (optional)")
ap.add_argument("--clingen", default=None)
ap.add_argument("--segdup", default=None, help="UCSC genomicSuperDups.txt.gz")
ap.add_argument("--common-bed", default=None, help="cohort-recurrent CNV regions from the PoN build (optional)")
ap.add_argument("--metrics", default=None, help="cnvkit metrics TSV covering the samples (optional)")
ap.add_argument("--exclude", default=None, help="samples excluded from the PoN (legacy PoN-size count)")
ap.add_argument("--outdir", default=None)
ap.add_argument("--out-xlsx", default="CNV_tiers.xlsx")
ap.add_argument("--max-bivar", type=float, default=0.25)
ap.add_argument("--min-batch-refs", type=int, default=10)
ap.add_argument("--min-probes", type=int, default=3)
ap.add_argument("--del-log2", type=float, default=-0.4)
ap.add_argument("--dup-log2", type=float, default=0.3)
ap.add_argument("--mosaic-probes", type=int, default=50)
ap.add_argument("--mosaic-log2", type=float, default=0.2)
ap.add_argument("--min-depth", type=float, default=20)
ap.add_argument("--max-spread", type=float, default=0.3)
ap.add_argument("--max-af", type=float, default=0.01)
ap.add_argument("--recip", type=float, default=0.5)
ap.add_argument("--haploid-x-reference", action="store_true",
                help="the PoN was built with -y (male/haploid-X reference); expected chrX log2 flips")
ap.add_argument("--small-exons", type=int, default=5, help="calls with <= this many exons need orthogonal confirmation")
a = ap.parse_args()
R, A = a.run_dir, a.annot
def pick(val, legacy):
    if val:
        return val
    return legacy if legacy and os.path.exists(legacy) else None
L = (lambda x: f"{R}/{x}") if R else (lambda x: None)
LA = (lambda x: f"{A}/{x}") if A else (lambda x: None)
a.sex = pick(a.sex, L("sex.tsv")); a.batch = pick(a.batch, L("batch.tsv"))
a.panel_bed = pick(a.panel_bed, LA("ACMG_SF_plus_HemOnc_MANE_exons_50bp.refseq.bed"))
a.gnomad_sv = pick(a.gnomad_sv, LA("gnomad_v4.1_sv_DEL_DUP_CNV.tsv"))
a.gnomad_exome_cnv = pick(a.gnomad_exome_cnv, LA("gnomad.v4.1.cnv.all.vcf.gz"))
a.clingen = pick(a.clingen, LA("ClinGen_gene_curation_list_GRCh38.tsv"))
a.segdup = pick(a.segdup, LA("genomicSuperDups.txt.gz"))
a.common_bed = pick(a.common_bed, L("common_cnv.bed"))
a.metrics = pick(a.metrics, L("loo_metrics.tsv"))
a.exclude = pick(a.exclude, L("exclude.txt"))
for k in ("sex", "panel_bed", "gnomad_sv", "clingen", "segdup"):
    if not getattr(a, k):
        ap.error(f"--{k.replace('_', '-')} is required")
sex = dict(l.split()[:2] for l in open(a.sex) if l.strip())
batch = dict(l.split()[:2] for l in open(a.batch) if l.strip()) if a.batch else {}
pons = {}
for x in a.pon:
    b, _, p = x.rpartition("=")
    pons[b or "default"] = p
if not pons and R:
    for p in glob.glob(f"{R}/pon_final*.cnn"):
        m = re.match(r"pon_final\.?(.*)\.cnn$", os.path.basename(p))
        pons[m.group(1) or "default"] = p
if not pons:
    ap.error("--pon is required")
if a.calls:
    call_files = sorted({p for g in a.calls for p in glob.glob(g)})
else:
    call_files = sorted(glob.glob(f"{R}/*/*.call.cns"))
outdir = a.outdir or R or "."
os.makedirs(outdir, exist_ok=True)
def batch_of(s):
    b = batch.get(s, "default")
    return b if b in pons else ("default" if "default" in pons else next(iter(pons)))


class Index:
    """Per-chromosome interval index: intervals sorted by start, with max span for safe lookback."""
    def __init__(self):
        self.d = defaultdict(list)
    def add(self, c, s, e, v=None):
        self.d[c].append((s, e, v))
    def build(self):
        self.starts, self.maxlen = {}, {}
        for c, iv in self.d.items():
            iv.sort(key=lambda x: x[0])
            self.starts[c] = [x[0] for x in iv]
            self.maxlen[c] = max(x[1] - x[0] for x in iv)
        return self
    def query(self, c, s, e):
        if c not in self.d:
            return []
        iv, st_ = self.d[c], self.starts[c]
        i = bisect.bisect_left(st_, s - self.maxlen[c])
        j = bisect.bisect_left(st_, e)
        return [x for x in iv[i:j] if x[1] > s]


def covered_fraction(hits, s, e):
    """Fraction of [s,e) covered by the union of hit intervals."""
    segs = sorted((max(h[0], s), min(h[1], e)) for h in hits)
    tot, cur_s, cur_e = 0, None, None
    for x, y in segs:
        if cur_e is None or x > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = x, y
        else:
            cur_e = max(cur_e, y)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot / max(1, e - s)


def recip(s1, e1, s2, e2, f):
    ov = min(e1, e2) - max(s1, s2)
    return ov > 0 and ov >= f * (e1 - s1) and ov >= f * (e2 - s2)


print("Loading annotations ...")
panel = Index()
for l in open(a.panel_bed):
    f = l.rstrip("\n").split("\t")
    panel.add(f[0], int(f[1]), int(f[2]), f[3].split("|")[0])
panel.build()

segdup = Index()
with gzip.open(a.segdup, "rt") as fh:
    for l in fh:
        f = l.split("\t")
        segdup.add(f[1], int(f[2]), int(f[3]))
segdup.build()

clingen = Index()
for l in open(a.clingen):
    if l.startswith("#"):
        continue
    f = l.rstrip("\n").split("\t")
    m = re.match(r"(chr[\w]+):(\d+)-(\d+)", f[3].replace(" ", ""))
    if m:
        clingen.add(m.group(1), int(m.group(2)), int(m.group(3)), (f[0], f[4], f[12]))
clingen.build()

gsv = Index()
with open(a.gnomad_sv) as fh:
    next(fh)
    for l in fh:
        c, s, e, t, fr, gm, filt, name = l.rstrip("\n").split("\t")
        if filt != "PASS":
            continue
        try:
            fr = float(fr)
        except ValueError:
            continue
        gsv.add(c, int(s), int(e), (t, fr))
gsv.build()

gex = Index()
with (gzip.open(a.gnomad_exome_cnv, "rt") if a.gnomad_exome_cnv else open(os.devnull)) as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.split("\t", 8)
        if f[6] != "PASS":
            continue
        inf = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        gex.add(f[0], int(f[1]), int(inf["END"]), (inf["SVTYPE"], int(float(inf["SC"])), float(inf["SC"]) / float(inf["SN"])))
gex.build()

common = Index()
if a.common_bed:
    for l in open(a.common_bed):
        f = l.rstrip("\n").split("\t")
        common.add(f[0], int(f[1]), int(f[2]), f[4])
common.build()

spread = {}
for b, p in pons.items():
    idx = Index()
    rows = csv.reader(open(p), delimiter="\t")
    h = next(rows); ic, isp, ig = h.index("chromosome"), h.index("spread"), h.index("gene")
    for r in rows:
        if r[ig] != "Antitarget":
            idx.add(r[ic], int(r[1]), int(r[2]), float(r[isp]))
    spread[b] = idx.build()

metrics = {}
if a.metrics:
    rows = list(csv.reader(open(a.metrics), delimiter="\t"))
    for r in rows[1:]:
        metrics[os.path.basename(r[0]).replace(".cnr", "")] = float(r[rows[0].index("bivar")])
batch_n = defaultdict(int)
excl = set(l.strip() for l in open(a.exclude)) if a.exclude else set()
for x in a.pon_size:
    b, _, n_ = x.rpartition("=")
    batch_n[b or "default"] = int(n_)
if not a.pon_size:
    for s, b in batch.items():
        if s not in excl:
            batch_n[b] += 1


def expected(chrom, sx):
    """(expected copy number, expected log2) for this sample against the given reference.

    A female-normalised PoN (the default) puts a normal male chrX at log2 ~ -1; a PoN built
    with -y is already haploid on X, so a normal male sits at 0 and a female at +1.
    """
    if chrom == "chrX":
        if a.haploid_x_reference:
            return (1, 0.0) if sx == "male" else (2, 1.0)
        return (1, -1.0) if sx == "male" else (2, 0.0)
    if chrom == "chrY":
        return (1, 0.0) if sx == "male" else (0, None)
    return (2, 0.0)


cols = ["sample", "sex", "batch", "tier", "chrom", "start", "end", "size_kb", "type", "cn", "log2", "delta_log2",
        "est_cell_fraction", "exons", "depth", "ref_spread", "panel_genes", "clingen_dosage", "gnomad_sv_max_af",
        "gnomad_exome_cnv_carriers", "cohort_recurrent", "segdup_frac", "flags", "reasons", "all_genes"]
out_all = []
print("Filtering ...")
for f in call_files:
    s = os.path.basename(f).replace(".call.cns", "")
    if s not in sex:
        continue
    sx, b = sex[s], batch_of(s)
    bivar = metrics.get(s, float("nan"))
    rows = list(csv.reader(open(f), delimiter="\t"))
    ix = {k: i for i, k in enumerate(rows[0])}
    res = []
    for r in rows[1:]:
        c, st_, en = r[0], int(r[1]), int(r[2])
        cn = int(r[ix["cn"]]); lg = float(r[ix["log2"]]); n = int(r[ix["probes"]]); dp = float(r[ix["depth"]])
        ecn, elg = expected(c, sx)
        if elg is None or cn == ecn:
            continue
        t = "DEL" if cn < ecn else "DUP"
        d = lg - elg
        base = 2 if c not in ("chrX", "chrY") or (c == "chrX" and sx == "female") else 1
        ratio = 2 ** d
        cf = min(1.0, abs(ratio - 1) * base)  # one-copy change in a fraction of cells
        reasons, flags = [], []
        # level 1
        if bivar > a.max_bivar:
            reasons.append(f"sample noisy (bivar {bivar:.2f})")
        if batch_n.get(b, 99) < a.min_batch_refs:
            flags.append(f"batch {b}: only {batch_n[b]} normals")
        # level 2
        germ = (t == "DEL" and d <= a.del_log2) or (t == "DUP" and d >= a.dup_log2)
        mos = base == 2 and n >= a.mosaic_probes and abs(d) >= a.mosaic_log2
        if n < a.min_probes:
            reasons.append(f"<{a.min_probes} exons")
        elif not (germ or mos):
            reasons.append(f"weak signal (delta log2 {d:+.2f})")
        if mos and not germ:
            flags.append(f"mosaic-level (~{cf:.0%} of cells)")
        if dp < a.min_depth:
            reasons.append(f"low depth ({dp:.0f}x)")
        sp = [x[2] for x in spread[b].query(c, st_, en)]
        msp = st.median(sp) if sp else float("nan")
        if sp and msp >= a.max_spread:
            reasons.append(f"noisy region (ref spread {msp:.2f})")
        # level 3
        com = common.query(c, st_, en)
        cfrac = covered_fraction(com, st_, en) if com else 0
        if cfrac >= 0.5:
            reasons.append("recurrent in cohort")
        g_af = max([v[1] for (gs, ge, v) in gsv.query(c, st_, en)
                    if v[0] in (t, "CNV") and recip(st_, en, gs, ge, a.recip)], default=0.0)
        if g_af > a.max_af:
            reasons.append(f"common in gnomAD-SV (AF {g_af:.3f})")
        g_ex = max([v[1] for (gs, ge, v) in gex.query(c, st_, en)
                    if v[0] == t and recip(st_, en, gs, ge, a.recip)], default=0)
        sd = covered_fraction(segdup.query(c, st_, en), st_, en)
        if sd > 0.5:
            flags.append(f"segdup {sd:.0%}")
        # level 4
        pg = sorted({v for (_, _, v) in panel.query(c, st_, en)})
        cg = []
        for (_, _, (g, hi, ts)) in clingen.query(c, st_, en):
            if t == "DEL" and hi == "3":
                cg.append(f"{g}:HI3")
            if t == "DUP" and ts == "3":
                cg.append(f"{g}:TS3")
        if n <= a.small_exons and not reasons:
            flags.append("small: confirm (MLPA/ddPCR)")
        if reasons:
            tier = "Tier 3"
        elif pg and sd <= 0.5:
            tier = "Tier 1"
        else:
            tier = "Tier 2"
        genes = sorted({x for x in r[ix["gene"]].split(",") if x not in ("-", "Antitarget")})
        res.append(dict(sample=s, sex=sx, batch=b, tier=tier, chrom=c, start=st_, end=en,
                        size_kb=round((en - st_) / 1000, 1), type=t, cn=cn, log2=round(lg, 2), delta_log2=round(d, 2),
                        est_cell_fraction=round(cf, 2), exons=n, depth=round(dp), ref_spread=round(msp, 2) if sp else "",
                        panel_genes=",".join(pg), clingen_dosage=",".join(sorted(set(cg))), gnomad_sv_max_af=round(g_af, 4),
                        gnomad_exome_cnv_carriers=g_ex, cohort_recurrent="yes" if cfrac >= 0.5 else "",
                        segdup_frac=round(sd, 2), flags="; ".join(flags), reasons="; ".join(reasons),
                        all_genes=",".join(genes[:25]) + (" ..." if len(genes) > 25 else "")))
    per = f"{R}/{s}/{s}.cnv_filtered.tsv" if (R and not a.calls and not a.outdir) else f"{outdir}/{s}.cnv_filtered.tsv"
    with open(per, "w") as out:
        out.write("\t".join(cols) + "\n")
        for x in res:
            out.write("\t".join(str(x[k]) for k in cols) + "\n")
    out_all += res

with open(f"{outdir}/cnv_filtered_all.tsv", "w") as out:
    out.write("\t".join(cols) + "\n")
    for x in out_all:
        out.write("\t".join(str(x[k]) for k in cols) + "\n")

# Excel
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
wb = openpyxl.Workbook()
HF, HFont = PatternFill("solid", fgColor="2F6F73"), Font(bold=True, color="FFFFFF")
LOSS, GAIN = PatternFill("solid", fgColor="E5ECFB"), PatternFill("solid", fgColor="FBE8E6")
widths = dict(sample=9, sex=8, batch=6, tier=8, chrom=7, start=12, end=12, size_kb=9, type=6, cn=4, log2=7,
              delta_log2=8, est_cell_fraction=9, exons=7, depth=7, ref_spread=8, panel_genes=28, clingen_dosage=20,
              gnomad_sv_max_af=10, gnomad_exome_cnv_carriers=10, cohort_recurrent=9, segdup_frac=8, flags=34,
              reasons=34, all_genes=50)


def table(ws, rows):
    ws.append(cols)
    for cell in ws[1]:
        cell.font, cell.fill = HFont, HF
        cell.alignment = Alignment(wrap_text=True, vertical="center")
    for x in rows:
        ws.append([x[k] for k in cols])
    for i, k in enumerate(cols, 1):
        ws.column_dimensions[get_column_letter(i)].width = widths[k]
    ti = cols.index("type")
    for row in ws.iter_rows(min_row=2):
        row[ti].fill = LOSS if row[ti].value == "DEL" else GAIN
        row[cols.index("start")].number_format = row[cols.index("end")].number_format = "#,##0"
    ws.freeze_panes = "E2"
    ws.auto_filter.ref = ws.dimensions


ws = wb.active
ws.title = "Summary"
ws.append(["Sample", "Sex", "Batch", "bivar", "Tier 1", "Tier 2", "Tier 3", "Tier 1 genes"])
for cell in ws[1]:
    cell.font, cell.fill = HFont, HF
by = defaultdict(list)
for x in out_all:
    by[x["sample"]].append(x)
for s in sorted(by):
    v = by[s]
    t1 = [x for x in v if x["tier"] == "Tier 1"]
    ws.append([s, sex[s], batch.get(s, ""), metrics.get(s), len(t1), sum(x["tier"] == "Tier 2" for x in v),
               sum(x["tier"] == "Tier 3" for x in v), "; ".join(f'{x["type"]} {x["panel_genes"]}' for x in t1)])
for i, w in enumerate([10, 9, 7, 8, 8, 8, 8, 80], 1):
    ws.column_dimensions[get_column_letter(i)].width = w
table(wb.create_sheet("Tier 1"), [x for x in out_all if x["tier"] == "Tier 1"])
table(wb.create_sheet("Tier 2"), [x for x in out_all if x["tier"] == "Tier 2"])
table(wb.create_sheet("Tier 3 (filtered)"), [x for x in out_all if x["tier"] == "Tier 3"])
ws = wb.create_sheet("Thresholds")
for k, v in vars(a).items():
    ws.append([k, ", ".join(map(str, v)) if isinstance(v, list) else v])
ws.append([])
ws.append(["Panel", f"{a.panel_name}: {a.panel_bed}"])
ws.append(["Annotations", f"gnomAD SV: {a.gnomad_sv}; gnomAD exome CNV: {a.gnomad_exome_cnv}; ClinGen: {a.clingen}; segdups: {a.segdup}; cohort common CNVs: {a.common_bed}"])
for b, p in pons.items():
    ws.append([f"PoN {b}", f"{p} ({batch_n.get(b, '?')} normals)"])
ws.column_dimensions["A"].width, ws.column_dimensions["B"].width = 22, 100
wb.save(a.out_xlsx)

n = defaultdict(int)
for x in out_all:
    n[x["tier"]] += 1
print(f"{len(by)} samples, {len(out_all)} calls: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())))
print(f"Wrote per-sample *.cnv_filtered.tsv, {outdir}/cnv_filtered_all.tsv, {a.out_xlsx}")
