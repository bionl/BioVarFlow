#!/usr/bin/env python3
"""
Per-variant strand bias from the RAW alignment.

Background
----------
The pipeline already emits two strand QC files, but both are per-EXON and
allele-agnostic (ForwardReverseRatio, R1R2Ratio): they count every read over a
region regardless of which allele it carries, so a skew confined to the alt
reads is invisible to them. GATK's INFO/FS *is* per-allele, but it is computed
on post-reassembly allele depths — for MSH2 chr2:47414420 T>G it reported
FS=43.4, under the standard FS>60 filter, while the raw pileup at that site was
0 forward / 14 reverse alt reads against a reference that was 86% forward
(Fisher p ~ 5e-08, i.e. FS ~ 73). Reassembly had also inflated alt support
better than two-fold (AD 8,31 vs a raw pileup of T:57 G:11).

So this recomputes the 2x2 table from `bcftools mpileup`'s ADF/ADR, which are
counted straight off the alignment and — unlike a plain samtools pileup string —
cover indels as well as SNVs.

Metric and threshold
--------------------
Two-sided Fisher exact on

        ref_fwd  ref_rev
        alt_fwd  alt_rev

flagged at p < 1e-4. This is the reference-aware form: it asks whether the alt
reads are more skewed than the reference reads at the same site, so a locus
whose coverage is simply lopsided (capture edge, one-directional tiling) does
not flag on that account alone.

A fixed ratio cut such as 60/40 was rejected during calibration: it ignores
depth (4 vs 2 reads is not evidence) and ignores the site's own baseline, and
on real data flagged 282/524 variants (53.8%) in IQMM. p < 1e-4 flagged 2/524
in IQMM and 2/539 in VDDY (0.4% in both) while still catching MSH2.

Output is advisory. The report adds ALT_FWD / ALT_REV / StrandBias_P columns
and a QC_Flags entry; nothing is filtered out, because a flagged variant may
still be real and the decision (typically orthogonal confirmation) is the
reviewer's.

Scope: SNVs only. mpileup realigns gaps with its own model and does not
reproduce the caller's indel alleles -- on IQMM, 12 of the 14 largest
pileup-vs-caller VAF disagreements were indels (ACTC1 TCACA>T read 0.24 against
the caller's 1.00). Indels are reported as NA rather than given a p-value that
would look like evidence.
"""

import argparse
import collections
import gzip
import math
import sys

# Flag threshold. FS=60, the GATK convention, is phred for p=1e-6; 1e-4 is
# deliberately more sensitive because it runs on raw counts, which are noisier
# than the caller's polished ones, and only raises a flag rather than filtering.
DEFAULT_P = 1e-4


def _lchoose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def fisher_exact_two_sided(a, b, c, d):
    """Two-sided Fisher exact p for [[a,b],[c,d]]. No scipy dependency."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    row1, col1 = a + b, a + c
    lo = max(0, col1 - (c + d))
    hi = min(row1, col1)
    denom = _lchoose(n, col1)

    def lp(x):
        return _lchoose(row1, x) + _lchoose(n - row1, col1 - x) - denom

    obs = lp(a)
    # Sum every table at most as probable as the observed one. The 1e-7 slack
    # absorbs floating-point noise between mathematically equal probabilities.
    total = 0.0
    for x in range(lo, hi + 1):
        v = lp(x)
        if v <= obs + 1e-7:
            total += math.exp(v)
    return min(1.0, total)


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_vcf_variants(path):
    """(chrom, pos, ref, alt) for every record. VCF is already normalized."""
    out = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            for alt in f[4].split(","):
                out.append((f[0], int(f[1]), f[3].upper(), alt.upper()))
    return out


def parse_mpileup(path):
    """{(chrom, pos): [(ref, [alts], [adf], [adr]), ...]} from the query TSV.

    Columns are CHROM POS REF ALT ADF ADR, with ADF/ADR comma-separated in
    allele order (reference first, then each ALT, including mpileup's <*>).

    A position maps to a LIST, not a single record: bcftools mpileup emits the
    SNV and the indel at a site as separate lines sharing the same POS. Keying
    one record per position dropped whichever came first -- at MSH2
    chr2:47414420, which carries both a T>G SNV and a TAA>T deletion, that left
    both VCF rows untested.
    """
    table = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            chrom, pos, ref, alt, adf, adr = f[0], int(f[1]), f[2].upper(), f[3].upper(), f[4], f[5]

            def nums(s):
                return [int(x) if x.isdigit() else 0 for x in s.split(",")]

            table[(chrom, pos)].append((ref, alt.split(","), nums(adf), nums(adr)))
    return table


def match_allele(v_ref, v_alt, m_ref, m_alts):
    """Index into mpileup's allele list (0 = reference) for this VCF ALT.

    SNVs only -- indels are filtered out before this is called, because
    mpileup's gap realignment does not reproduce the caller's indel alleles.
    """
    if len(v_ref) != 1 or len(v_alt) != 1 or len(m_ref) != 1:
        return None
    for i, a in enumerate(m_alts):
        if a == v_alt:
            return i + 1
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mpileup", required=True,
                    help="TSV: CHROM POS REF ALT ADF ADR from bcftools mpileup|query")
    ap.add_argument("--vcf", required=True, help="Normalized (single-ALT) VCF")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--p-threshold", type=float, default=DEFAULT_P,
                    help=f"Flag variants below this Fisher p (default {DEFAULT_P:g})")
    args = ap.parse_args()

    mp = parse_mpileup(args.mpileup)
    variants = read_vcf_variants(args.vcf)

    flagged = unmatched = skipped_indel = 0
    with open(args.out, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tREF_FWD\tREF_REV\tALT_FWD\tALT_REV\t"
                  "ALT_FWD_FRAC\tStrandBias_P\tStrandBias_Flag\n")
        for chrom, pos, ref, alt in variants:
            na = f"{chrom}\t{pos}\t{ref}\t{alt}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"

            # Indels are not tested. bcftools mpileup realigns gaps with its own
            # model and does not reproduce the caller's indel alleles: on IQMM,
            # 12 of the 14 largest pileup-vs-caller VAF disagreements were
            # indels, ACTC1 TCACA>T reading 0.24 against the caller's 1.00. A
            # Fisher p on counts that far off the called allele is noise, and
            # the p<1e-4 threshold was calibrated on SNVs only. Emitting NA says
            # "not tested"; emitting a number would imply it had been.
            if not (len(ref) == 1 and len(alt) == 1 and alt != "*"):
                skipped_indel += 1
                out.write(na)
                continue

            # A position can carry several mpileup records (SNV line + indel
            # line); take the first whose alleles match this VCF row.
            hit = None
            for rec in mp.get((chrom, pos), ()):
                idx = match_allele(ref, alt, rec[0], rec[1])
                if idx is not None and idx < len(rec[2]) and idx < len(rec[3]):
                    hit = (rec, idx)
                    break
            if hit is None:
                # No pileup evidence to judge on: emit NA rather than a
                # fabricated zero, so the report can tell "not tested" from
                # "tested and clean".
                unmatched += 1
                out.write(na)
                continue
            rec, idx = hit
            _, _, adf, adr = rec

            rf, rr = adf[0], adr[0]
            af, ar = adf[idx], adr[idx]
            p = fisher_exact_two_sided(rf, rr, af, ar)
            tot = af + ar
            frac = f"{af / tot:.3f}" if tot else "NA"
            flag = "PASS"
            if tot == 0:
                flag = "NA"
            elif p < args.p_threshold:
                flag = "STRAND_BIAS"
                flagged += 1
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{rf}\t{rr}\t{af}\t{ar}\t"
                      f"{frac}\t{p:.3g}\t{flag}\n")

    n = len(variants)
    pct = (100.0 * flagged / n) if n else 0.0
    print(f"strand_bias: {n} variants, {flagged} flagged ({pct:.1f}%) at "
          f"p<{args.p_threshold:g}, {unmatched} SNVs without pileup support, "
          f"{skipped_indel} indels not tested",
          file=sys.stderr)


if __name__ == "__main__":
    main()
