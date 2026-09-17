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

Indels
------
Indels are matched to the pileup's own alleles on net length change, because
bcftools and the caller frequently pad the same event differently. An indel is
then tested only when the pileup reproduces the caller's VAF within
--indel-max-vaf-delta: where the two disagree, a realignment mismatch cannot be
told apart from a real artifact, so it is reported as not tested rather than
given a misleading p-value. On IQMM this makes 98 of 134 indels testable, up
from none.

`bcftools call -C alleles` was tried and rejected. It matches more indels, but
`call` makes a genotype decision, and on weak one-sided alt support it calls
hom-ref and drops the ALT -- returning ALT='.' for MSH2, RUNX1 and DSG2, and
flagging nothing at all.

That gate is deliberately not applied to SNVs: there, a large VAF disagreement
IS the artifact signal -- it is what caught MSH2.
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
    """(chrom, pos, ref, alt, caller_vaf) per record. VCF is already normalized.

    caller_vaf comes from FORMAT/AD and is None when AD is absent. It is used
    only as a reconciliation check on indels -- see the gate in main().
    """
    out = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            vaf = None
            if len(f) >= 10:
                keys = f[8].split(":")
                if "AD" in keys:
                    vals = f[9].split(":")
                    try:
                        ad = [int(x) for x in vals[keys.index("AD")].split(",")]
                        tot = sum(ad)
                        if tot and len(ad) >= 2:
                            vaf = sum(ad[1:]) / tot
                    except ValueError:
                        pass
            for alt in f[4].split(","):
                out.append((f[0], int(f[1]), f[3].upper(), alt.upper(), vaf))
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
    """Index into the pileup's allele list (0 = reference) for this VCF ALT.

    Exact REF/ALT identity first. Failing that, indels are matched on net
    length change: bcftools and the caller frequently pad the same event
    differently -- chr1:55056226 GGGGCGGA>G comes back as
    GGGGCGGAGGGCGGAGGGCGGAGGG>GGGGCGGAGGGCGGAGGG, the same 7 bp deletion
    written against a longer repeat context. On IQMM this recovered 42 indels
    and every one of them matched a single candidate, no ties.

    A tie is left unmatched rather than guessed at.
    """
    if m_ref == v_ref:
        for i, a in enumerate(m_alts):
            if a == v_alt:
                return i + 1

    if len(v_ref) == 1 and len(v_alt) == 1:
        return None                      # SNVs: exact identity or nothing

    want = len(v_alt) - len(v_ref)
    hits = [i for i, a in enumerate(m_alts)
            if a not in (".", "<*>") and len(a) - len(m_ref) == want]
    return hits[0] + 1 if len(hits) == 1 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mpileup", required=True,
                    help="TSV: CHROM POS REF ALT ADF ADR from bcftools mpileup|query")
    ap.add_argument("--vcf", required=True, help="Normalized (single-ALT) VCF")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--p-threshold", type=float, default=DEFAULT_P,
                    help=f"Flag variants below this Fisher p (default {DEFAULT_P:g})")
    ap.add_argument("--indel-max-vaf-delta", type=float, default=0.15,
                    help="Indels only: skip the test unless the pileup VAF is within this "
                         "of the caller's, i.e. unless the two agree on how the indel is "
                         "represented (default 0.15). Not applied to SNVs.")
    args = ap.parse_args()

    mp = parse_mpileup(args.mpileup)
    variants = read_vcf_variants(args.vcf)

    flagged = unmatched = skipped_indel = underpowered = hom_alt = 0
    indel_unreconciled = 0
    with open(args.out, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tREF_FWD\tREF_REV\tALT_FWD\tALT_REV\t"
                  "ALT_FWD_FRAC\tStrandBias_P\tStrandBias_MinP\tStrandBias_Flag\n")
        for chrom, pos, ref, alt, caller_vaf in variants:
            def na(reason):
                return (f"{chrom}\t{pos}\t{ref}\t{alt}\t"
                        f"NA\tNA\tNA\tNA\tNA\tNA\tNA\t{reason}\n")

            is_indel = not (len(ref) == 1 and len(alt) == 1)
            if alt == "*":
                # Spanning deletion placeholder: no allele of its own to count.
                skipped_indel += 1
                out.write(na("NOT_TESTED_INDEL"))
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
                # `call -C alleles` writes ALT='.' when it cannot reconcile the
                # requested allele with its own view of the pileup. This must
                # NOT be read as "zero alt reads": chr11:119027725 AC>A is a
                # hom-ALT deletion with AD 0,235 that comes back ALT='.' against
                # 187 reference reads. Scoring those as alt=0 would flag real
                # homozygous variants as artifacts.
                recs = mp.get((chrom, pos), ())
                reason = ("NOT_TESTED_NO_ALT_ALLELE"
                          if recs and all(r[1] in (["."], ["<*>"]) for r in recs)
                          else "NOT_TESTED_NO_PILEUP")
                unmatched += 1
                out.write(na(reason))
                continue
            rec, idx = hit
            _, _, adf, adr = rec

            rf, rr = adf[0], adr[0]
            af, ar = adf[idx], adr[idx]
            tot = af + ar
            if tot == 0:
                unmatched += 1
                out.write(na("NOT_TESTED_NO_ALT_READS"))
                continue

            # Hom-ALT: no reference reads, which is simply what homozygous
            # means -- 39% of tested SNVs in IQMM, the expected Hardy-Weinberg
            # proportion for germline data.
            #
            # This test looks for alt reads that should not be there: a false
            # allele produced by reassembly or misalignment, betrayed by sitting
            # on one strand while the reference reads do not. A hom-ALT call at
            # depth has no reference row to compare against, and nothing for the
            # test to add -- such a call is already unambiguous. So this is the
            # test correctly standing aside, not a gap in coverage, and it is
            # distinct from UNDERPOWERED, which more reads would fix.
            #
            # The case worth watching -- a 1/1 call that unexpectedly retains
            # reference reads -- is covered separately by HOM_REF_SUPPORT.
            if rf + rr == 0:
                hom_alt += 1
                out.write(na("NOT_TESTED_HOM_ALT"))
                continue

            # Indels only: require the pileup to reproduce the caller's VAF
            # before trusting its ref/alt split.
            #
            # For a SNV, "which allele does this read carry" is one base and
            # never ambiguous. For an indel it depends on where the aligner puts
            # the gap, and bcftools and the caller use different realignment
            # models -- ACTC1 chr15:34791307 TCACA>T reads 0.24 here against the
            # caller's 1.00, with 47 reads changing sides. Constraining the
            # pileup with `bcftools call -C alleles` fixes most of this, but not
            # all, so where the two still disagree we cannot tell a realignment
            # mismatch from a real artifact and decline to guess.
            #
            # This gate is NOT applied to SNVs: there a large VAF disagreement
            # is the artifact signal itself. It is what caught MSH2, whose
            # caller VAF of 0.79 collapsed to 0.18 in the raw pileup.
            if is_indel and caller_vaf is not None:
                if abs((af + ar) / (rf + rr + af + ar) - caller_vaf) > args.indel_max_vaf_delta:
                    indel_unreconciled += 1
                    out.write(na("NOT_TESTED_INDEL_UNRECONCILED"))
                    continue

            p = fisher_exact_two_sided(rf, rr, af, ar)

            # Power: the smallest p this site could have produced, i.e. what a
            # maximally one-sided alt would score against this reference. If
            # even that cannot clear the threshold, the variant was never
            # testable and "PASS" would be an absence of evidence dressed up as
            # evidence of absence -- SEC23B chr20:18524919 has 9 alt reads and
            # tops out at p=0.01, so it would have passed whatever the reads did.
            min_p = min(fisher_exact_two_sided(rf, rr, 0, tot),
                        fisher_exact_two_sided(rf, rr, tot, 0))

            frac = f"{af / tot:.3f}"
            if p < args.p_threshold:
                flag = "STRAND_BIAS"
                flagged += 1
            elif min_p >= args.p_threshold:
                flag = "UNDERPOWERED"
                underpowered += 1
            else:
                flag = "PASS"
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{rf}\t{rr}\t{af}\t{ar}\t"
                      f"{frac}\t{p:.3g}\t{min_p:.3g}\t{flag}\n")

    n = len(variants)
    pct = (100.0 * flagged / n) if n else 0.0
    print(f"strand_bias: {n} variants, {flagged} flagged ({pct:.1f}%) at "
          f"p<{args.p_threshold:g}, {unmatched} SNVs without pileup support, "
          f"{underpowered} underpowered, {hom_alt} hom-alt (no reference reads), "
          f"{indel_unreconciled} indels unreconciled, {skipped_indel} skipped",
          file=sys.stderr)


if __name__ == "__main__":
    main()
