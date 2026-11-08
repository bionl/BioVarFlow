#!/usr/bin/env python3
"""Check a CNVkit panel of normals against the manifest.json the cnv_normal pipeline ships with it.

A PoN delivery is pon.<batch>.cnn + common_cnv.bed + manifest.json. The manifest records what the
PoN was built from, which lets the run fail before calling instead of producing quietly wrong
copy numbers. Two findings stop the run; the rest warn.

STOP:
  * the target BED's md5 differs from the PoN's -- the bins would differ, so every log2 is wrong.
    Stronger than cnv_pon_check.py's bin-count comparison, which only sees a size mismatch.
  * the reference sex model disagrees with --haploid-x-reference. Getting this wrong makes a
    normal male chrX read as cn = 0, which is a reportable-looking finding that is pure artifact.

WARN (the run continues; a warning is better than refusing a usable third-party PoN):
  * cnvkit version differs from the container's
  * no `pons` entry matches the .cnn given, so n_normals cannot be checked
  * n_normals below --min-pon-size (the low-confidence flag already covers the calls themselves)
  * common_cnv.bed md5 differs -- probably a BED from a different PoN, which is worse than none

Writes a two-column provenance TSV for the report, and prints every check to stdout for the log.
With no manifest, writes a provenance file saying so and exits 0.
"""
import argparse
import hashlib
import json
import os
import sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--manifest", help="manifest.json; omit or pass a missing path to skip the checks")
ap.add_argument("--pon", required=True, help="the .cnn actually being used")
ap.add_argument("--targets-bed", required=True, help="cnv_targets_bed, whose md5 must match the manifest")
ap.add_argument("--common-bed", help="cnv_common_bed, checked only when given")
ap.add_argument("--cnvkit-version", default="0.9.13", help="the version in the CNVKIT_POSTPROCESS container")
ap.add_argument("--min-pon-size", type=int, default=10)
ap.add_argument("--haploid-x-reference", action="store_true",
                help="set when cnv_pon_haploid_x is true, i.e. the PoN was built with -y")
ap.add_argument("--out", required=True, help="provenance TSV for the report")
a = ap.parse_args()

prov, warns, stops = [], [], []


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_and_exit():
    with open(a.out, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in prov:
            fh.write(f"{k}\t{v}\n")
    for w in warns:
        print(f"WARNING: {w}", file=sys.stderr)
    if stops:
        for s in stops:
            print(f"ERROR: {s}", file=sys.stderr)
        sys.exit(1)
    sys.exit(0)


pon_name = os.path.basename(a.pon)
prov.append(("pon_file", pon_name))

if not a.manifest or not os.path.exists(a.manifest):
    # An older or third-party PoN. Say so in the report rather than implying it was verified.
    msg = "no manifest.json found next to the PoN; its provenance is unverified"
    warns.append(msg + f". Looked for: {a.manifest or '<not given>'}")
    prov.append(("manifest", "absent -- PoN provenance unverified"))
    write_and_exit()

with open(a.manifest) as fh:
    m = json.load(fh)

prov.append(("manifest", os.path.basename(a.manifest)))
for k in ("pipeline", "pipeline_version", "created", "cnvkit_version", "reference_sex_model"):
    if k in m:
        prov.append((k, m[k]))

# ---- target BED: a mismatch means different bins, so stop ------------------------------------
want = (m.get("targets_bed") or {}).get("md5")
got = md5(a.targets_bed)
prov.append(("targets_bed", f"{os.path.basename(a.targets_bed)} (md5 {got})"))
if want and want != got:
    stops.append(
        f"the PoN was built from a different target BED. manifest {(m['targets_bed']).get('file')} "
        f"md5 {want}, but cnv_targets_bed {a.targets_bed} is md5 {got}. The bins would not match, "
        f"so every log2 ratio would be wrong. Supply the BED this PoN was built from.")
elif not want:
    warns.append("the manifest records no targets_bed md5, so the BED could not be verified")
else:
    print(f"OK   target BED md5 matches the manifest ({got})")

# ---- sex model: a mismatch silently corrupts chrX, so stop ------------------------------------
model = m.get("reference_sex_model", "")
is_female_norm = model.startswith("female-normalised")
if model:
    if is_female_norm and a.haploid_x_reference:
        stops.append(
            f"the PoN is female-normalised ({model!r}) but cnv_pon_haploid_x is true. A -y flag "
            f"against a female-normalised reference makes a normal male chrX read as cn = 0. "
            f"Set cnv_pon_haploid_x = false.")
    elif not is_female_norm and not a.haploid_x_reference:
        stops.append(
            f"the PoN's sex model is {model!r}, which is not female-normalised, but "
            f"cnv_pon_haploid_x is false. chrX expectations would be off by one copy. "
            f"Set cnv_pon_haploid_x = true, or use a female-normalised PoN.")
    else:
        print(f"OK   sex model agrees with cnv_pon_haploid_x={a.haploid_x_reference}")
else:
    warns.append("the manifest records no reference_sex_model, so the chrX model could not be verified")

# ---- cnvkit version: warn only ---------------------------------------------------------------
mv = m.get("cnvkit_version")
if mv and mv != a.cnvkit_version:
    warns.append(f"the PoN was built with cnvkit {mv} but this run uses {a.cnvkit_version}. "
                 f"Bin definitions and `fix` behaviour can differ between versions.")
elif mv:
    print(f"OK   cnvkit version matches ({mv})")

# ---- this PoN's entry: normals count and provenance ------------------------------------------
entry, batch = None, None
for b, p in (m.get("pons") or {}).items():
    if p.get("file") == pon_name:
        entry, batch = p, b
        break

if entry is None:
    warns.append(f"no entry for {pon_name} in the manifest's `pons` "
                 f"({', '.join(sorted((p.get('file') or '?') for p in (m.get('pons') or {}).values())) or 'none'}), "
                 f"so its normals count could not be checked")
    prov.append(("pon_normals", "unknown -- no manifest entry for this .cnn"))
else:
    n = entry.get("n_normals")
    prov.append(("pon_batch", batch))
    prov.append(("pon_normals", f"{n} ({entry.get('n_female', '?')} female, {entry.get('n_male', '?')} male)"))
    prov.append(("pon_target_bins", entry.get("target_bins", "?")))
    prov.append(("pon_normals_list", ", ".join(entry.get("normals") or []) or "not recorded"))
    if entry.get("md5") and entry["md5"] != md5(a.pon):
        warns.append(f"{pon_name} does not match its manifest md5; the file may have been "
                     f"modified or swapped since the manifest was written")
    if n is not None and n < a.min_pon_size:
        warns.append(f"PoN {batch} has only {n} normals, below cnv_min_pon_size={a.min_pon_size}. "
                     f"Its calls carry the low-confidence flag.")
    else:
        print(f"OK   PoN {batch}: {n} normals")

if m.get("excluded"):
    prov.append(("excluded_from_pon", ", ".join(m["excluded"])))

# ---- common_cnv.bed: warn only ---------------------------------------------------------------
if a.common_bed:
    cwant = (m.get("common_cnv_bed") or {}).get("md5")
    cgot = md5(a.common_bed)
    prov.append(("common_cnv_bed", f"{os.path.basename(a.common_bed)} (md5 {cgot})"))
    if cwant and cwant != cgot:
        warns.append(f"cnv_common_bed md5 {cgot} does not match the manifest's {cwant}. This BED "
                     f"is probably from a different PoN, which is worse than supplying none: its "
                     f"regions do not describe this cohort's recurrent artifacts.")
    elif cwant:
        print(f"OK   common_cnv.bed md5 matches the manifest ({cgot})")
else:
    prov.append(("common_cnv_bed", "not supplied"))

write_and_exit()
