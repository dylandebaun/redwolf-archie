#!/usr/bin/env python3
"""Combine per-donor IBDmix runs: max LOD per site across donors.

A site is called coyote-shared when its max LOD over the chosen donors is
>= --lod, i.e. when it lies in a LOD >= --lod segment from at least one donor.
Per individual and chromosome the qualifying segments of all donors are merged
into non-overlapping intervals; each keeps the max LOD and the donors that
support it.

    python3 scripts/ibdmix/combine_donors.py --in-dir results/ibdmix_coyotes_53g \
        --donors ALcoy,ILcoy --out results/ibdmix_coyotes_53g/combined_ALcoy_ILcoy

With --mask-dir, each donor's segments are first cut at that donor's mask
(<mask-dir>/<donor>.bed, integer chromosomes): ibdmix -r only sets LOD to 0
inside the mask, so a segment can still run across it. After cutting, a
donor's runs of homozygosity contribute nothing; there, only the other donors
can call a segment.

Writes <out>/ibdmix.chrNN.out (ID chrom start end slod donors), readable by
results/ibdmix_alcoy_53g/plot_map.py, and <out>/per_sample_LOD4.tsv.
Run from the repository root."""
import argparse
import glob
import os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--in-dir", required=True, help="folder with <donor>/ibdmix.chrNN.out")
ap.add_argument("--donors", required=True, help="comma-separated donors to combine")
ap.add_argument("--out", required=True)
ap.add_argument("--lod", type=float, default=4.0)
ap.add_argument("--mask-dir", default="", help="folder with <donor>.bed masks")
ap.add_argument("--fai", default="/mendel-nas1/ddebaun/RedWolf/canfam31_plusY.fa.gz.fai")
a = ap.parse_args()
donors = a.donors.split(",")
os.makedirs(a.out, exist_ok=True)

masks = {}
if a.mask_dir:
    for d in donors:
        m = pd.read_csv(f"{a.mask_dir}/{d}.bed", sep="\t", header=None, names=["chrom", "s", "e"])
        masks[d] = {c: g.sort_values("s")[["s", "e"]].to_numpy() for c, g in m.groupby("chrom")}


def cut(seg, d):
    """remove the parts of every segment that fall in donor d's mask"""
    if d not in masks:
        return seg
    out = []
    for r in seg.itertuples(index=False):
        pieces = [(r.start, r.end)]
        for ms, me in masks[d].get(r.chrom, []):
            if me <= r.start or ms >= r.end:
                continue
            pieces = [(x, y) for a0, b0 in pieces for x, y in ((a0, min(b0, ms)), (max(a0, me), b0)) if y > x]
        out += [(r.ID, r.chrom, x, y, r.slod, r.donor) for x, y in pieces]
    return pd.DataFrame(out, columns=seg.columns)


chroms = sorted({os.path.basename(f) for d in donors for f in glob.glob(f"{a.in_dir}/{d}/ibdmix.chr*.out")})
missing = [(d, c) for d in donors for c in chroms if not os.path.exists(f"{a.in_dir}/{d}/{c}")]
if missing or len(chroms) != 38:
    raise SystemExit(f"expected 38 chromosomes for every donor; found {len(chroms)}, missing {missing[:5]}")

totals = {}
for c in chroms:
    parts = []
    for d in donors:
        x = pd.read_csv(f"{a.in_dir}/{d}/{c}", sep="\t", usecols=[0, 1, 2, 3, 4], dtype={"ID": str}).assign(donor=d)
        parts.append(cut(x[x.slod >= a.lod], d))
    seg = pd.concat(parts).sort_values(["ID", "start"])
    rows = []
    for sid, g in seg.groupby("ID", sort=False):
        cur = None
        for r in g.itertuples(index=False):
            if cur and r.start <= cur[1]:                       # overlaps/abuts: extend
                cur[1] = max(cur[1], r.end); cur[2] = max(cur[2], r.slod); cur[3].add(r.donor)
            else:
                if cur:
                    rows.append((sid, r.chrom, cur[0], cur[1], cur[2], ",".join(sorted(cur[3]))))
                cur = [r.start, r.end, r.slod, {r.donor}]
        rows.append((sid, g.chrom.iloc[0], cur[0], cur[1], cur[2], ",".join(sorted(cur[3]))))
    out = pd.DataFrame(rows, columns=["ID", "chrom", "start", "end", "slod", "donors"])
    out.to_csv(f"{a.out}/{c}", sep="\t", index=False)
    for sid, bp in (out.end - out.start).groupby(out.ID).sum().items():
        totals[sid] = totals.get(sid, 0) + bp

fai = pd.read_csv(a.fai, sep="\t", header=None, usecols=[0, 1])
G = fai[fai[0].str.fullmatch(r"chr\d\d")][1].sum()
pct = pd.Series(totals).div(G).mul(100).round(1).sort_values(ascending=False)
pct.rename("pct_autosomes_LOD4").to_csv(f"{a.out}/per_sample_LOD4.tsv", sep="\t", index_label="sample")
print(f"combined {donors} over {len(chroms)} chromosomes -> {a.out}")
print(pct.to_string())
