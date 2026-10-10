#!/usr/bin/env python3
"""Expected red wolf ancestry per individual under the pooled ("all fractions")
model, by panel. Run from the repository root after summarize.py.

  filled dot   expected red wolf ancestry = mean P(red wolf) over all haplotype-
               windows (pooled model; well calibrated on simulated test data)
  bar          range of that value across the fixed-fraction models 0.30-0.60
  open circle  windows called red wolf at 5% FPR (pooled model; conservative,
               recall at that threshold is 0.67-0.77)

Writes summary/pooled_ancestry.{png,pdf} (table: summary/pooled_ancestry.tsv).
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

S = "sensitivity/admixture/summary"
t = pd.read_csv(f"{S}/pooled_ancestry.tsv", sep="\t", dtype={"sample": str})
PANELS = [("historical", "Historical 1975 + SAFE"), ("museum_early", "Museum early (1935–56)"),
          ("museum_late", "Museum late (1965–79)"),
          ("hist1975", "Historical 1975 only (no SAFE)"), ("safe2", "SAFE only")]
PANELS = [(p, l) for p, l in PANELS if p in set(t.panel)]
BLUE, INK, MUTED = "#2a78d6", "#1f1f1f", "#8a8a8a"

sizes = [len(t[t.panel == p]) for p, _ in PANELS]
ncol = min(len(PANELS), 3); nrow = -(-len(PANELS) // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 6.8 * nrow), sharex=True, squeeze=False)
for ax in axes.flat[len(PANELS):]:
    ax.set_visible(False)
for ax, (p, lab) in zip(axes.flat, PANELS):
    g = t[t.panel == p].sort_values("mean_prob").reset_index(drop=True)
    y = np.arange(len(g))
    ax.hlines(y, g.mp_lo, g.mp_hi, color=BLUE, lw=3, alpha=0.35)
    ax.plot(g.mean_prob, y, "o", color=BLUE, ms=7)
    ax.plot(g.frac_called, y, "o", mfc="white", mec=MUTED, mew=1.4, ms=6)
    med = g.mean_prob.median()
    ax.axvline(med, color=BLUE, lw=1, ls=":")
    ax.text(med, len(g) - 0.4, f" median {med:.2f}", color=INK, fontsize=8, va="bottom")
    ax.set_yticks(y, g["sample"], fontsize=9)
    ax.set_ylim(-0.7, len(g) + 0.3)
    ax.set_title(lab, fontsize=11, loc="left")
    ax.set_xlim(0, 1)
    ax.grid(axis="x", alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.supxlabel("red wolf ancestry fraction", fontsize=10)
handles = [Line2D([0], [0], marker="o", color=BLUE, ls="none", ms=7, label="expected (mean P(red wolf), pooled model)"),
           Line2D([0], [0], color=BLUE, lw=3, alpha=0.35, label="range across fixed fractions 0.30–0.60"),
           Line2D([0], [0], marker="o", mfc="white", mec=MUTED, mew=1.4, ls="none", ms=6,
                  label="windows called red wolf at 5% FPR (pooled; conservative)")]
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.995))
fig.suptitle("Red wolf ancestry per individual under the all-fractions (pooled) ArchIE model", y=0.94, fontsize=12)
fig.tight_layout(rect=(0, 0.02, 1, 0.91))
for ext in ("png", "pdf"):
    fig.savefig(f"{S}/pooled_ancestry.{ext}", dpi=150)
print("wrote", f"{S}/pooled_ancestry.png")
