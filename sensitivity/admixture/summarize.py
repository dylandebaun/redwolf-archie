#!/usr/bin/env python3
"""Summarise the admixture-fraction sensitivity grid (run from the repo root).

For every panel x fraction in grid.tsv:
  model    test AUC, best-F1, threshold and recall at 5% FPR   (thresholds.*.tsv)
  fit      mean KS (7 summary features), SFS L1, adversarial AUC (fit/fit_summary.txt)
  calls    per-individual fraction of haplotype-windows called red wolf at 5% FPR
  pooled   if pooled.tsv exists (pool.py), the "all fractions" model per panel is
           reported as an extra row/point labelled "all"
  vs 0.48  per panel, against the f = 0.48 model: Jaccard of called
           haplotype-windows, share of windows whose call flips, Spearman of the
           individuals' called fractions

Writes to sensitivity/admixture/summary/:
  models.tsv, per_sample.tsv, vs_baseline.tsv, table.md
  calls_vs_fraction.{png,pdf}, fit_vs_fraction.{png,pdf}
"""
import os, re
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "sensitivity/admixture"; OUT = f"{H}/summary"; BASE_F = 0.48; FPR = "0.05"
os.makedirs(OUT, exist_ok=True)
grid = pd.read_csv(f"{H}/grid.tsv", sep="\t", dtype={"fraction": float})
POOL_X = 0.9                                    # plot position of the pooled model
if os.path.exists(f"{H}/pooled.tsv"):
    pooled = pd.read_csv(f"{H}/pooled.tsv", sep="\t").assign(fraction=np.nan)
    grid = pd.concat([grid, pooled[grid.columns]], ignore_index=True)
PANEL_LABEL = {"historical": "Historical 1975 + SAFE", "museum_early": "Museum early (1935–56)",
               "museum_late": "Museum late (1965–79)",
               "hist1975": "Historical 1975 only (no SAFE)", "safe2": "SAFE only (Cru301, Cru792)"}

models, persample, calls = [], [], {}
for r in grid.itertuples():
    tag = os.path.basename(r.config)[:-5]; d = r.model_dir
    row = dict(panel=r.panel, fraction=r.fraction, tag=tag, pooled=tag.endswith("_pooled"))
    try:
        s = pd.read_csv(f"{d}/thresholds.summary.tsv", sep="\t").iloc[0]
        t = pd.read_csv(f"{d}/thresholds.recommended_by_fpr.tsv", sep="\t")
        t5 = t[np.isclose(t.target_fpr, float(FPR))].iloc[0]
        row.update(auc=s.auc, best_f1=s.best_f1_value, thr_fpr05=t5.threshold, recall_fpr05=t5.recall)
    except Exception as e:
        row["note"] = f"no model ({type(e).__name__})"
    try:
        txt = open(f"{d}/fit/fit_summary.txt").read()
        row["mean_ks"] = float(re.search(r"mean KS ([\d.]+)", txt).group(1))
        row["sfs_l1"] = float(re.search(r"frequency-spectrum L1 \([^)]*\): ([\d.]+)", txt).group(1))
        row["adv_auc_all"] = float(re.search(r"all features\s+([\d.]+)", txt).group(1))
    except Exception:
        row["note"] = (row.get("note", "") + " no fit").strip()
    try:
        p = pd.read_csv(f"{d}/calls/{tag}.per_sample.tsv", sep="\t", dtype={"sample": str})
        for q in p.to_dict("records"):          # (itertuples renames "call_fpr0.05")
            persample.append(dict(panel=r.panel, fraction=r.fraction, tag=tag, sample=q["sample"],
                                  frac_called=q[f"call_fpr{FPR}"], mean_prob=q["mean_prob"]))
        row.update(called_min=p[f"call_fpr{FPR}"].min(), called_median=p[f"call_fpr{FPR}"].median(),
                   called_max=p[f"call_fpr{FPR}"].max())
        w = pd.read_csv(f"{d}/calls/{tag}.windows.tsv.gz", sep="\t", usecols=["sample", "hap", "chrom", "start", f"call_fpr{FPR}"],
                        dtype={"sample": str})
        calls[(r.panel, tag)] = w
    except Exception:
        row["note"] = (row.get("note", "") + " no calls").strip()
    models.append(row)
models = pd.DataFrame(models)
persample = pd.DataFrame(persample)

# ---- agreement with the f = 0.48 model
vs = []
for panel in grid.panel.unique():
    btag = f"{panel}_f{BASE_F:.2f}"
    base = calls.get((panel, btag))
    if base is None:
        continue
    b = base[f"call_fpr{FPR}"].to_numpy().astype(bool)
    bps = persample[persample.tag == btag].set_index("sample").frac_called
    for (p2, t2), w in calls.items():
        if p2 != panel:
            continue
        same = (w[["sample", "hap", "chrom", "start"]].to_numpy() == base[["sample", "hap", "chrom", "start"]].to_numpy()).all()
        if not same:
            raise SystemExit(f"window keys differ between {t2} and {btag}")
        c = w[f"call_fpr{FPR}"].to_numpy().astype(bool)
        ps = persample[persample.tag == t2].set_index("sample").frac_called
        vs.append(dict(tag=t2, jaccard_vs_048=(b & c).sum() / max((b | c).sum(), 1),
                       flipped_vs_048=(b != c).mean(),
                       spearman_samples_vs_048=spearmanr(bps.loc[ps.index], ps).statistic))
vs = pd.DataFrame(vs)
models = models.merge(vs, on="tag", how="left")
models.to_csv(f"{OUT}/models.tsv", sep="\t", index=False, float_format="%.4f")
persample.to_csv(f"{OUT}/per_sample.tsv", sep="\t", index=False, float_format="%.4f")
vs.to_csv(f"{OUT}/vs_baseline.tsv", sep="\t", index=False, float_format="%.4f")

# ---- expected red wolf ancestry under the pooled model, with the range across
#      fixed fractions 0.30-0.60 (read by plot_pooled_ancestry.py)
if "tag" in persample:
    pool = persample[persample.tag.str.endswith("_pooled")].set_index(["panel", "sample"])
    fx = persample[persample.fraction.between(0.295, 0.605)]
    rng = fx.groupby(["panel", "sample"]).agg(mp_lo=("mean_prob", "min"), mp_hi=("mean_prob", "max"),
                                              fc_lo=("frac_called", "min"), fc_hi=("frac_called", "max"))
    pool[["mean_prob", "frac_called"]].join(rng).reset_index().to_csv(
        f"{OUT}/pooled_ancestry.tsv", sep="\t", index=False, float_format="%.4f")

# ---- markdown table in the layout of the earlier comparison
lines = []
for panel, g in models.groupby("panel", sort=False):
    lines += [f"### {PANEL_LABEL.get(panel, panel)}", "",
              "| admixture fraction | model AUC | recall @5% FPR | mean KS / adversarial AUC | SFS L1 | windows called red wolf (5% FPR) | Jaccard vs 0.48 | ranking vs 0.48 (Spearman) |",
              "|---|---|---|---|---|---|---|---|"]
    for r in g.sort_values("fraction", na_position="last").itertuples():
        fmt = lambda v, f="{:.3f}": "–" if pd.isna(v) else f.format(v)
        called = "–" if pd.isna(getattr(r, "called_min", np.nan)) else f"{r.called_min*100:.0f}–{r.called_max*100:.0f}% (median {r.called_median*100:.0f}%)"
        flab = "all (pooled 0.10–0.80)" if r.pooled else (f"**{r.fraction:.2f}**" if np.isclose(r.fraction, BASE_F) else f"{r.fraction:.2f}")
        lines.append(f"| {flab} "
                     f"| {fmt(getattr(r, 'auc', np.nan))} | {fmt(getattr(r, 'recall_fpr05', np.nan), '{:.2f}')} "
                     f"| {fmt(getattr(r, 'mean_ks', np.nan))} / {fmt(getattr(r, 'adv_auc_all', np.nan), '{:.2f}')} "
                     f"| {fmt(getattr(r, 'sfs_l1', np.nan))} | {called} "
                     f"| {fmt(getattr(r, 'jaccard_vs_048', np.nan), '{:.2f}')} | {fmt(getattr(r, 'spearman_samples_vs_048', np.nan), '{:.2f}')} |")
    lines.append("")
open(f"{OUT}/table.md", "w").write("\n".join(lines))

# ---- plots
INK, MUTED, ACCENT = "#1f1f1f", "#9a9a9a", "#2a78d6"
panels = list(grid.panel.unique())
fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.6), sharey=True)
persample["x"] = persample.fraction.fillna(POOL_X)
models["x"] = models.fraction.fillna(POOL_X)
XT = sorted(models.x.dropna().unique())
XL = ["all" if np.isclose(x, POOL_X) else f"{x:.2f}".rstrip("0").rstrip(".") for x in XT]
for ax, panel in zip(np.atleast_1d(axes), panels):
    g = persample[persample.panel == panel]
    for s, x in g.groupby("sample"):
        fx = x[x.fraction.notna()].sort_values("fraction")
        ax.plot(fx.fraction, fx.frac_called * 100, color=MUTED, lw=1, alpha=0.8)
        px = x[x.fraction.isna()]
        ax.plot(px.x, px.frac_called * 100, ls="none", marker="_", ms=10, color=MUTED)
    m = g[g.fraction.notna()].groupby("fraction").frac_called.median()
    ax.plot(m.index, m.values * 100, color=ACCENT, lw=2.5, marker="o", label="median individual")
    pm = g[g.fraction.isna()].frac_called.median()
    if np.isfinite(pm):
        ax.plot([POOL_X], [pm * 100], color=ACCENT, marker="D", ms=8, ls="none", label="all-fractions model")
        ax.axvline((0.8 + POOL_X) / 2, color=MUTED, lw=0.8)
    ax.set_xticks(XT, XL)
    ax.axvline(BASE_F, color=INK, lw=0.8, ls=":")
    ax.set_title(PANEL_LABEL.get(panel, panel), fontsize=11)
    ax.set_xlabel("simulated coyote admixture fraction")
    ax.grid(alpha=0.25)
np.atleast_1d(axes)[0].set_ylabel("windows called red wolf at 5% FPR (%)")
np.atleast_1d(axes)[0].legend(frameon=False, fontsize=9)
fig.suptitle("Called red wolf ancestry vs the admixture fraction the model was trained on (grey = individuals; dotted = 0.48)", fontsize=11)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/calls_vs_fraction.{ext}", dpi=150)

metrics = [("auc", "model AUC (simulated test data)"), ("mean_ks", "mean KS vs real (lower = better fit)"),
           ("adv_auc_all", "adversarial AUC vs real (lower = better)"), ("jaccard_vs_048", "Jaccard of called windows vs 0.48")]
fig, axes = plt.subplots(1, len(metrics), figsize=(4.6 * len(metrics), 4))
styles = {"historical": ("-", "o"), "museum_early": ("--", "s"), "museum_late": (":", "^"),
          "hist1975": ("-.", "D"), "safe2": ((0, (5, 1, 1, 1, 1, 1)), "v")}
for ax, (col, lab) in zip(axes, metrics):
    for panel in panels:
        g = models[models.panel == panel].sort_values("fraction")
        if col in g:
            ls, mk = styles.get(panel, ("-", "o"))
            fx = g[g.fraction.notna()]
            ax.plot(fx.fraction, fx[col], color=INK, ls=ls, marker=mk, lw=1.6, label=PANEL_LABEL.get(panel, panel))
            px = g[g.fraction.isna()]
            ax.plot(px.x, px[col], color=INK, marker=mk, ls="none", ms=8)
    ax.axvline(BASE_F, color=MUTED, lw=0.8, ls=":")
    ax.set_title(lab, fontsize=10); ax.set_xlabel("simulated coyote admixture fraction"); ax.grid(alpha=0.25)
    ax.set_xticks(XT, XL)
    if (models.pooled).any():
        ax.axvline((0.8 + POOL_X) / 2, color=MUTED, lw=0.8)
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/fit_vs_fraction.{ext}", dpi=150)

print(open(f"{OUT}/table.md").read())
print(f"wrote {OUT}/")
