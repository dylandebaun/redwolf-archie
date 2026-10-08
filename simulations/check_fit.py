#!/usr/bin/env python3
"""
check_fit.py — how well do simulated ArchIE features match the real windows?

    python3 simulations/check_fit.py --config simulations/config/final_53g.yaml \
        --real results/rw_hist_safe_53g/data/window_stats

Compares the config's training_data.txt with the real per-window feature files
(chr*.stats.tsv.gz; chrX skipped). Reports, lower = closer to real:
  * KS distance for the 7 summary features (min dist to reference, S*,
    private SNPs, distance mean/var/skew/kurtosis) and their mean
  * L1 distance between the normalized focal-haplotype frequency spectra
    (all-derived bin excluded when the config drops it)
  * adversarial AUC: a gradient-boosting classifier trying to tell simulated
    from real windows (0.5 = indistinguishable). One random haplotype per
    window, equal window counts, 4-fold CV, mean over 3 seeds.

Writes to --out-dir (default <out_dir>/fit/):
  fit_summary.txt                       everything printed below
  train_vs_real_features.tsv            per feature: medians, 5/95 %, means, KS
  train_vs_real_sfs.tsv                 frequency spectrum, mean sites per bin
  train_vs_real_distance_by_rank.tsv    mean haplotype distance by sorted rank
  adversarial_auc.tsv                   AUC per feature group and seed
  train_vs_real_{features,sfs,distance}.{pdf,png}   plots (needs matplotlib)
"""
import argparse, glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import ks_2samp
try:
    from sklearn.ensemble import HistGradientBoostingClassifier
except ImportError:          # scikit-learn < 1.0
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa: F401
    from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import cross_val_score
import yaml

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--config", required=True)
ap.add_argument("--real", required=True, help="window_stats dir with chr*.stats.tsv.gz")
ap.add_argument("--train", help="override training_data.txt path")
ap.add_argument("--max-windows", type=int, default=2000, help="windows per side for the AUC")
ap.add_argument("--out-dir", help="where to write tables and plots (default <out_dir>/fit)")
args = ap.parse_args()

cfg = yaml.safe_load(open(args.config))
meta = json.load(open(os.path.join(cfg["out_dir"], "sim_meta.json")))
n = meta["n_focal_haps"]; nf = meta["feature_count"]
drop = bool(meta.get("drop_all_derived_bin", 0))
out_dir = args.out_dir or os.path.join(cfg["out_dir"], "fit")
os.makedirs(out_dir, exist_ok=True)
_log = []


def say(msg=""):
    print(msg); _log.append(msg)

files = sorted(f for f in glob.glob(os.path.join(args.real, "chr*.stats.tsv.gz")) if "chrX" not in f)
if not files:
    sys.exit("no autosomal chr*.stats.tsv.gz in " + args.real)
real = pd.concat([pd.read_csv(f, sep="\t", header=None) for f in files], ignore_index=True)
if real.shape[1] != nf + 4:
    sys.exit("real files have %d columns, expected %d (n_focal mismatch?)" % (real.shape[1], nf + 4))
rgroup = (real[nf + 1].astype(str) + ":" + real[nf + 2].astype(str)).values
real = real.iloc[:, :nf].to_numpy(float)

sim = pd.read_csv(args.train or os.path.join(cfg["out_dir"], "training_data.txt"), sep=r"\s+", header=None)
sim = sim.iloc[:, :nf].to_numpy(float)
sgroup = np.repeat(np.arange(len(sim) // n + 1), n)[:len(sim)]     # n rows per replicate
ok = np.isfinite(sim).all(1); sim, sgroup = sim[ok], sgroup[ok]

F = {"min_dist_ref": 2*n+6, "s_star": 2*n+7, "private_snps": 2*n+8, "dist_mean": 2*n+2,
     "dist_var": 2*n+3, "dist_skew": 2*n+4, "dist_kurt": 2*n+5}
say("windows: sim %d, real %d (autosomes)" % (len(np.unique(sgroup)), len(np.unique(rgroup))))
say("\nper-feature KS (sim median / real median):")
ks, rows = [], []
for nm, c in F.items():
    r = real[:, c-1][np.isfinite(real[:, c-1])]; s_ = sim[:, c-1]
    ks.append(ks_2samp(s_, r).statistic)
    say("  %-14s KS %.3f   (%.4g / %.4g)" % (nm, ks[-1], np.median(s_), np.median(r)))
    rows.append(dict(feature=nm, sim_median=np.median(s_), real_median=np.median(r),
                     sim_q05=np.quantile(s_, .05), real_q05=np.quantile(r, .05),
                     sim_q95=np.quantile(s_, .95), real_q95=np.quantile(r, .95),
                     sim_mean=s_.mean(), real_mean=r.mean(), ks=ks[-1]))
say("  mean KS %.3f" % np.mean(ks))
pd.DataFrame(rows).to_csv(os.path.join(out_dir, "train_vs_real_features.tsv"), sep="\t", index=False, float_format="%.6g")

bins = [c for c in range(1, n + 1) if not (drop and c == n)]
ss = np.nanmean(sim[:, bins], 0); rs = np.nanmean(real[:, bins], 0)
say("\nfrequency-spectrum L1 (%s): %.3f" % ("all-derived bin excluded" if drop else "all bins",
    np.abs(ss/ss.sum() - rs/rs.sum()).sum()))
all_s = np.nanmean(sim[:, 1:n+1], 0); all_r = np.nanmean(real[:, 1:n+1], 0)
pd.DataFrame(dict(derived_count=np.arange(1, n+1), sim_mean_sites=all_s, real_mean_sites=all_r,
                  used_by_model=[not (drop and c == n) for c in range(1, n+1)])
             ).to_csv(os.path.join(out_dir, "train_vs_real_sfs.tsv"), sep="\t", index=False, float_format="%.6g")
rank = list(range(n+2, 2*n+1))                     # sorted distances, self (rank 1) dropped
pd.DataFrame(dict(rank=np.arange(1, n), sim_mean=np.nanmean(sim[:, rank], 0), real_mean=np.nanmean(real[:, rank], 0))
             ).to_csv(os.path.join(out_dir, "train_vs_real_distance_by_rank.tsv"), sep="\t", index=False, float_format="%.6g")

feat = [c for c in range(1, nf) if not (drop and c == n)]      # col 0 (IFS bin 0) is always 0


def one_per(X, g, rng):
    idx = pd.Series(np.arange(len(g))).groupby(g).apply(lambda x: rng.choice(x.values)).values
    return X[idx]


groups = {"7 summary": [c-1 for c in F.values()], "distance vector": list(range(n+2, 2*n+1)), "all features": feat}
res = {k: [] for k in groups}
for seed in (0, 1, 2):
    rng = np.random.default_rng(seed)
    S1, R1 = one_per(sim, sgroup, rng), one_per(real, rgroup, rng)
    k = min(len(S1), len(R1), args.max_windows)
    S1 = S1[rng.choice(len(S1), k, replace=False)]; R1 = R1[rng.choice(len(R1), k, replace=False)]
    X = np.vstack([S1, R1]); y = np.r_[np.zeros(k), np.ones(k)]
    for nm, cols in groups.items():
        res[nm].append(cross_val_score(HistGradientBoostingClassifier(max_iter=150), X[:, cols], y,
                                       cv=4, scoring="roc_auc").mean())
say("\nadversarial AUC (0.5 = indistinguishable; %d windows per side):" % k)
for nm, v in res.items():
    say("  %-16s %.3f" % (nm, np.mean(v)))
pd.DataFrame([dict(features=nm, seed=i, auc=a) for nm, v in res.items() for i, a in enumerate(v)]
             ).to_csv(os.path.join(out_dir, "adversarial_auc.tsv"), sep="\t", index=False, float_format="%.4f")
with open(os.path.join(out_dir, "fit_summary.txt"), "w") as fh:
    fh.write("config: %s\nreal: %s\n\n" % (args.config, args.real) + "\n".join(_log) + "\n")

# ---- plots ------------------------------------------------------------------
try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
except ImportError:
    say("\n(matplotlib not available: tables written, plots skipped)")
    sys.exit(0)
rng = np.random.default_rng(1)
sub = lambda x, cap=60000: (lambda y: y if len(y) <= cap else rng.choice(y, cap, replace=False))(x[np.isfinite(x)])
SETS = [("real (autosomes)", real, "#222222", "-"), ("simulated (training)", sim, "#e6550d", "--")]

fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, (nm, c), kv in zip(axes.flat, F.items(), ks):
    v_all = np.concatenate([sub(m[:, c-1]) for _, m, _, _ in SETS]); hi = np.quantile(v_all, .99); lo = v_all.min()
    xs = np.linspace(lo, hi if hi > lo else lo + 1, 300)
    for lab, m, col, ls in SETS:
        v = sub(m[:, c-1]); v = v[v <= hi]
        if np.ptp(v) > 0:
            ax.plot(xs, gaussian_kde(v)(xs), color=col, ls=ls, lw=2, label=lab)
    ax.set_title("%s (KS %.3f)" % (nm, kv)); ax.legend(fontsize=7, frameon=False)
axes.flat[-1].axis("off")
fig.suptitle("%s: simulated vs real features" % cfg["experiment"]["name"]); fig.tight_layout()
for ext in ("pdf", "png"): fig.savefig(os.path.join(out_dir, "train_vs_real_features." + ext), dpi=110)

fig, axes = plt.subplots(1, 2, figsize=(13, 5)); x = np.arange(1, n+1)
for lab, m, col, ls in SETS:
    v = np.nanmean(m[:, 1:n+1], 0)
    axes[0].plot(x, v, "o" + ls, color=col, lw=2, label=lab)
    keep = np.array([not (drop and c == n) for c in range(1, n+1)])
    axes[1].plot(x[keep], v[keep] / v[keep].sum(), "o" + ls, color=col, lw=2, label=lab)
axes[0].set(xlabel="derived-allele count among focal haps", ylabel="mean sites / window", title="frequency spectrum (all bins)")
axes[1].set(xlabel="derived-allele count among focal haps", ylabel="proportion",
            title="shape, bins used by the model" + (" (all-derived bin dropped)" if drop else ""))
for ax in axes: ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"): fig.savefig(os.path.join(out_dir, "train_vs_real_sfs." + ext), dpi=110)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
v_all = np.concatenate([sub(m[:, rank].ravel()) for _, m, _, _ in SETS]); xs = np.linspace(v_all.min(), np.quantile(v_all, .99), 300)
for lab, m, col, ls in SETS:
    axes[0].plot(xs, gaussian_kde(sub(m[:, rank].ravel()))(xs), color=col, ls=ls, lw=2, label=lab)
    axes[1].plot(np.arange(1, n), np.nanmean(m[:, rank], 0), "o" + ls, color=col, lw=2, label=lab)
axes[0].set(title="pairwise haplotype distance (pooled)", xlabel="distance")
axes[1].set(title="mean distance by sorted rank", xlabel="neighbour rank (nearest -> farthest)", ylabel="mean distance")
for ax in axes: ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"): fig.savefig(os.path.join(out_dir, "train_vs_real_distance." + ext), dpi=110)
say("\nwrote tables and plots to " + out_dir)
