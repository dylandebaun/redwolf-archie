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
args = ap.parse_args()

cfg = yaml.safe_load(open(args.config))
meta = json.load(open(os.path.join(cfg["out_dir"], "sim_meta.json")))
n = meta["n_focal_haps"]; nf = meta["feature_count"]
drop = bool(meta.get("drop_all_derived_bin", 0))

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
print("windows: sim %d, real %d (autosomes)" % (len(np.unique(sgroup)), len(np.unique(rgroup))))
print("\nper-feature KS (sim median / real median):")
ks = []
for nm, c in F.items():
    r = real[:, c-1][np.isfinite(real[:, c-1])]; s = sim[:, c-1]
    ks.append(ks_2samp(s, r).statistic)
    print("  %-14s KS %.3f   (%.4g / %.4g)" % (nm, ks[-1], np.median(s), np.median(r)))
print("  mean KS %.3f" % np.mean(ks))

bins = [c for c in range(1, n + 1) if not (drop and c == n)]
ss = np.nanmean(sim[:, bins], 0); rs = np.nanmean(real[:, bins], 0)
print("\nfrequency-spectrum L1 (%s): %.3f" % ("all-derived bin excluded" if drop else "all bins",
      np.abs(ss/ss.sum() - rs/rs.sum()).sum()))

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
print("\nadversarial AUC (0.5 = indistinguishable; %d windows per side):" % k)
for nm, v in res.items():
    print("  %-16s %.3f" % (nm, np.mean(v)))
