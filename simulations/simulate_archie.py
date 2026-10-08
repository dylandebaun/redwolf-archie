#!/usr/bin/env python3
"""
simulate_archie.py — ArchIE training/test simulations for the red-wolf panels.

    python3 simulations/simulate_archie.py meta --config simulations/config/final_53g.yaml
    python3 simulations/simulate_archie.py run  --config simulations/config/final_53g.yaml --set train
    python3 simulations/simulate_archie.py run  --config simulations/config/final_53g.yaml --set test

Run from the repository root. Each replicate is one 50 kb window:

  1. draw this window's theta = 4*N0*mu*L * Gamma(k, 1/k)          (rate variation)
  2. run the modified ms (templates/rw_collapse_outgroup.tmpl): the focal red
     wolf haps (pop 1), the red-wolf ghost (2), the reference coyote haps (3)
     and one dog/wolf outgroup hap (4)
  3. labels: red-wolf-major ancestry of each focal hap from the TRUE ancestry
     track (out.ADMIXED.anc), as in calc_stats_ms_rw.py
  4. polarization: recode 0/1 so that 0 = the outgroup hap's allele, as the
     real data are coded against the canFam3.1 (boxer) reference
  5. features: calc_stats_ms_rw.py feature functions on the focal and
     reference haps

Output rows have the same layout as the old create_sims.sh output
(2N+8 features, 3 one-hot labels, prop red-wolf), so train.R and
calibrate_thresholds.R read them unchanged. Dropping the all-derived IFS bin
is done at model-fit time by scripts/sim/train.R and calibrate_thresholds.R
(see sim_meta.json: drop_all_derived_bin).
"""
import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import numpy as np

try:
    import yaml
except ImportError:
    sys.exit("PyYAML not found: pip install pyyaml")

HERE = os.path.dirname(os.path.abspath(__file__))
# feature functions + labels: ArchIE's calc_stats_ms_rw.py, created by
# archie/setup_archie.sh (patches/calc_stats_ms_rw.patch); ArchIE isn't vendored.
CALC = os.path.join(HERE, "..", "archie", "ArchIE", "simulations", "calc_stats_ms_rw.py")
if not os.path.exists(CALC):
    sys.exit("missing %s: run bash archie/setup_archie.sh first" % os.path.normpath(CALC))
_spec = importlib.util.spec_from_file_location("calc", CALC)
calc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(calc)


# --------------------------------------------------------------------- config
def load_config(path):
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    with open(cfg["real_config"]) as fh:
        real = yaml.safe_load(fh)
    cfg["n_focal"] = 2 * len(real["admixed_samples"])
    cfg["n_ref"] = 2 * len(real["reference_samples"])
    cfg["window_size"] = int(real["window_size"])
    # ms runs inside a per-replicate temp dir, so resolve paths now
    cfg["sim"]["ms_bin"] = os.path.abspath(cfg["sim"]["ms_bin"])
    cfg["sim"]["template"] = os.path.abspath(cfg["sim"]["template"])
    if not os.access(cfg["sim"]["ms_bin"], os.X_OK):
        sys.exit("ms binary not found/executable: %s (see simulations/README.md)" % cfg["sim"]["ms_bin"])
    return cfg


def sim_meta(cfg):
    n = cfg["n_focal"]
    archaic_col = 2 * n + 9
    return {
        "n_focal_haps": n,
        "n_ref_haps": cfg["n_ref"],
        "window_size": cfg["window_size"],
        "nreps_train": int(cfg["sim"]["nreps_train"]),
        "nreps_test": int(cfg["sim"]["nreps_test"]),
        "feature_count": 2 * n + 8,
        "archaic_col": archaic_col,
        "not_archaic_col": archaic_col + 1,
        "inbetween_col": archaic_col + 2,
        "prop_archaic_col": archaic_col + 3,
        "total_cols": 2 * n + 12,
        "drop_all_derived_bin": int(bool(cfg["features"]["drop_all_derived_bin"])),
        "model_dir": cfg["out_dir"],
    }


def render_template(cfg, theta):
    """ms command for one replicate; sizes -> units of N0, times -> units of 4*N0 gens."""
    s, d = cfg["sim"], cfg["demography"]
    n0, L = float(s["N0"]), cfg["window_size"]
    g = lambda t: float(t) / (4 * n0)
    sz = lambda ne: float(ne) / n0
    with open(s["template"]) as fh:
        tmpl = fh.read().strip()
    n_rw, n_coy, n_out = cfg["n_focal"], cfg["n_ref"], 1      # one outgroup hap = the "reference genome"
    vals = {
        "MS_BIN": s["ms_bin"], "N_TOTAL": n_rw + 1 + n_coy + n_out, "L": L,
        "THETA": theta, "RHO": 4 * n0 * float(s["recomb_rate"]) * L,
        "N_RW": n_rw, "N_COY": n_coy, "N_OUT": n_out,
        "SIZE_RW_RECENT": sz(d["Ne_rw_recent"]), "SIZE_RW_HIST": sz(d["Ne_rw_hist"]),
        "SIZE_RWEXTINCT": sz(d["Ne_rwextinct"]), "SIZE_ANC": sz(d["Ne_anc"]),
        "SIZE_OUTGROUP": sz(d["Ne_outgroup"]),
        "ADMIX_STAY": 1.0 - float(d["admixture_fraction"]),
        "T_COLLAPSE": g(d["t_collapse"]), "T_RWEXTINCT": g(d["t_rwextinct"]),
        "T_ADMIX": g(d["t_admix"]), "T_RWMERGE": g(d["t_rwmerge"]),
        "T_RW_COYOTE": g(d["t_rw_coyote"]), "T_OUTGROUP": g(d["t_outgroup"]),
    }
    if float(d["Ne_coyote"]) != n0:
        sys.exit("Ne_coyote must equal sim.N0 (the template fixes the coyote size at 1 x N0)")
    for k, v in vals.items():
        tmpl = tmpl.replace("{{%s}}" % k, str(v))
    if "{{" in tmpl:
        sys.exit("unfilled template token in: " + tmpl)
    return tmpl


# ------------------------------------------------------------------ replicate
def read_geno(path):
    with open(path) as fh:
        rows = [l.rstrip("\n") for l in fh]
    if not rows:
        return np.zeros((0, 0), dtype=np.int8)
    return np.array([[int(c) for c in r] for r in rows], dtype=np.int8)   # sites x haps


def one_replicate(job):
    rep, cfg, set_offset, workdir = job
    s = cfg["sim"]
    rng = np.random.default_rng([int(s["seed"]), set_offset, rep])
    L = cfg["window_size"]
    theta = 4 * float(s["N0"]) * float(s["mu"]) * L
    k = float(s.get("theta_gamma_shape", 0) or 0)
    if k > 0:
        theta *= rng.gamma(k, 1.0 / k)
    n, n_ref = cfg["n_focal"], cfg["n_ref"]

    with tempfile.TemporaryDirectory(dir=workdir) as wd:
        with open(os.path.join(wd, "seedms"), "w") as fh:          # ms seeds itself from ./seedms
            fh.write(" ".join(str(x) for x in rng.integers(1, 2**31 - 1, 3)))
        cmd = render_template(cfg, theta)
        # pipefail: an ms crash must stop the run, not leave stale/missing out.* files
        subprocess.run(["bash", "-c", "set -o pipefail; " + cmd + " | tail -n1"], cwd=wd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rw = read_geno(os.path.join(wd, "out.ADMIXED.geno"))       # pop 1
        coy = read_geno(os.path.join(wd, "out.2.geno"))            # pop 3
        out = read_geno(os.path.join(wd, "out.3.geno"))            # pop 4
        pos_all = np.array(calc.parse_snp(os.path.join(wd, "out.snp")))
        with open(os.path.join(wd, "out.ADMIXED.anc")) as fh:
            anc = [l.rstrip("\n") for l in fh]

    if rw.shape[0] == 0:
        return []
    # labels from the true ancestry track
    labels = [calc.label_rw_major([row[i] for row in anc], list(pos_all), L, 0.7, 0.3) for i in range(n)]

    if cfg["polarization"]["by_outgroup"]:
        flip = (out[:, 0] == 1)[:, None]          # the "reference genome" carries the derived allele
        rw, coy = np.where(flip, 1 - rw, rw), np.where(flip, 1 - coy, coy)
    adm, ref, pos = rw, coy, pos_all

    genotypes, ref_geno = adm.tolist(), ref.tolist()
    priv = ref.sum(1) == 0                     # sites not carried by any reference hap (S*, private)
    s_haps, ppos = adm[priv].T, pos[priv]
    rows = []
    for i in range(n):
        rows.append(calc.N_ton(genotypes, n, i) + calc.distance_vector(genotypes, i)
                    + [calc.min_dist_ref(genotypes, ref_geno, i)]
                    + [calc.s_star_ind(s_haps.copy(), ppos, i)]
                    + [int(s_haps[i].sum())] + labels[i])
    return rows


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["meta", "run"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--set", choices=["train", "test"], default="train")
    ap.add_argument("--nreps", type=int, help="override nreps_train / nreps_test")
    ap.add_argument("--workers", type=int, help="override sim.workers")
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(cfg["out_dir"], exist_ok=True)
    meta = sim_meta(cfg)
    with open(os.path.join(cfg["out_dir"], "sim_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    if args.command == "meta":
        for k, v in meta.items():
            print("%s=%s" % (k.upper(), v))
        print("TARGET_FPR=" + ",".join(str(x) for x in cfg["calibration"]["target_fpr"]))
        print("INCLUDE_INBETWEEN=" + str(cfg["calibration"]["include_inbetween"]).lower())
        print("EXPERIMENT_NAME=" + cfg["experiment"]["name"])
        return

    nreps = args.nreps or int(cfg["sim"]["nreps_%s" % args.set])
    workers = args.workers or int(cfg["sim"]["workers"])
    out_path = os.path.join(cfg["out_dir"], "%s_data.txt" % ("training" if args.set == "train" else "test"))
    print("%s: %d replicates, n_focal=%d, %d workers -> %s" % (cfg["experiment"]["name"], nreps,
          cfg["n_focal"], workers, out_path), file=sys.stderr)
    print("example ms command (theta at its mean): " + render_template(
        cfg, 4 * float(cfg["sim"]["N0"]) * float(cfg["sim"]["mu"]) * cfg["window_size"]), file=sys.stderr)

    set_offset = 0 if args.set == "train" else 1          # train and test never share seeds
    workdir = os.path.abspath(cfg["out_dir"])
    jobs = [(r, cfg, set_offset, workdir) for r in range(nreps)]
    n_rows = 0
    with Pool(workers) as pool, open(out_path, "w") as fh:
        for rows in pool.imap(one_replicate, jobs, chunksize=4):
            for row in rows:
                fh.write("\t".join(str(x) for x in row) + "\n")
            n_rows += len(rows)
    print("wrote %d rows (%d expected columns each)" % (n_rows, meta["total_cols"]), file=sys.stderr)


if __name__ == "__main__":
    main()
