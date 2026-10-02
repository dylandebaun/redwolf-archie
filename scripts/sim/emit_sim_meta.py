#!/usr/bin/env python3
"""
emit_sim_meta.py — single source of the column contract for the sim track.

Reads experiment.yaml, writes sim_meta.json (consumed by train.R /
calibrate_thresholds.R so they never parse ms.sh again), and prints
shell-sourceable KEY=VALUE lines so create_sims.sh can:

    eval "$(python3 scripts/sim/emit_sim_meta.py --config config/x.yaml)"   # -> <model_dir>/sim_meta.json

Column contract (N = n_focal_haps):
    feature_count = 2N + 8        # cols 1 .. 2N+8
    archaic_col   = 2N + 9        # one-hot label start
    total_cols    = 2N + 12       # features + 3 labels + prop_archaic
"""

import argparse
import json
import os
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML not found. `pip install pyyaml`.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--meta-out", default=None, help="path to write sim_meta.json (default <model_dir>/sim_meta.json)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    n_focal = 2 * len(cfg["admixed_samples"])
    n_ref = 2 * len(cfg["reference_samples"])

    # assert against hand-written values if present
    for key, val in (("n_focal_haps", n_focal), ("n_ref_haps", n_ref)):
        if key in cfg and int(cfg[key]) != val:
            sys.exit(f"Config mismatch: {key}={cfg[key]} but lists imply {val}.")

    # simulations, trained model and thresholds live in model_dir
    model_dir = str(cfg.get("model_dir") or os.path.join(cfg["out_dir"], "model"))
    meta_out = args.meta_out or os.path.join(model_dir, "sim_meta.json")
    os.makedirs(os.path.dirname(meta_out) or ".", exist_ok=True)

    feature_count = 2 * n_focal + 8
    archaic_col = 2 * n_focal + 9
    meta = {
        "n_focal_haps": n_focal,
        "n_ref_haps": n_ref,
        "window_size": int(cfg["window_size"]),
        "nreps_train": int(cfg["sim"]["nreps_train"]),
        "nreps_test": int(cfg["sim"]["nreps_test"]),
        "feature_count": feature_count,
        "archaic_col": archaic_col,
        "not_archaic_col": archaic_col + 1,
        "inbetween_col": archaic_col + 2,
        "prop_archaic_col": archaic_col + 3,
        "total_cols": 2 * n_focal + 12,
        "model_dir": model_dir,
    }

    with open(meta_out, "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")

    # shell-sourceable output (uppercased keys)
    for k, v in meta.items():
        print(f"{k.upper()}={v}")

    # calibration vars — not stored in sim_meta.json (not part of column contract)
    target_fpr_str = ",".join(str(x) for x in cfg["calibration"]["target_fpr"])
    print(f"TARGET_FPR={target_fpr_str}")
    print(f"INCLUDE_INBETWEEN={str(cfg['calibration']['include_inbetween']).lower()}")
    print(f"EXPERIMENT_NAME={cfg['experiment']['name']}")
    print(f"SIM_META={meta_out}")
    print(f"ENV_ARCHIE={cfg.get('conda_envs', {}).get('archie', 'archie')}")


if __name__ == "__main__":
    main()
