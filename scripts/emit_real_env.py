#!/usr/bin/env python3
"""
emit_real_env.py — emit shell vars for the real (application) ArchIE track,
and validate the panels against the trained model's column contract.

    eval "$(python3 emit_real_env.py --config config/x.yaml [--sim-meta sim_meta.json])"

Validation gate (fails loud, before any compute):
  * 2 * len(admixed_samples)   == n_focal_haps (if n_focal_haps is set)
  * 2 * len(reference_samples) == n_ref_haps   (if n_ref_haps is set)
  * if --sim-meta is given: its n_focal_haps must equal the config's,
    so the real feature vectors match what the model was trained on.
"""

import argparse
import json
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML not found. `pip install pyyaml`.")


def sh_quote(s):
    return "'" + str(s).replace("'", "'\\''") + "'"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--sim-meta", default=None,
                    help="optional sim_meta.json to cross-check against the model")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    # str(): YAML reads numeric IDs (e.g. museum catalogue numbers) as integers
    admixed = [str(x) for x in cfg["admixed_samples"]]
    reference = [str(x) for x in cfg["reference_samples"]]
    n_focal = 2 * len(admixed)
    n_ref = 2 * len(reference)

    # --- validation gate ---
    for key, val in (("n_focal_haps", n_focal), ("n_ref_haps", n_ref)):
        if key in cfg and int(cfg[key]) != val:
            sys.exit(f"GATE FAIL: {key}={cfg[key]} but sample list implies {val}.")

    if args.sim_meta:
        try:
            with open(args.sim_meta) as fh:
                meta = json.load(fh)
        except FileNotFoundError:
            sys.exit(f"GATE FAIL: sim_meta.json not found at {args.sim_meta}.")
        if int(meta["n_focal_haps"]) != n_focal:
            sys.exit(
                f"GATE FAIL: model trained on n_focal={meta['n_focal_haps']} "
                f"but real panel has {n_focal} haps. Feature vectors would not align."
            )

    envs = cfg.get("conda_envs", {})
    species = cfg["species"]            # {name: taxa_file}

    def emit(k, v):
        print(f"{k}={sh_quote(v)}")

    emit("EXPERIMENT_NAME", cfg.get("experiment", {}).get("name", "archie_run"))
    emit("MASTER_VCF_RAW", cfg.get("master_vcf_raw", ""))
    emit("MASTER_VCF", cfg["master_vcf"])
    emit("REFERENCE_FASTA", cfg.get("reference_fasta", ""))
    emit("TAXA_DIR", cfg["taxa_dir"])
    emit("RECOMB_MAP_DIR", cfg["recomb_map_dir"])
    emit("OUT_DIR", cfg["out_dir"])
    emit("MODEL_DIR", cfg.get("model_dir") or cfg["out_dir"] + "/model")

    emit("WINDOW_SIZE", cfg["window_size"])
    emit("STEP", cfg.get("step", cfg["window_size"]))
    emit("N_FOCAL_HAPS", n_focal)
    emit("N_REF_HAPS", n_ref)

    emit("ADMIXED_CSV", ",".join(admixed))
    emit("REFERENCE_CSV", ",".join(reference))

    emit("SPECIES_LIST", " ".join(species.keys()))
    for name, taxa in species.items():
        emit(f"SP_TAXA__{name}", taxa)
    emit("MERGE_ORDER", " ".join(cfg["merge_order"]))

    for name, taxa in (cfg.get("phasing_only_taxa") or {}).items():
        emit(f"SP_PHASE__{name}", taxa)

    emit("CHROMS", " ".join(str(c) for c in cfg["chromosomes"]))

    emit("ENV_BCFTOOLS", envs.get("bcftools", "genome-assembly"))
    emit("ENV_BEAGLE", envs.get("beagle", "redwolf"))
    emit("ENV_ARCHIE", envs.get("archie", "archie"))
    emit("ENV_PLOT", envs.get("plot", "tspop"))


if __name__ == "__main__":
    main()
