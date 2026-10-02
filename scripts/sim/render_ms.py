#!/usr/bin/env python3
"""
render_ms.py — turn experiment.yaml + an ms template into a ready-to-run ms.sh.

Usage:
  python3 scripts/sim/render_ms.py --config config/x.yaml --out <model_dir>/ms.sh
  (run from the repository root; sim.ms_bin and the template path are repo-relative)

The config is the single source of truth. Counts derive from the sample lists,
times are given in generations and scaled by 4*Ne, and the user supplies the
ADMIXED (introgressed) fraction f — the template's "staying" proportion is 1-f.
"""

import argparse
import os
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML not found. `pip install pyyaml` (or conda install pyyaml).")


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def derive(cfg):
    """Compute every value the template needs from the config."""
    sim = cfg["sim"]
    demo = cfg["demography"]

    Ne = float(sim["Ne"])
    mu = float(sim["mu"])
    rho_rate = float(sim["recomb_rate"])
    L = int(cfg["window_size"])

    # --- master quantities (counts) ---
    n_focal = 2 * len(cfg["admixed_samples"])
    n_ref = 2 * len(cfg["reference_samples"])
    n_total = n_focal + 1 + n_ref  # +1 for the always-singleton RWEXTINCT pop

    # --- assert against any hand-written values, if present ---
    for key, val in (("n_focal_haps", n_focal), ("n_ref_haps", n_ref)):
        if key in cfg and int(cfg[key]) != val:
            sys.exit(
                f"Config mismatch: {key}={cfg[key]} but sample lists imply {val}. "
                "Fix the list or remove the hand-written value."
            )

    # --- ms-scaled mutation / recombination ---
    theta = 4 * Ne * mu * L
    rho = 4 * Ne * rho_rate * L

    # --- admixture: user gives introgressed fraction f; ms -es keeps (1-f) ---
    f = float(demo["admixture_fraction"])
    if not 0.0 <= f <= 1.0:
        sys.exit(f"admixture_fraction must be in [0,1], got {f}")
    admix_stay = 1.0 - f

    # --- times: generations -> ms coalescent units (4*Ne) ---
    tg = demo["times_in_generations"]
    scale = 4 * Ne
    t = {k: tg[k] / scale for k in ("t_rwextinct", "t_admix", "t_rwmerge",
                                    "t_collapse", "t_anc")}

    # modified ms binary (ArchIE msmodified); absolute so ms.sh runs from any dir
    ms_bin = os.path.abspath(sim.get("ms_bin", "archie/ArchIE/msmodified/ms"))

    return {
        "MS_BIN": ms_bin,
        "N_TOTAL": n_total,
        "N_FOCAL": n_focal,
        "N_REF": n_ref,
        "L": L,
        "THETA": theta,
        "RHO": rho,
        "SIZE_HYBRID_RECENT": demo["size_hybrid_recent"],
        "SIZE_HYBRID_HIST": demo["size_hybrid_hist"],
        "SIZE_RWEXTINCT": demo["size_rwextinct"],
        "SIZE_ANC": demo["size_anc"],
        "ADMIX_STAY": admix_stay,
        "T_RWEXTINCT": t["t_rwextinct"],
        "T_ADMIX": t["t_admix"],
        "T_RWMERGE": t["t_rwmerge"],
        "T_COLLAPSE": t["t_collapse"],
        "T_ANC": t["t_anc"],
    }


def render(template_text, values):
    out = template_text
    for key, val in values.items():
        token = "{{" + key + "}}"
        if token not in out:
            sys.exit(f"Template is missing placeholder {token}")
        out = out.replace(token, repr(val) if isinstance(val, float) else str(val))
    leftover = [w for w in out.split() if w.startswith("{{")]
    if leftover:
        sys.exit(f"Unsubstituted placeholders remain: {leftover}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True, help="output ms.sh path")
    ap.add_argument("--template", help="override sim.demography_template")
    args = ap.parse_args()

    cfg = load_config(args.config)
    tmpl_path = args.template or cfg["sim"]["demography_template"]
    with open(tmpl_path) as fh:
        template_text = fh.read().strip()

    values = derive(cfg)
    ms_cmd = render(template_text, values)

    # write a self-contained ms.sh: no $1/$2 args needed anymore
    script = (
        "set -eo pipefail\n"
        f"# rendered from {args.config} + {tmpl_path} — do not edit by hand\n"
        f"# N_total={values['N_TOTAL']} (focal={values['N_FOCAL']}, "
        f"ref={values['N_REF']}, +1 rwextinct)\n"
        f"{ms_cmd}\n"
    )
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(script)

    # echo the key derived numbers so the user can eyeball them
    print(f"Wrote {args.out}")
    print(f"  N_total={values['N_TOTAL']}  N_focal={values['N_FOCAL']}  N_ref={values['N_REF']}")
    print(f"  theta={values['THETA']:.6g}  rho={values['RHO']:.6g}")
    print(f"  admix_stay(1-f)={values['ADMIX_STAY']:.6g}")
    print("  times(ms units): " + ", ".join(
        f"{k}={values[k]:.6g}" for k in
        ("T_RWEXTINCT", "T_ADMIX", "T_RWMERGE", "T_COLLAPSE", "T_ANC")))


if __name__ == "__main__":
    main()
