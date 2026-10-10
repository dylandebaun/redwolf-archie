#!/usr/bin/env python3
"""Build an "all admixture fractions" training/test set per panel from the grid.

The grid models share simulation seeds, so replicate i has the same random
draws at every fraction. To avoid near-duplicate windows, the replicate
indices are split into equal consecutive blocks, one per fraction: with 8
fractions and 10,000 training replicates, replicates 0-1249 come from f=0.10,
1250-2499 from f=0.20, ... (likewise for the 1,000 test replicates). The
result is a 10,000-window training set spread evenly over all fractions.

    python3 sensitivity/admixture/pool.py            # from the repository root
    python3 sensitivity/admixture/pool.py safe2      # one panel only

Writes per panel sensitivity/admixture/results/<panel>_pooled/{training,test}_data.txt,
configs/<panel>_pooled.yaml (for train.sh / check_fit.py), and pooled.tsv
(panel, config, real window_stats, model_dir) for submit/summarize.
"""
import json, os, sys
import pandas as pd
import yaml

H = "sensitivity/admixture"
grid = pd.read_csv(f"{H}/grid.tsv", sep="\t", dtype={"fraction": float})
only = set(sys.argv[1:])
# keep rows of panels that are not rebuilt this time
rows = ["panel\tconfig\treal_window_stats\tmodel_dir"]
if only and os.path.exists(f"{H}/pooled.tsv"):
    rows += [l for l in open(f"{H}/pooled.tsv").read().splitlines()[1:] if l.split("\t")[0] not in only]

for panel, g in grid.groupby("panel", sort=False):
    if only and panel not in only:
        continue
    g = g.sort_values("fraction").reset_index(drop=True)
    base_cfg = yaml.safe_load(open(g.config[0]))
    meta = json.load(open(os.path.join(g.model_dir[0], "sim_meta.json")))
    n = int(meta["n_focal_haps"])
    tag = f"{panel}_pooled"
    out = f"{H}/results/{tag}"
    os.makedirs(out, exist_ok=True)

    for kind, nreps_key in (("training", "nreps_train"), ("test", "nreps_test")):
        nreps = int(meta[nreps_key]); k = len(g); block = nreps // k
        with open(f"{out}/{kind}_data.txt", "w") as fo:
            for i, r in g.iterrows():
                lines = open(os.path.join(r.model_dir, f"{kind}_data.txt")).read().splitlines()
                if len(lines) != nreps * n:
                    raise SystemExit(f"{r.model_dir}/{kind}_data.txt has {len(lines)} rows, expected {nreps}*{n}; "
                                     "replicate blocks would not line up")
                a, b = i * block * n, (i + 1) * block * n
                fo.write("\n".join(lines[a:b]) + "\n")
        print(f"{tag} {kind}: {k} fractions x {block} replicates x {n} haps = {k * block * n} rows")

    cfg = dict(base_cfg)
    cfg["experiment"] = {"name": tag}
    cfg["out_dir"] = out
    cfg["demography"] = dict(base_cfg["demography"])
    cfg["demography"]["admixture_fraction"] = "pooled: " + ", ".join(f"{f:.2f}" for f in g.fraction)
    path = f"{H}/configs/{tag}.yaml"
    with open(path, "w") as fh:
        fh.write(f"# {tag}: training/test data pooled by {H}/pool.py from the fixed-fraction grid\n"
                 "# (equal replicate blocks per fraction). Not meant for simulate_archie.py run;\n"
                 "# used by train.sh (meta) and check_fit.py.\n")
        yaml.safe_dump(cfg, fh, sort_keys=False)
    rows.append(f"{panel}\t{path}\t{g.real_window_stats[0]}\t{out}")

open(f"{H}/pooled.tsv", "w").write("\n".join(rows) + "\n")
print(f"wrote {H}/pooled.tsv")
