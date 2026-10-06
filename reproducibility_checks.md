## Reproducibility checks

- `calc_stats_window_data.py`, `calc_stats_ms_rw.py` and `ms` produced by
  `setup_archie.sh` match the versions used for the original analysis (files
  byte-identical; `ms` output identical over 40 replicates).
- `render_ms.py` with `config/rw_hist_safe_53g.yaml` reproduces the original
  `ms` command exactly.
- Retraining with `scripts/sim/train.sh` on the original simulations gives the
  same thresholds (0.911 / 0.881 / 0.806 at FPR 1 / 2 / 5%).
