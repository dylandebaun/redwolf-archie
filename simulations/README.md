# ArchIE training simulations

Simulates 50 kb windows under a red wolf / coyote demography, computes ArchIE's
per-haplotype features and red-wolf ancestry labels, and trains the logistic
classifier that stage 8 of `scripts/run_archie_real.sh` applies to the real
window features.

This replaces the model track in `scripts/sim/` for the 53-genome panels. It
fixes a bug in that track's templates (the ancestral population size was never
applied), refits the demography to the 53-genome data, and codes alleles the
same way as the real VCF (0 = dog reference allele). The model fit and threshold
calibration reuse `scripts/sim/train.R` and `scripts/sim/calibrate_thresholds.R`.

Two models, one per real-data panel:

| config | real-data config | focal panel | focal haps | coyote admixture |
|---|---|---|---|---|
| `config/final_53g.yaml` | `config/rw_hist_safe_53g.yaml` | 13 historical (1975) + 2 SAFE red wolves | 30 | 0.48 |
| `config/safe2_53g.yaml` | `config/rw_safe2_53g.yaml` | Cru301, Cru792 (SAFE) | 4 | 0.11 |

The sample lists and window size are read from the real-data config, so the
number of focal haplotypes always matches the real features. Everything else
(demography, rate variation, allele coding) is shared by both models.

## Quick start

Run everything from the repository root.

```bash
# one-time: ArchIE (pinned + patches, incl. simulations/calc_stats_ms_rw.py)
# and its modified ms, which tracks introgressed ancestry
bash archie/setup_archie.sh

# Python: numpy, scipy, scikit-learn, pyyaml (+ pandas for check_fit.py)
# R: data.table, jsonlite

CFG=simulations/config/final_53g.yaml            # or safe2_53g.yaml
python3 simulations/simulate_archie.py run --config $CFG --set train   # 10,000 windows
python3 simulations/simulate_archie.py run --config $CFG --set test    #  1,000 windows
bash simulations/train.sh $CFG                    # train.R + calibrate_thresholds.R

# optional: how closely do the simulated features match the real ones?
python3 simulations/check_fit.py --config $CFG --real results/rw_hist_safe_53g/data/window_stats
#   (use results/rw_safe2_53g/data/window_stats for safe2_53g)
```

To score the real data with these models, set `model_dir` in the real-data
config to the simulation `out_dir` and run stage 8:

| real-data config | `model_dir` |
|---|---|
| `config/rw_hist_safe_53g.yaml` | `simulations/results/final_53g` (currently `results/rw_hist_safe_53g/model`, the old track) |
| `config/rw_safe2_53g.yaml` | `simulations/results/safe2_53g` (already set) |

```bash
sbatch --export=ALL,CONFIG=config/rw_safe2_53g.yaml,STAGES=8 scripts/run_archie_real.sh
```

If you build `ms` by hand instead (e.g. on macOS), K&R C needs
`cc -O2 -std=gnu89 -Wno-implicit-function-declaration -Wno-return-type -o ms ms.c streec.c rand1.c -lm`
in `archie/ArchIE/msmodified/`.

Outputs go to the config's `out_dir` (`simulations/results/<name>/`):

| file | contents |
|---|---|
| `training_data.txt`, `test_data.txt` | one row per focal haplotype per simulated window |
| `sim_meta.json` | column layout read by `train.R`, `calibrate_thresholds.R`, `apply_model.R` |
| `trained_model.Rdata` | the fitted logistic model (`model`) |
| `thresholds.{summary,metrics,recommended_by_fpr}.tsv` | calibrated probability thresholds |

Use `--nreps 200` for a quick test run. Runtime on 8–12 cores is roughly
20–40 min for 10,000 final_53g windows and a few minutes for safe2_53g. Each
replicate's `ms` seeds come from `(sim.seed, train/test, replicate)`, so runs
are reproducible and the training and test sets never share seeds.

## Row layout

N = number of focal haplotypes (30 or 4). Same layout as the real feature files
from stage 6, plus labels:

| columns | content |
|---|---|
| 1 … N+1 | individual frequency spectrum: number of sites at which the focal hap is derived, by derived count 0 … N among the focal haps |
| N+2 … 2N+1 | sorted distances from the focal hap to every focal hap (first = itself, 0) |
| 2N+2 … 2N+5 | mean, variance, skew, kurtosis of those distances |
| 2N+6 | minimum distance to a reference (coyote) hap |
| 2N+7 | S\* on sites absent from the reference |
| 2N+8 | private SNPs (derived sites absent from the reference) |
| 2N+9 … 2N+11 | one-hot label: red wolf (>70 % red-wolf ancestry), not red wolf (<30 %), in between |
| 2N+12 | proportion red-wolf ancestry |

Column N+1 (derived in **all** N focal haps) is written but **not used by the
classifier**: `scripts/sim/train.R` and `calibrate_thresholds.R` drop it when
`sim_meta.json` has `drop_all_derived_bin: 1` (set from
`features.drop_all_derived_bin` in the config; see "Dropped feature" below).
`apply_model.R` then ignores that column automatically, because the fitted
model has no term for it.

## The demography

`templates/rw_collapse_outgroup.tmpl`, filled in by `simulate_archie.py`. Sizes
are in individuals, times in generations (the script converts them to `ms`
units with N0 = 60,000).

```
                           outgroup (dog/wolf, 120k)
                                 |
   ------------------------------+  t = 50,000  (outgroup joins)
   ancestor of red wolf + coyote (84k; stays 84k above 50,000)
   ---------+-------------------  t = 18,750
            |                  \
   red wolf (15k)            coyote (60k)
            |  t = 100: collapse to 2k
            |  t = 24: unadmixed red-wolf ghost merges in
            |  t = 19: 48 % (final) / 11 % (safe2) coyote ancestry enters
   red wolf (2k) -- sampled
```

| parameter | value | source / reasoning |
|---|---|---|
| mutation rate μ | 4 × 10⁻⁹ | vonHoldt et al. 2016 calibration; kept fixed |
| recombination rate | 8 × 10⁻⁹ /bp/gen | as before |
| coyote Ne (= N0) | 60,000 | as before |
| red wolf Ne, historical | 15,000 | fit to the real windows (both panels agree); see below |
| red wolf Ne, recent (last 100 gen) | 2,000 | fit; collapse kept at 100 gen as the more plausible timing |
| red wolf + coyote ancestor Ne | 84,000 (1.4 × N0) | fit; vonHoldt fig. S5 θ_COY_RED is ~110k |
| red wolf / coyote split | 18,750 gen | as before; within fig. S5 τ_COY_RED range |
| admixture time, ghost lineage | 19 / 19 / 24 gen | as before |
| coyote admixture fraction | 0.48 (final), 0.11 (safe2) | as before |
| outgroup split | 50,000 gen | ~150 kya at 3 yr/gen: vonHoldt 2016 fig. S5 τ_COY_ANC (coyote vs wolf/dog) |
| outgroup Ne | 120,000 | fig. S5 model 1 θ_NA_WOLF_ANC; outgroup size barely affects the fit |

The outgroup is not part of the red wolf demography; it only exists to mimic
how the real data were coded (next section). No gray-wolf ancestry is simulated.

## Matching how the real data were produced

1. **Window-to-window rate variation.** Each window's θ is the mean
   (4·N0·μ·L = 48) times a Gamma(6, 1/6) draw (CV ≈ 0.4). Real windows vary in
   mutation rate and callable sites; with a constant θ the simulated feature
   distributions were too narrow. Shape 6 fit better than 4 or 8.
2. **Polarization.** In the real data, 0 = the canFam3.1 (boxer) allele, not the
   ancestral allele. Each replicate simulates one outgroup haplotype as the
   stand-in for the reference genome and recodes each site so that 0 = its
   allele. Labels come from the true ancestry track, so the recoding changes the
   features but not the truth.

Each replicate simulates exactly what ArchIE needs: the focal haplotypes, the
reference coyote haplotypes, the one-haplotype red-wolf ghost, plus the single
reference-genome haplotype used for the coding above. The real VCF also holds
sites that vary only in the other genomes of the 53-genome set and is filtered
to MAF ≥ 1 % across them. Neither is simulated: those genomes are not a clean
outgroup (most of the 12 LSU museum samples are red wolves or red wolf hybrids),
and with the S\* position fix (below) sites the focal haplotypes don't vary at
no longer change any feature.

## Dropped feature

The frequency-spectrum bin "derived in all N focal haplotypes" can't be matched:
the simulations produce ~2.5× too many such sites (final: ~23 vs 8.75 per
window). The cause is the dog lineage's own mutations (every wild haplotype
looks derived after recoding); real red wolves share part of that history
through ancestry this model does not include. Sites derived in every focal hap
do not differ between haplotypes, so they carry no information about which
haplotype has which ancestry. Dropping the column halves the frequency-spectrum
mismatch (L1 0.32 → 0.14 for final, 0.22 → 0.09 for safe2).

## How well the simulations match the real data

Real = autosomal windows from stage 6 of the 53-genome data. Lower is better
for every number; adversarial AUC = how well a classifier separates simulated
from real windows (0.5 = indistinguishable; one haplotype per window, matched
window counts, 3 seeds). safe2: full run (10,000 training windows) against real
features recomputed with the S\* fix. final: pending the stage 6 re-run with the
fix; the value shown is from the run before the fix (S\* misaligned on both
sides).

| | final: previous model | final: this model | safe2: previous | safe2: this model |
|---|---|---|---|---|
| mean KS, 7 summary features | 0.200 | (0.048, pre-fix) | 0.240 | 0.119 |
| S\* KS | 0.169 | (0.045, pre-fix) | 0.241 | 0.107 |
| frequency-spectrum L1 (bin dropped) | 0.136 | (0.166, pre-fix) | 0.116 | 0.057 |
| adversarial AUC, all features (bin dropped) | 0.93 | (0.88, pre-fix) | 0.94 | 0.815 |
| model test AUC (simulated test set) | | (0.961, pre-fix) | | 0.958 |

"Previous model" = the models used in the original analysis
(`rw_3pop_collapse.tmpl`, historical 10k → recent 5k at 100 generations), run
with the ancestral-size bug below.

## What was tested and why the settings are what they are

All fits were scored against the real autosomal windows of both panels.

- **Ancestral-size bug (fixed here; still present in `scripts/sim/ms_templates/`,
  kept unchanged so that track reproduces the original analysis).** The old
  templates set the ancestral size with
  `-en T_ANC 1 SIZE_ANC` immediately before `-ej T_ANC 1 3`. After the join
  pop 1 is empty, so `size_anc: 1.66` was silently ignored and the ancestor kept
  size 1 × N0. Simulations with and without it were byte-identical. The size
  is now set on pop 3 after the join. With the fix the simulations gain ~45 %
  more segregating sites, which closes most of the diversity gap that the old
  model had (real windows had ~25–30 % more pairwise diversity).
- **S\* position bug (fixed in both ArchIE patches).** S\* scoring dropped the
  sites where the focal haplotype is ancestral from the genotypes but not from
  the positions. It surfaced because simulated S\* changed with the number of
  simulated haplotypes (safe2 median 69k with 8 haplotypes vs 47k with 107)
  although the genealogy of the focal haplotypes doesn't; with the fix both give
  ~90k. All fits above that involve S\*, and the earlier comparisons of a
  simulated MAF filter / ascertainment panel, were made with the bug; the real
  safe2 S\* median moves from 39k to 70k with the fix.
- **Mutation rate × 1.35 or all sizes × 1.35** fit well only under the bug
  (they were compensating for it); with the fixed template both overshoot.
- **Ancestral size**: 1.66 × N0 overshoots divergence from the coyote; 1.4 fits
  both panels.
- **Older red wolf / coyote split** (× 1.35, × 1.7): worse.
- **Collapse vs none**: a 10k → 5k collapse is indistinguishable from no
  collapse; a strong collapse (→ 2k) helps both panels. Historical size 15k
  beats 10k (distances too low) and 30k (distance skew/kurtosis break). 100 vs
  200 generations makes little difference.
- **Outgroup**: a closer split (~20k gen) fit the spectrum slightly better but
  contradicts the published ~150 kya split; with the published split the
  outgroup size has no effect, and the remaining spectrum mismatch is handled by
  dropping the all-derived bin.
- **Gray-wolf ancestry pulse into red wolves**: reproduced the top of the
  spectrum but is not in this model by choice.
- **MAF ≥ 1 % filter / 53-genome ascertainment panel**: tested (before the S\*
  fix) with a simulated 106-haplotype panel; small effect on everything except
  S\*, and the S\* effect turned out to come from the position bug. Not
  included: the extra genomes aren't a clean outgroup and, with the fix, sites
  the focal haplotypes don't vary at don't change the features. Turning off the
  dog-reference coding mainly worsens S\* for the 15-wolf panel (KS 0.09 → 0.14,
  pre-fix).
- **Phasing switch errors / genotype errors** added to the simulated haplotypes:
  they matched the real rate of identical haplotypes but made every other
  measure worse or no better, so they are not included.
- **Removing rw75114** (carries long dog-like blocks on chr01/chr06, ~29 % of
  the near-fixed dog alleles): barely changed the mismatch, so the sample stays.

## Known limitations

- **chrX** is not modeled (autosomal Ne, recombination and sample sizes).
  5 of the 13 historical wolves and the ALcoy reference are male. Treat chrX
  calls with caution, or exclude non-PAR chrX (> ~6.6 Mb).
- **SAFE inbreeding.** The SAFE wolves have long runs of homozygosity (many
  near-identical haplotype pairs within an individual) that a random-mating
  model does not produce; safe2's distance variance/skew fit worst (KS ~0.2).
- **Simulated windows are still distinguishable** from real ones (AUC ~0.83),
  mostly through the haplotype-distance vector. Errors and the demographic
  parameters above do not remove it; individual structure (inbreeding,
  relatedness) is the likely remaining cause.
- **Dog-like segments.** ArchIE only knows "red wolf vs coyote"; segments of
  other ancestry (e.g. rw75114 on chr01/chr06) may be called either way.

## Files

| file | role |
|---|---|
| `config/*.yaml` | demography, simulation and polarization settings per model |
| `templates/rw_collapse_outgroup.tmpl` | `ms` command template |
| `simulate_archie.py` | `meta` (writes `sim_meta.json`) and `run` (simulates a train/test set); uses ArchIE's `simulations/calc_stats_ms_rw.py` from `archie/setup_archie.sh` |
| `train.sh` | runs `scripts/sim/train.R` + `scripts/sim/calibrate_thresholds.R` on a config's output |
| `check_fit.py` | simulated vs real feature comparison (KS, spectrum, adversarial AUC) |
