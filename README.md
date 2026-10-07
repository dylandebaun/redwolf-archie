# redwolf-archie

Local-ancestry inference with [ArchIE](https://github.com/sriramlab/ArchIE) for
historical hybrid and SAFE red wolves, using one Alabama coyote as the reference panel.

Dataset: 53-genome red wolf / coyote dataset (`RedWolfProject_53genomes_geno10_maf1.vcf.gz`).

There are two tracks, both driven by one config file (`config/*.yaml`):
- **Model track** (`scripts/sim/`): simulates 50 kb windows under the red wolf
  demography with a modified `ms` program, computes ArchIE features with labels, trains
  a logistic regression, and calibrates probability thresholds for target
  false-positive rates.
- **Real-data track** (`scripts/`): cleans up the master VCF, phases each species
  with Beagle, builds the admixed (red wolf) and reference (coyote) panels,
  computes the same features in 50 kb windows, and applies the trained model.

**Updated model track (`simulations/`).** A refit of the simulations to the
53-genome data, with configs for the 15-wolf panel (`final_53g`) and the 2 SAFE
wolves (`safe2_53g`): fixed ancestral population size, red wolf 15k → 2k
100 generations ago, per-window mutation-rate variation, and the real data's
MAF filter and dog-reference allele coding simulated. See
[`simulations/README.md`](simulations/README.md); point `model_dir` at its
output to use it in stage 8. `config/rw_safe2_53g.yaml` already does.

**What the model calls.** The simulation label (`calc_stats_ms_rw.py`) is
*red-wolf-major ancestry*: a haplotype-window is positive when more than 70% of
it does not descend from the coyote/donor lineage, and negative below 30%.
`prob` and `call_fpr*` in the output = "this window is red wolf",
not "this window is introgressed".

## Setup (once)

```bash
conda env create -f envs/bcftools.yml   # env names must match conda_envs in the config
conda env create -f envs/beagle.yml
conda env create -f envs/archie.yml
conda env create -f envs/plot.yml       # optional (stage 7 plots)
bash archie/setup_archie.sh             # clone ArchIE @ 00820bc, apply patches/, build ms
```

`setup_archie.sh` produces `archie/ArchIE/`, the ArchIE folder this project
uses:

| file | change from upstream |
|---|---|
| `data/calc_stats_window_data.py` | S\* is computed on the positions of the reference-filtered sites, so positions and genotypes line up (`patches/calc_stats_window_data.patch`) |
| `simulations/calc_stats_ms_rw.py` | new: red-wolf-major labels with `--rw-thresh`/`--not-rw-thresh`, the same S\* fix, skips empty replicates (`patches/calc_stats_ms_rw.patch`) |
| `msmodified/ms` | compiled on Linux from the unmodified upstream source; upstream ships a macOS binary |

Download ArchIE. The patches apply to the pinned commit and reproduce the files exactly.
Then edit the input paths at the top of the config (`master_vcf_raw`,
`reference_fasta`, `recomb_map_dir`).

## Running

From the repository root, on the login node:

```bash
export CONDA_BASE=/path/to/miniconda3        # if conda isn't on PATH inside batch jobs
CFG=config/rw_hist_safe_53g.yaml #experiment information

# ArchIE model creation: render ms.sh -> simulate train + test (parallel) -> train + calibrate
bash scripts/sim/submit_model.sh $CFG          # prints the train job id

# real-data : prep data: subsample/phase/eigenstrat formatting -> feature extraction per chromosome -> plots -> apply ArchIE model
MODEL_JOB=<train job id> bash scripts/submit_archie_real.sh $CFG
```

The two tracks are independent until stage 8: the simulations don't use the
real data. If the model is already trained, leave out `MODEL_JOB` and stage 8
is submitted automatically. Otherwise, stage 8 waits for model track completion.
Extra `sbatch` options (partition, mail, account) go in `SBATCH_EXTRA="..."`.

Each stage writes a checkpoint, so a failed job can be resubmitted and resumes
where it stopped. To run stages by hand:
`sbatch --export=ALL,CONFIG=$CFG,STAGES="5" scripts/run_archie_real.sh`
(add `FORCE=1` to ignore checkpoints).

## Outputs

Everything goes under `out_dir` (`results/<experiment>/`; git-ignored):

| path | contents |
|---|---|
| `model/ms.sh` | rendered `ms` command for this config |
| `model/sim_meta.json` | column layout of the features (`n_focal_haps`, `feature_count`, ...) |
| `model/{training,test}_data.txt` | simulated features + labels (one row per haplotype × replicate) |
| `model/trained_model.Rdata` | fitted `glm` (`model`) |
| `model/thresholds.{summary,metrics,recommended_by_fpr}.tsv` | AUC, full threshold sweep, threshold per target FPR |
| `vcfs/` | per-species subsample, Beagle output per chromosome, merged phased VCF |
| `data/rw_coy_archie.*`, `data/by_chr/` | phased Eigenstrat panels (admixed / reference), per chromosome |
| `data/window_stats/chrNN.stats.tsv.gz` | real features: 2N+8 feature columns, then sample, chrom, start, end |
| `data/window_stats/plots/` | S\*, min. distance to the reference, private SNPs along the genome |
| `calls/<experiment>.windows.tsv.gz` | per haplotype × window: `prob`, `call_fpr0.01/0.02/0.05` |
| `calls/<experiment>.per_sample.tsv` | fraction of windows called per sample |

## Pipeline detail

### Model track (`scripts/sim/`)

| script | what |
|---|---|
| `render_ms.py` | fills `ms_templates/rw_3pop_collapse.tmpl` from `sim:` and `demography:` (times in generations → 4Ne units; θ = 4Neμ·L, ρ = 4Ner·L) |
| `emit_sim_meta.py` | derives the feature/label column layout from the panel size; writes `sim_meta.json` |
| `create_sims.sh` (+ `create_training.sh`, `create_test.sh`) | `nreps` × (`ms` → `calc_stats_ms_rw.py`) |
| `train.sh` → `train.R`, `calibrate_thresholds.R` | logistic regression on features 1..2N+8; threshold sweep on the test set |
| `submit_model.sh` | renders `ms.sh` and submits the three jobs with dependencies |

The template simulates HYBRID (the focal red wolves, N = 2 × admixed
samples), a one-haplotype extinct red wolf lineage, and COYOTE (the reference,
2 × reference samples). A ghost lineage carries a pulse of
`admixture_fraction` from COYOTE into HYBRID at `t_admix`. The real panel size
must match the model: `n_focal_haps` = 2 × `len(admixed_samples)`. The scripts
derive it from the config and check it at every step.

ms_templates/rw_3pop_collapse.tmpl models a population collapse in red wolf hybrids
ms_templates/rw_3pop.tmpl has typical demography

### Real-data track (`scripts/run_archie_real.sh`)

| stage | what | tool |
|---|---|---|
| 0 | clean up the PLINK master VCF (below) | `prepare_master_vcf.sh` |
| 1 | subsample each species panel (+ phasing-only samples) | bcftools |
| 2 | phase each species per chromosome with the dog recombination map  | Beagle 5.5 |
| 3 | concatenate chromosomes per species | bcftools |
| 4 | drop phasing-only samples, merge species | bcftools |
| 5 | VCF → phased Eigenstrat; split admixed/reference; split by chromosome | `archie_prep/*.py` |
| 6 | `calc_stats_window_data.py`, one chromosome per array task | ArchIE (patched) |
| 7 | plots | matplotlib |
| 8 | score windows with the trained model; call at calibrated thresholds | `apply_model.R` |

## Samples

`taxa/` holds one-ID-per-line sample lists, using the 53-genome IDs.

| file | samples |
|---|---|
| `redwolf_1975.txt` | 13 historical red wolves |
| `redwolf_safe.txt` | Cru301, Cru792 |
| `rw1975safe.txt` | the two above combined (15) |
| `redwolf_hybrids.txt` | 19 current Louisiana/Texas canids |
| `rwcurrentsafe.txt` | SAFE + current (21) |
| `coyote.txt` / `coyote_bama.txt` | 4 coyotes (phasing) / ALcoy (reference panel) |
| `wolf.txt` | AKwolf, WIwolf |

Compared with the earlier 63-genome VCF, the 53-genome set no longer has
CL13859, CL13869, CL14408 or the golden jackal / Caucasian outgroups. It adds
samples (1022–1025, 13750, 15308, 17018, 17033, 17679, 20904, 7006, 9044, GSHP)
that this analysis doesn't use.

To run another panel (e.g. current red wolves), copy the config and change
`experiment.name`, `out_dir`, `model_dir`, `species.redwolf` and
`admixed_samples`. A different number of admixed samples needs its own
simulations and model.

## Caveats

- **Site ascertainment.** The master VCF is filtered to ≤10% missingness and
  MAF ≥ 1% across all 53 genomes, which removes most singletons. The IFS bins,
  S\* and private-SNP counts are sensitive to rare variants, and the `ms`
  simulations are unascertained. Compared with the 63-genome run, median
  private SNPs per window dropped 25 → 20 and S\* 18.0k → 16.2k. Applying the
  same MAF filter to the simulated data would make the training data match.
- **`size_anc` is ignored by `scripts/sim/ms_templates/*.tmpl`.** They set the
  ancestral size on pop 1 just before pop 1 joins pop 3, so the root keeps
  size 1 × Ne (simulations with and without `size_anc` are identical). Kept as
  is to reproduce the original analysis; `simulations/` uses the corrected model.
- **chrX** is phased and scored like the autosomes, as in the earlier analysis.
  Males are effectively haploid there.

## Other files

- `scripts/make_beagle_maps.sh`: builds `chrNN.beagle.map` from the canFam3.1
  maps of [Auton et al. 2013](https://github.com/auton1/dog_recomb).
- `envs/*.yml`: conda environments (tool versions as used in this project).
