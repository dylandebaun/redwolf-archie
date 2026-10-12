# Red wolf local ancestry with ArchIE: results so far

Status as of 2026-10-10. All paths are relative to the repository root.
This file summarises the ArchIE runs and what they found; methods are in
`README.md` (real data) and `simulations/README.md` (models).

**Contents**
1. Data
2. ArchIE models and calls (historical, SAFE only, museum early, museum late)
3. Museum samples: quality control and panel choice
4. Open issues and next steps
5. Files

---

## 1. Data

**Master VCF:** `RedWolfProject_53genomes_geno10_maf1.vcf.gz`
- **REF/ALT:** PLINK had put the minor allele in ALT at 2,556,745 of
  13,833,055 sites (18.5%). These were swapped back to the canFam3.1
  reference allele with genotypes recoded (`bcftools +fixref -m flip-all`).
  No strand flips were needed, all REF alleles match afterwards, and a spot
  check of 500 swapped sites showed every genotype mirrored. This matters
  because ArchIE reads allele `1` as derived.

## ArchIE models and calls

**What the model calls.** The simulation label is *red-wolf-major ancestry*:
a haplotype-window is positive when > 70% of it is not from the coyote
lineage. "Called" below means called red wolf.

### Panels and models

| panel | focal samples | haplotypes | sim config | time shift |
|---|---|---|---|---|
| historical | 13 historical (1975) + Cru301, Cru792 (SAFE) | 30 | `final_53g` | 0 (model time 0 = 1975) |
| SAFE only | Cru301, Cru792 (SAFE) | 4 | `safe2_53g` | 0 |
| museum early | 1022, 1023, 1024, 1025 (1935) + 7006 (1956) | 10 | `museum_early_53g` | 8.95 generations |
| museum late | 13750, 15308 (1965–66), 17018, 17033, 17679 (1971–73), 20904 (1979) | 12 | `museum_late_53g` | 1.0 generation |

- **Reference panel:** ALcoy (2 haplotypes) for all four.
- **Simulation set-up:** `simulations/README.md`. Red wolf Ne 15,000
  collapsing to 2,000 at 100 generations; ancestor 84,000; red wolf/coyote
  split at 18,750 generations; a dog/wolf outgroup haplotype used only to code
  alleles like the canFam3.1 reference; per-window mutation-rate variation
  (Gamma(6)); no MAF ascertainment; the "derived in all focal haplotypes"
  bin dropped from the classifier.
- **Admixture:** a coyote pulse at 19 generations before 1975 (≈ 1900,
  4-year generations). Coyote fraction per panel:
  - historical **0.48**, the published estimate;
  - museum early **0.30** and museum late **0.40**, the fractions that fit the
    real data best in the admixture sensitivity analysis
    (`sensitivity/admixture/`; previously both used 0.48 as a placeholder);
  - SAFE only **0.11**, as in the earlier SAFE-only analysis
    (`rw_safe2_11pct`), with the same timing.
- **SAFE-only real data:** phasing is shared with the historical panel (the 13
  historical wolves are phased alongside); only the panel split and features
  (stages 5–6) differ.
- **Museum time shifts:** every event time is moved closer by
  (1975 − mean collection year) / 4.
- **Data:** 10,000 training and 1,000 test windows per model.

### Model performance and fit to the real data

| | historical | SAFE only | museum early | museum late |
|---|---|---|---|---|
| coyote admixture fraction | 0.48 | 0.11 | 0.30 | 0.40 |
| model AUC (simulated test data) | 0.961 | 0.958 | 0.965 | 0.954 |
| best F1 (threshold) | 0.900 (0.478) | 0.966 (0.468) | 0.936 (0.497) | 0.902 (0.502) |
| threshold / recall at 5% FPR | 0.671 / 0.82 | 0.940 / 0.80 | 0.813 / 0.84 | 0.751 / 0.79 |
| fit to real data: mean KS (7 summary features) | 0.060 | 0.119 | 0.064 | 0.058 |
| frequency-spectrum L1 (all-derived bin dropped) | 0.166 | 0.057 | 0.094 | 0.122 |
| windows called red wolf (5% FPR) | 33–60% | 57–59% | 51–69% | 32–67% |

Fit plots and tables: `simulations/results/<model>/fit/`.

### Called red wolf ancestry per individual (fraction of haplotype-windows, 5% FPR)

| historical | | SAFE only | | museum early | | museum late | |
|---|---|---|---|---|---|---|---|
| Cru792 | 0.60 | Cru792 | 0.59 | 1022 | 0.69 | 17679 | 0.67 |
| Cru301 | 0.58 | Cru301 | 0.57 | 1023 | 0.68 | 17018 | 0.66 |
| CL5480 | 0.57 | | | 7006 | 0.55 | 17033 | 0.56 |
| rw75114 | 0.51 | | | 1025 | 0.55 | 20904 | 0.38 |
| CL5487, rw76047 | 0.46 | | | 1024 | 0.51 | 15308 | 0.32 |
| rw75157, CL5485, rw75092, rw75014 | 0.41–0.43 | | | | | 13750 | 0.32 |
| rw75154, rw75032, rw76046, CL5494 | 0.38–0.39 | | | | | | |
| CL5479 | 0.33 | | | | | | |

Calls: `results/rw_<panel>_53g/calls/` (per window and per sample).

- **Consistent ranking.** The SAFE wolves, CL5480 and rw75114 carry the most
  red wolf ancestry among the historical animals. In the late museum panel,
  17018, 17033 and 17679 (the samples closest to SAFE genetically and by
  mitogenome) have the most. The catalogued RW × coyote specimens 13750 and
  15308 have the least of any animal.
- **SAFE wolves agree across models.** Scored by the SAFE-only model (4
  haplotypes, admixture 0.11), Cru301 and Cru792 get 0.57 and 0.59; scored
  within the 15-wolf historical model (30 haplotypes, admixture 0.48), 0.58
  and 0.60.

## 3. Museum samples: quality control and panel choice

- **9044 excluded.** It's catalogued RW × dog, but genetically almost a dog:
  allele-sharing distance to GSHP 0.088, vs ≥ 0.17 to every other group. Its
  mitogenome is in the dog/wolf cluster.
- **17018 and 20904** (catalogued RW × dog) show no clear dog signal: 17018 is
  closest to SAFE, 20904 to the current LA/TX animals.
- **1013** is not in the 53-genome VCF.
- **1022** has about 15% missing genotypes (others < 2%) and low
  heterozygosity; it is kept but flagged.
- **1022–1025 and 7006** are closer to gray wolf (0.172–0.178) than to the
  1975 red wolves. That could be older red wolf ancestry or reference bias; it
  hasn't been resolved.
- **7006 moved to the early panel.** Collected in 1956, it groups with the
  1935 wolves in both the nuclear distances and the mitogenome network.

## 5. Files

| path | contents | in git |
|---|---|---|
| `simulations/results/{final,safe2,museum_early,museum_late}_53g/` | thresholds, `sim_meta.json`, `fit/` | yes |
| `results/rw_{hist_safe,safe2,museum_early,museum_late}_53g/calls/` | ArchIE calls per window and per sample | yes |
| `config/rw_museum_{early,late}_53g.yaml`, `simulations/config/museum_*` | museum panels | yes |
| `results/museum_qc/` | museum distance heatmap and tables | no |
| trained models, training/test data | `trained_model.Rdata`, `*_data.txt` | no (reproducible from configs) |
