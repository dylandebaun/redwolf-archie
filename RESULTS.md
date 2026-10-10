# Red wolf local ancestry with ArchIE: results so far

Status as of 2026-10-10. All paths are relative to the repository root.
This file summarises the ArchIE runs and what they found; methods are in
`README.md` (real data) and `simulations/README.md` (models).

**Contents**
1. Data
2. ArchIE models and calls (historical, museum early, museum late)
3. Museum samples: quality control and panel choice
4. Open issues and next steps
5. Files

---

## 1. Data

**Master VCF:** `RedWolfProject_53genomes_geno10_maf1.vcf.gz` (PLINK 1.9
export; ≤ 10% missingness, MAF ≥ 1% across 53 genomes). Stage 0
(`scripts/prepare_master_vcf.sh`) makes four changes:

- **Sample IDs:** PLINK `FID_IID` (`ALcoy_ALcoy`) → IID (`ALcoy`).
- **Chromosomes:** PLINK dog codes (1–38, 39 = X) → `chr01`…`chr38`, `chrX`;
  codes 40/41 dropped.
- **REF/ALT:** PLINK had put the minor allele in ALT at 2,556,745 of
  13,833,055 sites (18.5%). These were swapped back to the canFam3.1
  reference allele with genotypes recoded (`bcftools +fixref -m flip-all`).
  No strand flips were needed, all REF alleles match afterwards, and a spot
  check of 500 swapped sites showed every genotype mirrored. This matters
  because ArchIE reads allele `1` as derived.
- **Damage-type sites already absent:** like the earlier 63-genome VCF, this
  one contains no C>T or G>A sites relative to canFam3.1. That removes the
  substitutions museum-DNA damage creates, but also about a third of all SNPs.
  The simulations don't mimic it.

**Samples relative to the earlier 63-genome VCF:**
- **Missing:** CL13859, CL13869, CL14408 (current LA/TX), 1013 (museum) and
  the jackal / Caucasian outgroups.
- **Added:** 12 LSUMZ museum specimens and GSHP (German shepherd).

**S\* fix:** commit `3dbfda1` corrected a position misalignment in ArchIE's S\*
for both real data and simulations. Every result below uses the corrected
code. The earlier 63- vs 53-genome comparison
(`results/compare_63g_vs_53g/`) predates the fix and the new models; its
numbers are superseded.

## 2. ArchIE models and calls

**What the model calls.** The simulation label is *red-wolf-major ancestry*:
a haplotype-window is positive when > 70% of it is not from the coyote
lineage. "Called" below means called red wolf.

### Panels and models

| panel | focal samples | haplotypes | sim config | time shift |
|---|---|---|---|---|
| historical | 13 historical (1975) + Cru301, Cru792 (SAFE) | 30 | `final_53g` | 0 (model time 0 = 1975) |
| museum early | 1022, 1023, 1024, 1025 (1935) + 7006 (1956) | 10 | `museum_early_53g` | 8.95 generations |
| museum late | 13750, 15308 (1965–66), 17018, 17033, 17679 (1971–73), 20904 (1979) | 12 | `museum_late_53g` | 1.0 generation |

- **Reference panel:** ALcoy (2 haplotypes) for all three.
- **Simulation set-up:** `simulations/README.md`. Red wolf Ne 15,000
  collapsing to 2,000 at 100 generations; ancestor 84,000; red wolf/coyote
  split at 18,750 generations; a dog/wolf outgroup haplotype used only to code
  alleles like the canFam3.1 reference; per-window mutation-rate variation
  (Gamma(6)); no MAF ascertainment; the "derived in all focal haplotypes"
  bin dropped from the classifier.
- **Admixture:** coyote fraction 0.48 at 19 generations before 1975 (≈ 1900,
  4-year generations) for all three. 0.48 is the published estimate for the
  historical panel and **a placeholder for the museum panels**.
- **Museum time shifts:** every event time is moved closer by
  (1975 − mean collection year) / 4.
- **Data:** 10,000 training and 1,000 test windows per model.

### Model performance and fit to the real data

| | historical | museum early | museum late |
|---|---|---|---|
| model AUC (simulated test data) | 0.961 | 0.917 | 0.923 |
| best F1 (threshold) | 0.900 (0.478) | 0.843 (0.410) | 0.851 (0.444) |
| threshold / recall at 5% FPR | 0.671 / 0.82 | 0.758 / 0.65 | 0.752 / 0.67 |
| fit to real data: mean KS (7 summary features) | 0.060 | 0.097 | 0.073 |
| frequency-spectrum L1 (all-derived bin dropped) | 0.166 | 0.096 | 0.160 |
| adversarial AUC: summary / distance vector / all | 0.79 / 0.73 / 0.88 | 0.76 / 0.68 / 0.84 | 0.77 / 0.71 / 0.88 |
| windows called red wolf (5% FPR) | 33–60% | 44–60% | 27–62% |

Fit plots and tables: `simulations/results/<model>/fit/`.

- **Haplotype-distance shape** is the largest mismatch in every panel:
  simulated windows have more strongly negative skew and higher kurtosis than
  real ones.
- **Frequency spectrum:** the real data has fewer singletons than the
  simulations and more variants at high frequency among the focal
  haplotypes. Two likely causes:
  - the MAF filter and missing C>T/G>A sites, which the simulations don't
    mimic;
  - for the late panel, relatedness among the SAFE-like wolves.
- **Overall:** simulated and real windows remain distinguishable (adversarial
  AUC 0.84–0.88), as in the README's earlier assessment.

### Called red wolf ancestry per individual (fraction of haplotype-windows, 5% FPR)

| historical | | museum early | | museum late | |
|---|---|---|---|---|---|
| Cru792 | 0.60 | 1023 | 0.60 | 17679 | 0.62 |
| Cru301 | 0.58 | 1022 | 0.58 | 17018 | 0.61 |
| CL5480 | 0.57 | 7006 | 0.48 | 17033 | 0.52 |
| rw75114 | 0.51 | 1025 | 0.48 | 20904 | 0.33 |
| CL5487, rw76047 | 0.46 | 1024 | 0.44 | 13750 | 0.27 |
| rw75157, CL5485, rw75092, rw75014 | 0.41–0.43 | | | 15308 | 0.27 |
| rw75154, rw75032, rw76046, CL5494 | 0.38–0.39 | | | | |
| CL5479 | 0.33 | | | | |

Calls: `results/rw_<panel>_53g/calls/` (per window and per sample).

- **Consistent ranking.** The SAFE wolves, CL5480 and rw75114 carry the most
  red wolf ancestry among the historical animals. In the late museum panel,
  17018, 17033 and 17679 (the samples closest to SAFE genetically and by
  mitogenome) have the most. The catalogued RW × coyote specimens 13750 and
  15308 have the least of any animal.
- **Not directly comparable across panels.** Each panel has its own model and
  panel size, and the museum models use a placeholder admixture fraction.

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
- **Damage:** the museum samples' transition fraction (0.529–0.540) is in the
  range of the modern dog, wolves and coyotes, with no excess.
- **Files:** `results/museum_qc/` (distance heatmap, tables; not tracked in git).

## 4. Open issues and next steps

1. **Admixture fraction for the museum panels:** 0.48 is a placeholder; the
   amount of called red wolf ancestry depends on it, so museum percentages
   need a justified fraction or a sensitivity range.
2. **Close the simulation–real gap:**
   - simulate the missing C>T/G>A sites (drop them after the outgroup
     recoding) and the MAF filter;
   - investigate the haplotype-distance shape mismatch.
3. **chrX:** not modeled in the simulations (5 of the 13 historical wolves and
   ALcoy are male).
4. **Museum early panel:** only 5 samples (10 haplotypes), including 1022
   with 15% missing data. Its calls are the least certain.

## 5. Files

| path | contents | in git |
|---|---|---|
| `simulations/results/{final,museum_early,museum_late}_53g/` | thresholds, `sim_meta.json`, `fit/` | yes |
| `results/rw_{hist_safe,museum_early,museum_late}_53g/calls/` | ArchIE calls per window and per sample | yes |
| `config/rw_museum_{early,late}_53g.yaml`, `simulations/config/museum_*` | museum panels | yes |
| `results/museum_qc/` | museum distance heatmap and tables | no |
| trained models, training/test data | `trained_model.Rdata`, `*_data.txt` | no (reproducible from configs) |
