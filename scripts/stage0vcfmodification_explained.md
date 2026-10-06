### Stage 0: why the master VCF is modified

`RedWolfProject_53genomes_geno10_maf1.vcf.gz` was exported by PLINK 1.9, so it
doesn't match what the rest of the pipeline expects:

1. **Sample IDs** are PLINK `FID_IID` (`ALcoy_ALcoy`). They are stripped to the
   IID (`ALcoy`); the mapping is written to `*.sample_map.tsv`.
2. **Chromosomes** use PLINK dog codes (`1`–`38`, `39` = X, `40`/`41` =
   MT/unplaced). They are renamed to canFam3.1 names (`chr01`…`chr38`, `chrX`)
   to match the Beagle maps and the reference FASTA (`config/chr_rename.tsv`).
   `40`/`41` are dropped.
3. **REF/ALT are PLINK's major/minor alleles, not the reference allele**
   (`INFO/PR`, "provisional reference"). About 18% of sites had the canFam3.1
   allele in ALT. The genotypes were consistent with that labelling, but ArchIE
   reads allele `1` as derived. That makes the IFS, the reference filter
   behind S\*, and the private-SNP count depend on the coding; distance-based
   features don't change. Sites are swapped back to canFam3.1 with
   `bcftools +fixref -m flip-all`, which exchanges REF/ALT and recodes GTs.
   The script aborts if any strand complement is applied or any mismatch
   remains, and writes counts to `*.prep_report.txt`. In the 53-genome run,
   2,556,745 of 13,833,055 sites were swapped and 0 strand flips were needed.
   With this coding `0` = canFam3.1 (dog) allele, a proxy for the ancestral
   state that the earlier 63-genome (GATK) analysis also used.
