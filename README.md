# Phlag analysis of the avian phylogeny
This repository contains the Phlag analysis results, gene trees, species trees, and supporting data for the avian dataset (experiment E3) presented in:

> Şapcı AOB, Arasti S, Braun EL, Mirarab S. **Phlag: Scalable detection of genomic regions with unexplained phylogenetic heterogeneity.** Bioinformatics, ISMB issue (2026).

Phlag is available at [github.com/bo1929/phlag](https://github.com/bo1929/phlag).

## Files

- `sorted_genetrees/`: Compressed gene trees for each chromosome, sorted by genomic position.
  * `gene_trees-Stiller2024-chr{CHR}-sorted.nwk.gz`: Gene trees from the 363-taxon avian dataset by [Stiller et al.](https://doi.org/10.1038/s41586-024-07323-1), selecting 1Kbp subalignments with minimum missing data from each 10Kbp segment. Gene trees were estimated using IQ-TREE under the GTR+G4 model with approximate Bayesian support.
  * Covers 28 autosomes (chr1–chr28) and the Z chromosome (chrZ).

- `qqs/`: Precomputed quartet quartet site (QQS) frequencies, used as input to Phlag.
  * `gene_trees-Stiller2024-chr{CHR}-sorted-qqs.txt.gz`: Per-chromosome QQS frequencies.
  * `qqs-chr1-5.txt.gz`: Concatenated QQS for the five macrochromosomes (chr1–chr5).
  * `qqs-chr1-5-Z.txt.gz`: Concatenated QQS for the five macrochromosomes plus the Z chromosome.

- `coordinates/`: Coordinate files for each chromosome assembly.
  * `coordinates-Stiller2024-chr{CHR}.txt`: Genomic coordinates for gene tree loci.

- `positions-gene_trees-Stiller2024-chr{CHR}-sorted.txt`: Genomic positions (chromosome and base pair coordinate) for each gene tree window.

- Species trees:
  * `main.tre`: Main avian species tree with CU branch lengths (363 taxa; Stiller et al.).
  * `main_alternative_mod.tre`: Modified species tree.
  * `main-num_generations.tre`: Species tree with branch lengths in number of generations.
  * `63K_dated.tre` / `2023-04-dated.tre`: Dated species trees with branch lengths in time (million years).
  * `castlespro_stiller.rooted.tre`: Species tree with branch lengths in substitution units, estimated by CASTLES-pro.

- Flagged regions (Z chromosome analysis):
  * `flagged_regions-Stiller2024-chrZ-species_tree.nwk`: Species tree used for flagged Z chromosome analysis.
  * `flagged_gene_trees-Stiller2024-chrZ-sorted.nwk`: Gene trees from flagged regions of the Z chromosome.
  * `main_species_tree-Stiller2024-chrZ_flagged.nwk.support`: Species tree with support values from flagged Z chromosome gene trees.

- `mapping_stiller_fig2a.tsv` / `mapping_stiller_main.tsv` / `mapping_stiller_merged.tsv`: Mapping from internal node labels (`N{ID}`) to branch numbers used in Stiller et al. Fig. 2a and in the paper figures.

- `list_chr.txt`: List of all 34 chromosome/scaffold identifiers in the assembly.

- `plot.R`: R script (using ggplot2) for generating visualization plots from Phlag predictions.

## Phlag predictions

Phlag was run on individual chromosomes and on concatenations of gene trees across the five macrochromosomes (chr1–chr5), with and without the Z chromosome. Each prediction directory contains:
  * `distances_chr{CHR}.txt` / `distances_concat.txt`: Hellinger distance between the null and alternative emission distributions for each branch.
  * `preds-chr{CHR}-{PARAMS}.txt` / `preds-concat-{PARAMS}.txt`: Tab-separated prediction matrix with one row per gene tree and one column per branch. Values are binary (0 = null, 1 = alternative, nan = branch not present in the gene tree).

The naming convention encodes hyperparameters as follows:
  * `eap{beta}`: expected number of anomalies (`--expected-num-anomalies`), e.g., `eap40`: beta = 40.
  * `ep{1-rho}`: expected anomaly proportion (`--expected-anomaly-proportion`), e.g., `ep005`: 1-rho = 0.05.
  * `penalty{lambda}` / `npenalty{lambda}`: prior penalty strength (`lambda`), with `n` prefix indicating a negative value.
  * `noprior`: segment mode (no MSC-based prior on emissions).

Available parameter combinations (per-chromosome predictions):
  * `eap40_ep005_penalty15/` (prior-updated, the main paper analysis)
  * `eap40_ep005_noprior/` (segment mode)
  * `eap40_ep005_npenalty15/` (negative penalty)

Concatenated analyses:
  * `preds-concat-eap100_ep005_penalty15-chr_1_5.txt` and `distances_concat-chr_1_5.txt`: Five macrochromosomes concatenated, beta = 100 (the main paper analysis).
  * `preds-concat-eap100_ep005_penalty15-chr_1_5_Z.txt` and `distances_concat-chr_1_5_Z.txt`: Five macrochromosomes + Z chromosome, beta = 100.
  * `preds-merged-eap100_ep005_penalty15.txt` and `distances_concat_merged.txt`: Merged analysis.
  * `eap40_ep005_noprior-concat_1_5/` and `eap40_ep005_noprior-concat_1_5_Z/`: Concatenated analyses without prior.

## Description
Phlag was applied to 39,849 gene trees sampled from five macrochromosomes (chr1–chr5) of the 363-taxon avian dataset by Stiller et al. Gene trees were estimated by selecting 1Kbp subalignments with minimum missing data from each 10Kbp segment and running IQ-TREE under GTR+G4. We focused on 20 branches near the base of Neoaves, selected from 33 basal branches (60–68 Mya) highlighted by Stiller et al. that had quartet support below 0.5. Each branch was analyzed individually (single focal branch) using the prior-updated mode with topology-order emissions.