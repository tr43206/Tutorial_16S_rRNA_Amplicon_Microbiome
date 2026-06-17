# Microbiome_Tutorial_16S-rRNA-amplicon

Output fastq files of Gut-Lung Axis project from MGI DNBSEQ-G99.



# Step 0 : Preprocessing

1) Demultiplex fastq files by generating `manifest file`. Use `01. manifest_auto.ipynb`.
2) Low frequency samples `trimming`, `denoising`, `merging` forward and reverse reads, and `chimera removal` using DADA2. Normally, `Q20` or `Q30` is used for quality threshold.
3) Visualize output files of DADA2 and choose which `depth threshold` is appropriate.



# Step 1 : QIIME 2 (v qiime2-amplicon-2024.10)

1) Generate a tree for phylogenetic diversity analyses (e.g., `Faith PD`, `UniFrac`).
   ```bash
   qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences merged_rep_seq.qza \
--output-dir mafft-fasttree \
--p-n-threads 10
   ```
3) Alpha and beta diversity analysis by generating `core-metrics`.
4) Visualize (Alpha Div. - `qza files` / Beta Div. - `distance qza files`, export `pcoa files` to folder, and check pc proportion with `emperor qzv` files).

5) ASVs classification using `classify-learn`.
6) Collapse to genus level (or family, phylum, whatever).
7) Convert raw read frequency to relative frequency.
8) Export frequency tables (raw and relative, both) to folder, and convert inner `biom files` to `tsv files`.

9) Visualization



# Step 2 : Key features identification

1) More strict abundance and prevalence filtering (to minimize overfitting in machine learning-based models caused by rare taxa)
2) Machine learning-based models (e.g., Random forest, Support vector machine (SVM)), or Linear regression-based key features identification.
