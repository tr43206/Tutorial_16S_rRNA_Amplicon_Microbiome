# Microbiome_Tutorial_16S-rRNA-amplicon

Output fastq files of Gut-Lung Axis project from MGI DNBSEQ-G99.



# Step 0 : Preprocessing

1) Demultiplex fastq files by generating `manifest file`.
2) Low freauency samples `trimming`, `denoising`, `merging` forward and reverse reads, and `chimera removal` using DADA2.
3) Visualize output files of DADA2 and choose which `depth threshold` is appropriate.



# Step 1 : QIIME 2 (v qiime2-amplicon-2024.10)

1) Generate a tree for phylogenetic diversity analyses (e.g., `Faith PD`, `UniFrac`)
2) Alpha and beta diversity analysis by generating `core-metrics`.
3) Visualize (Alpha Div. - `qza files` / Beta Div. - `distance qza files`, export `pcoa files` to folder, and check pc proportion with `emperor qzv` files).

4) ASVs classification using `classify-learn`.
5) Collapse to genus level (or family, phylum, whatever).
6) Convert raw read frequency to relative frequency.
7) Export frequency tables (raw and relative, both) to folder, and convert inner `biom files` to `tsv files`.



# Step 2 : Visualization

1) Alpha diversity metrices to box plots - Statistics: Kruskal-Wallis
2) Beta diversity metrices to PCoA plots - Statistics: PERMANOVA
3) Relative frequency tables to bar plots - Statistics: Differential abundance (DA) methods (e.g., LEfSe, ALDEs2, ANCOM-BC2, MaAsLin3)



# Step 3 : Key features identification

1) More strict abundance and prevalence filtering (to minimize overfitting in machine learning-based models caused by rare taxa)
2) Machine learning-based models (e.g., Random forest, Support vector machine (SVM)), or Linear regression-based key features identification.
