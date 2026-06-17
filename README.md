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
   --p-n-threads 32
   ```
2) Alpha and beta diversity analysis by generating `core-metrics`.
   ```bash
   qiime diversity core-metrics-phylogenetic \
   --i-phylogeny mafft-fasttree/rooted-tree.qza \
   --i-table filtered_table.qza \
   --p-sampling-depth 1816 \
   --m-metadata-file Aim1_metadata.txt \
   --output-dir core-metrics-results-1816
   ```
3) Visualize (Alpha Div. - `qza files` / Beta Div. - `distance qza files`, export `pcoa files` to folder, and check pc proportion with `emperor qzv` files).
   ```bash
   qiime diversity alpha-group-significance \
   --i-alpha-diversity core-metrics-phylogenetic/faith_pd_vector.qza \
   --m-metadata-file Aim1_metadata.txt \
   --o-visualization core-metrics-phylogenetic/faith_pd_vector_group.qzv
   ```
   ```bash
   qiime diversity beta-group-significance \
   --i-distance-matrix core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
   --m-metadata-file Aim1_metadata.txt \
   --m-metadata-column Sample_name \
   --o-visualization core-metrics-phylogenetic/unweighted_unifrac_distance_group_significance.qzv \
   --p-pairwise
   ```
4) ASVs classification using `classify-learn`.
   ```bash
   qiime feature-classifier classify-sklearn \
   --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
   --i-reads merged_rep-seqs.qza \
   --o-classification gg2_taxonomy.qza
   ```
5) Collapse to genus level (or family, phylum, whatever).
    ```bash
    qiime taxa collapse \
    --i-table merged_table.qza \
    --i-taxonomy gg2_taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table gg2_genus_table.qza
    ```
6) Convert raw read frequency to relative frequency.
    ```bash
    qiime feature-table relative-frequency \
    --i-table gg2_genus_table.qza \
    --o-relative-frequency-table rel_gg2_genus_table.qza
    ```
7) Export frequency tables (raw and relative, both) to folder, and convert inner `biom files` to `tsv files`.
    ```bash
    qiime tools export \
    --input-path rel_gg2_genus_table.qza \
    --output-path rel_gg2_genus_table
    ```
    ```bash
    biom convert \
    -i feature-table.biom \
    -o rel_gg2_genus_table.tsv \
    --to-tsv
    ```
8) Visualization



# (Optional) Feature table filtering and combining

1) Check file format
   ```bash
   qiime tools peek ancombc-5-antrum_case1-case2.qza
   ```
2) Feature table subgrouping
   ```bash
   qiime feature-table filter-samples \
   --i-table merged_table.qza \
   --m-metadata-file Cheese_mapping.txt \
   --p-where "[Sample_name] IN ('C.NT_2','c_fos','c_inulin','c_kestose','c_sucrose')" \
   --output-dir filter_merge \
   --p-no-exclude-ids
   ```
3) Taxonomy-based filtering of tables and sequences
   ```bash
   qiime taxa filter-table \
   --i-table table.qza \
   --i-taxonomy taxonomy.qza \
   --p-include p__ \
   --o-filtered-table table-with-phyla.qza
   ```
   ```bash
   qiime taxa filter-table \
   --i-table table.qza \
   --i-taxonomy taxonomy.qza \
   --p-include p__ \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table table-with-phyla-no-mitochondria-no-chloroplast.qza
   ```bash
   qiime taxa filter-table \
   --i-table table.qza \
   --i-taxonomy taxonomy.qza \
   --p-mode exact \
   --p-exclude "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__mitochondria" \
   --o-filtered-table table-no-mitochondria-exact.qza
   ```
4) Filtering sequences
   ```bash
   qiime taxa filter-seqs \
   --i-sequences sequences.qza \
   --i-taxonomy taxonomy.qza \
   --p-include p__ \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-sequences sequences-with-phyla-no-mitochondria-no-chloroplast.qza
   ```
5) Combine collections of feature tables
   ```bash
   qiime feature-table merge-seqs \
   --i-data 17_Trim_6_6_235_232/representative_sequences.qza \
   --i-data 18_Trim_6_6_242_230/representative_sequences.qza \
   --o-merged-data merged_representative_sequences.qza
   ```




# (Optional) Utilizing NCBI data

1) NCBI data download
   ```bash
   conda create -n sra
   ```
   ```bash
   mamba install -c bioconda entrez-direct sra-tools pigz
   ```
   ```bash
   esearch -db sra -query PRJEB27662 \
   | efetch -format runinfo \
   | awk -F',' 'NR>1{print $1}' \
   | grep -E '^(ERR|SRR|DRR)' \
   | xargs -n 1 -P 4 fasterq-dump -e 16 -p -t ./tmp -O fastq
   ```
2) Importing with NCBI downloaded data
   ```bash
   qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path manifest.txt \
   --output-path paired-end-demux.qza \
   --input-format PairedEndFastqManifestPhred33V2
   ```



# (Optional) fastqc

1) fastqc install
   ```bash
   conda env create -n fastqc -y
   ```
   ```bash
   conda activate fastqc
   ```
   ```bash
   cd ~/Desktop/TaehunRoh
   ```
   ```bash
   wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
   ```
   ```bash
   unzip fastqc_v0.12.1.zip
   ```
   ```bash
   cd FastQC
   ```
   ```bash
   sudo apt install default-jre
   ```
   ```bash
   chmod +x fastqc
   ```
2) Run fastqc
   ```bash
   cat *R1*.fastq.gz > Run45_forward_R1.fastq.gz
   ```
   ```bash
   cat *R2*.fastq.gz > Run45_reverse_R2.fastq.gz
   ```
   ```bash
   ./fastqc ../00.Raw_data/Dry2014_1.fastq.gz ../00.Raw_data/Dry2014_2.fastq.gz ../00.Raw_data/Wet2014_1.fastq.gz ../00.Raw_data/Wet2014_2.fastq.gz
   ```



# (Optional) Rarefaction

1) Rarefaction
   ```bash
   qiime diversity alpha-rarefaction \
   --i-table table-week6-UVBvs38A1High.qza \
   --i-phylogeny mafft-fasttree/rooted_tree.qza \
   --p-max-depth 91316 \
   --p-metrics shannon \
   --m-metadata-file ~/Desktop/SolheeKim/Run43/metadata_week6_UVB.txt \
   --o-visualization rarefaction_shannon
   ```



# (Optional) Alpha and beta diversity analysis (not in core-metrics)

1) Alpha
   ```bash
   qiime diversity alpha-phylogenetic \
   --i-table filtered_table.qza \
   --i-phylogeny rooted-tree.qza \
   --p-metric faith_pd \
   --o-alpha-diversity faith_pd_vector.qza \
   --verbose
   ```
   ```bash
   qiime diversity alpha \
   --i-table filtered_table.qza \
   --p-metric simpson \
   --o-alpha-diversity simpson_vectors.qza \
   --verbose
   ```  
2) Beta
   ```bash
   qiime diversity beta-phylogenetic \
   --i-table filtered_table.qza \
   --i-phylogeny rooted-tree.qza \
   --p-metric unweighted_unifrac \
   --o-distance-matrix unweighted_unifrac_distance_matrix.qza \
   --verbose
   ```
   ```bash
   qiime diversity beta \
   --i-table filtered_table.qza \
   --p-metric braycurtis \
   --o-distance-matrix unweighted_unifrac_distance_matrix.qza \
   --verbose
   ```



# (Optional) Filtering distance matrices

1) Beta
   ```bash
   qiime diversity filter-distance-matrix \
   --i-distance-matrix distance-matrix.qza \
   --m-metadata-file samples-to-keep.tsv \
   --o-filtered-distance-matrix identifier-filtered-distance-matrix.qza
   ```
   ```bash
   qiime diversity filter-distance-matrix \
   --i-distance-matrix distance-matrix.qza \
   --m-metadata-file sample-metadata.tsv \
   --p-where "[subject]='subject-2'" \
   --o-filtered-distance-matrix subject-2-filtered-distance-matrix.qza
   ```



# (Optional) Download outdated Greengenes database

1) Download
   ```bash
   wget https://ftp.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz
   ```
   ```bash
   tar -xzf gg_13_5_otus.tar.gz
   ```
2) Import to QIIME 2
   ```bash
   qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path gg_13_5_otus/rep_set/99_otus.fasta \
   --output-path gg13_5-99-seqs.qza
   ```
   ```bash
   qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --input-path gg_13_5_otus/taxonomy/99_otu_taxonomy.txt \
   --input-format HeaderlessTSVTaxonomyFormat \
   --output-path gg13_5-99-tax.qza
   ```
3) Classifier training after specific target region extraction (e.g., V4 region - 515F/806R)
   ```bash
   qiime feature-classifier extract-reads \
   --i-sequences gg13_5-99-seqs.qza \
   --p-f-primer GTGCCAGCMGCCGCGGTAA \
   --p-r-primer GGACTACHVGGGTWTCTAAT \
   --o-reads gg13_5-99-515F806R.qza
   ```
   ```bash
   qiime feature-classifier fit-classifier-naive-bayes \
   --i-reference-reads gg13_5-99-515F806R.qza \
   --i-reference-taxonomy gg13_5-99-tax.qza \
   --o-classifier gg13_5-99-515F806R-classifier.qza
   ```
4) Application
   ```bash
   qiime feature-classifier classify-sklearn \
   --i-classifier gg13_5-99-515F806R-classifier.qza \
   --i-reads rep-seqs.qza \
   --o-classification taxonomy.qza
   ```



# (Optional) txt to qza

1) txt to qza
   ```bash
   qiime tools import \
   --input-path old_distance_matrix.txt \
   --output-path new_distance_matrix.qza \
   --type DistanceMatrix
   ```



# Step 2 : Pathways prediction

1) PICRUSt2
   ```bash
   picrust2_pipeline.py \
   -s sequences.fasta \
   -i otu_for_picrust.biom \
   -o picrust-MPGA/
   ```
2) MICOM
   ```bash
   qiime micom build \
   --i-abundance filtered_table_min-freq-1000.qza \
   --i-taxonomy gtdb207_taxonomy.qza \
   --i-models ../agora201_gtdb207_genus_1.qza \
   --p-cutoff 0 \
   --p-threads 4 \
   --o-community-models models.qza \
   --verbose
   ```
   ```bash
   qiime micom tradeoff \
   --i-models models.qza \
   --i-medium ../western_diet_gut_agora.qza \
   --p-threads 4 \
   --o-results tradeoff.qza \
   --verbose
   ```
   ```bash
   qiime micom plot-tradeoff \
   --i-results tradeoff.qza \
   --o-visualization tradeoff.qzv
   ```
   ```bash
   qiime micom grow \
   --i-models models.qza \
   --i-medium ../western_diet_gut_agora.qza \
   --p-tradeoff 1.0 \
   --p-threads 32 \
   --o-results growth.qza \
   --verbose
   ```
   ```bash
   qiime micom plot-growth \
   --i-results growth.qza \
   --o-visualization growth.qzv
   ```
   ```bash
   qiime micom exchanges-per-sample \
   --i-results growth.qza \
   --o-visualization exchanges_gtdb207.qzv
   ```
   ```bash
   qiime micom exchanges-per-taxon \
   --i-results growth.qza \
   --o-visualization niche.qzv
   ```
   ```bash
   qiime micom association \
   --i-results growth.qza \
   --p-fdr-threshold 0.05 \
   --m-metadata-file ../RUN45_Metadata_GutLungAxis.txt \
   --m-metadata-column Group \
   --o-visualization associations.qzv
   ```
3) metnet
   ```bash
   qiime metnet generateFeatures \
   --i-frequency filtered_table_min-freq-1000_transposed.qza \
   --i-taxa silva_taxonomy.qza \
   --p-selection AGORAv103 \
   --p-level g \
   --o-reactions rxns_scores.qza \
   --o-subsystems subs_scores.qza \
   --o-xmatrix Xmatrix.qza
   ```
   ```bash
   qiime metnet plotPCA \
   --i-table rxns_scores.qza \
   --m-sample-metadata-file RUN45_Metadata_GutLungAxis_min-freq-1000_Week0.txt \
   --m-sample-metadata-column Group \
   --o-visualization rxns_pca.qzv
   ```
   ```bash
   qiime metnet plotClusteMap \
   --i-table rxns_scores.qza \
   --m-sample-metadata-file RUN45_Metadata_GutLungAxis_min-freq-1000_Week0.txt \
   --m-sample-metadata-column Group \
   --o-visualization rxns_hclust.qzv
   ```
   ```bash
   qiime metnet differentialReactions \
   --i-reactions rxns_scores.qza \
   --m-metadata-file RUN45_Metadata_GutLungAxis_min-freq-1000_Week0.txt \
   --m-metadata-column Group \
   --p-condition-name VNAM \
   --p-control-name Control \
   --p-selection-model AGORAv103 \
   --o-differential-analysis diff_reactions.qza
   ```
   ```bash
   qiime metnet differentialSubSystems \
   --i-subsystems subs_scores.qza \
   --m-metadata-file RUN45_Metadata_GutLungAxis_min-freq-1000_Week0.txt \
   --m-metadata-column Group \
   --p-condition-name VNAM \
   --p-control-name Control \
   --o-differential-analysis diff_subs.qza
   ```
   ```bash
   qiime metadata tabulate \
   --m-input-file diff_subs.qza \
   --o-visualization diff_subs.qzv
   ```
   ```bash
   qiime metnet differentialExchanges \
   --i-reactions rxns_scores.qza \
   --m-metadata-file RUN45_Metadata_GutLungAxis_min-freq-1000_Week0.txt \
   --m-metadata-column Group \
   --p-condition-name VNAM \
   --p-control-name Control \
   --p-selection-model AGORAv103 \
   --o-differential-analysis diff_ex.qza
   ```
   ```bash
   qiime metadata tabulate \
   --m-input-file diff_ex.qza \
   --o-visualization diff_ex.qzv
   ```
   ```bash
   qiime metnet plotBoxplot \
   --i-table subs_scores.qza \
   --i-differentialresults diff_subs.qza \
   --m-sample-metadata-file RUN45_Metadata_GutLungAxis_min-freq-1000_Week0.txt \
   --m-sample-metadata-column Group \
   --p-condition-name VNAM \
   --p-control-name Control \
   --p-namefeature "S65 | Tryptophan metabolism" \
   --o-visualization sub_boxplot.qzv
   ```



# Step 3 : Key features identification

1) More strict abundance and prevalence filtering (to minimize overfitting in machine learning-based models caused by rare taxa)
2) Machine learning-based models (e.g., Random forest, Support vector machine (SVM)), or Linear regression-based key features identification.



# (Optional) Nested-cross validation (NCV)-based machine learning on QIIME 2

1) classify samples (NCV)
   ```bash
   qiime sample-classifier classify-samples-ncv \
   --i-table Taxa/all_collapse.qza \
   --m-metadata-file qiime2_metadata.txt \
   --m-metadata-column <예측하고자 하는 column명> \
   --p-estimator RandomForestClassifier \
   --p-n-estimators 2000 \
   --p-random-state 42 \
   --output-dir q2-sample-classifier/classify-samples-ncv \
   --p-parameter-tuning
   ```
2) Confution-matrix
   ```bash
   qiime sample-classifier confusion-matrix \
   --i-predictions /classify-samples-ncv/predictions.qza \
   --i-probabilities classify-samples-ncv/probabilities.qza \
   --m-truth-file qiime2_metadata.txt \
   --m-truth-column <예측하고자 하는 column명> \
   --o-visualization classify-samples-ncv/ncv-confusion-matrix.qzv
   ```
3) Feature importance scores
   ```bash
   qiime metadata tabulate \
   --m-input-file classify-samples-ncv/ feature_importance.qza \
   --o-visualization classify-samples-ncv/ feature_importance.qzv
   ```
