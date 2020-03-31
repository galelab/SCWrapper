
# Gale Lab single-cell analysis pipeline

Goal of this software is to provide a pipeline for lab members to analyze gene expression data


## Install 
First load library('devtools'). If you do not have devtools run install.packages('devtools). To download in R please use install_github('galelab/SCWrapper', dependencies=TRUE)

## How to execute pipeline
SCWrapper currently has 3 primary functions
1. s1_load_data - function will load and normalize single cell gene expression data using Seurat.  Data can be loaded in two ways: 1) import downsampled data from CellRanger agg function or 2) load individual analyses from CellRanger count.  Seurat parameters for loading data are all defualt except for min.cells which is set equal to 1 instead of 3.  All results will be outputed into the defualt results folder s1_quality_control_results/ in your working directory. Examples below:

    downsample example = s1_load_data("path/to/Fresh_Aggregate/outs/filtered_feature_bc_matrix", downsampled=TRUE)

    Individual sample example = s1_load_data(c("path/to/sample1/outs/filtered_feature_bc_matrix/", "path/to/sample2/outs/filtered_feature_bc_matrix"), sample_names=c("sample1", "sample2"), downsampled=FALSE)

2. classify_pbmc_cells - function will perform cell-type classification with garnett.  Will need a marker file and user to use a pre-compiled classifier or generate their own from their data.  Example marker files and compiled classifiers can be found at https://cole-trapnell-lab.github.io/garnett/classifiers/. Examples below:

    downsampled example = 

5. s4_gene_enrichment_analysis - function will perform over representation analysis (ORA) on genes significantly expressed identified in step s3_DE_anlysis.  Step 3 clusters these genes according to expression patterns (see heatmap generated in step s3_DE_anlysis) into modules.  ORA is done on each module.  If user wants to look at ORA of a specific comparison specify the comparison in the comparison option of the s4_gene_enrichment_analysis function.  Also the user can input their own rank file and perform ORA on that instead of gene lists in step s3_DE_anlysis (see examples below). For each comparison an output folder will be generated with a number of output files such as:

    over_enrich_barplotall_ge.png - bar plot of significantly enriched GO terms using all differentially expressed genes

    over_enrich_barplotup_ge.png - bar plot of significantly enriched GO terms using up regulated differentially expressed genes

    over_enrich_barplotdown_ge.png - bar plot of significantly enriched GO terms using down regulated differentially expressed genes

    over_enrich_cneplotall_ge.png - cnet plot (network) of significantly enriched GO terms using all differentially expressed genes

    over_enrich_cneplotup_ge.png - cnet plot (network) of significantly enriched GO terms using up regulated differentially expressed genes

    over_enrich_cneplotdown_ge.png - cnet plot (network) of significantly enriched GO terms using down regulated differentially expressed genes

    *all_ora_genes.csv tables of genes and corresponding log fold changes in GO terms when all differentially expresed genes were used for analysis

    *up_ora_genes.csv tables of genes and corresponding log fold changes in GO terms when up regulated differentially expresed genes were used for analysis

    *down_ora_genes.csv tables of genes and corresponding log fold changes in GO terms when down regulated differentially expresed genes were used for analysis

## Example of pipeline run
Commands for running pipeline

1. p1_modify_count_matrix(countfile='count_matrix.txt', targetfile='target.csv', samples_to_remove_count_matrix='samplestoremove.csv')

2. s1_normalize_raw_counts(countfile='./p1_modified_count_matrix_results/count_matrix_mod.txt', targetfile='./p1_modified_count_matrix_results/targets_mod.csv', gene_conversion_file='rhesus2human.csv', target_column=3, batch_column=2,filter_genes_below_counts=50)

3. s2_feature_reduction(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/targets_mod.csv', target_columns=c(2,5))

4. s3_DE_analysis(countfile='./s1_norm_raw_counts_results/1.norm_matrix.txt', targetfile='./p1_modified_count_matrix_results/target_file.csv', gene_conversion_file='rhesus2human.csv',  'matrixfile=MATRIXEXAMPLE.txt', pvalue=0.05, logfoldchange=1.5)

* This shows couple examples of how s4_gene_enrichment_analysis could be run
5. s4_gene_enrichment_analysis(go_enrich_type='BP', universe=TRUE, modules=TRUE, NumTopGoTerms=30) - This will perform ORA on genes classified in specific modules from step s3_DE_analysis ONLY as well as use all expressed genes as background (collected from ./s1_norm_raw_counts_results/1.norm_matrix_HGNC.txt) instead of whole genome (because universe is specified as true) and create tables and figures including the top 30 enriched go terms.

5. s4_gene_enrichment_analysis(go_enrich_type='BP', comparison='treatmentInfected_w19-treatmentUninfected_w19' universe=TRUE, modules=TRUE, NumTopGoTerms=30) - This will perform ORA on genes classified in specific modules from step s3_DE_analysis and genes significantly differentially expressed between the comparison 'treatmentInfected_w19-treatmentUninfected_w19' as well as use all expressed genes as background (collected from ./s1_norm_raw_counts_results/1.norm_matrix_HGNC.txt) instead of whole genome (because universe is specified as true) and create tables and figures including the top 30 enriched go terms.

5. s4_gene_enrichment_analysis(go_enrich_type='BP', comparison='treatmentInfected_w19-treatmentUninfected_w19' universe=FALSE, modules=FALSE, NumTopGoTerms=30) - This will perform ORA on genes classified in specific modules from step s3_DE_analysis and genes significantly differentially expressed between the comparison 'treatmentInfected_w19-treatmentUninfected_w19' as well as use whole genome (because universe is specified as false) as background genes and create tables and figures including the top 30 enriched go terms.

## Other important notes 
1. matrix file must be layed out in the exact same manner in as the example file (MATRIXEXAMPLE.txt).  The word 'treatment' must be in front of each comparison as in the example file.
2. different conversion file can be used but it must be in csv format layed out in the same manner as rhesus2human.csv
3. target file must be layed out in a similar manner as the example file (targetexample.csv) but amount of information in the target file will vary with each experiment