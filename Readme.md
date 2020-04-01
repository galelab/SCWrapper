
# Gale Lab single-cell analysis pipeline

Goal of this software is to provide a pipeline for lab members to analyze gene expression data


## Install 
First load library('devtools'). If you do not have devtools run install.packages('devtools). To download in R please use install_github('galelab/SCWrapper', dependencies=TRUE)

## How to execute pipeline
SCWrapper currently has 3 primary functions
1. s1_load_data - function will load and normalize single cell gene expression data using Seurat.  Data can be loaded in two ways: 1) import downsampled data from CellRanger agg function or 2) load individual analyses from CellRanger count.  Seurat parameters for loading data are all defualt except for min.cells which is set equal to 1 instead of 3 and names.delim which is set to "-".  All results will be outputed into the defualt results folder s1_quality_control_results/ in your working directory. Examples below

- downsample example: pbmc <- s1_load_data("/path/to/Fresh_Aggregate/outs/filtered_feature_bc_matrix", downsampled=TRUE)
- Individual sample example: pbmc <- s1_load_data(c("path/to/sample1/outs/filtered_feature_bc_matrix/", "path/to/sample2/ outs/filtered_feature_bc_matrix"), sample_names=c("sample1", "sample2"), downsampled=FALSE)

2. classify_pbmc_cells - function will perform cell-type classification with garnett.  Will need a marker file and user to use a pre-compiled classifier or generate their own from their data. All results will be outputed into the defualt results folder s2_garnett_classifier_results/ in your working directory.  Example marker files and compiled classifiers can be found at https://cole-trapnell-lab.github.io/garnett/classifiers/. Note that parameters that are altered from the defualt parameters in garnett are:

- In train_cell_classifier function cluster_extend_max_frac_incorrect=0, cluster_extend_max_frac_unknown=1, and verbose=TRUE 

- Examples below:

    downsampled example: pbmc <- s2_classify_cells(pbmc, "/path/to/Fresh_Aggregate/" marker_file_path="/path/to/marker/markerfile.txt", classifier="/path/to/classifier/classifier.rds")

    not downsampled example: pbmc <- s2_classify_cells(pbmc, c("/path/to/sample1/", "/path/to/sample2"), sample_names=c("sample1", "sample2"), marker_file_path="/path/to/marker/markerfile.txt", classifier="/path/to/classifier/classifier.rds")  - samples have to be in the same order as step 1.  


3. s3_DE_analysis - function will perform over differential gene expression.  Will perform DE analysis on individual clusters between conditions of option is specified. Genes are considered DE if log fold change is >= 1.2 and adjusted P value < 0.05 (both of these are adjustable with parameters LFC and pvalue). All results will be outputed into the defualt results folder s3_DE_results/ in your working directory. Examples below (same no matter if downsampling was used or not): 

- run DE between to conditions: s3_DE_analysis(pbmc, "condition1", "condition2")
- run DE between clusters in two different conditions: s3_DE_analysis(pbmc, "condition1", "condition2", cluster_DE=TRUE)
- run DE between conditions and generate violin plots of DE genes: s3_DE_analysis(pbmc, "condition1", "condition2", individual_plots=TRUE)


## Other important notes 
1.  <a href="https://satijalab.org/seurat/" target="_blank">  The single cell R library Seuratis used in this Wrapper.  
2.  <a href="https://cole-trapnell-lab.github.io/garnett/docs/" target="_blank"> The cell type classifier garnett is used in this Wrapper.

## References 
* Integrating single-cell transcriptomic data across different conditions, technologies, and species Butler et al. 2018 <a href="https://www.nature.com/articles/nbt.4096" target="_blank">
* Supervised classification enables rapid annotation of cell atlases Pilner et al. 2019 <a href="https://www.nature.com/articles/s41592-019-0535-3" target="_blank">