#' Initially load data into Seurat
#'
#' This generates new folder (internal function only)
#' @param path_10x_data path to 10x folder
#' @param sample_names list of sample names (used when downsampling is FALSE) (default is FALSE)
#' @param lower_gene_filter_num minimum number of genes that has to be detected in a cell 
#' @param upper_gene_filter_num maximum number of genes that can be detected in a cell
#' @param percent_mito allowable percenatage of mitochondria DNA allowed
#' @param downsampled if downsampled feature was used when cell ranger was run
#' @param top_variable_genes how many genes to label in variablity plots (defualt is 10)
#' @param run_umap run umap for feature reduction
#' @param run_tsne run tsne for feature reduction
#' @param save_object save seurat normalized object in RDS file so this step doesn't have to be redone everytime (defualt is TRUE)
#' @param results_folder folder that results appear in (default: s1_quality_control_results)
#' @import Seurat
#' @import stringr
#' @import rlist
#' @import ggplot2
#' @keywords load data
#' @export

s1_load_data <- function(path_10_data, sample_names=FALSE,
                         lower_gene_filter_num=200,
                         upper_gene_filter_num=2500,
                         percent_mito=10, downsampled=TRUE,
                         top_variable_genes=10,
                         run_umap=TRUE, run_tsne=FALSE, save_object=TRUE,
                         results_folder="s1_quality_control_results") {
    results_path <- generate_folder(results_folder)
    unlink(paste0(results_folder, "/*"))

    if (isTRUE(downsampled)) {
        print("STATUS: Read 10x downsampled object")
        pbmc_data <- Read10X(data.dir = path_10_data)

        path_to_library <- str_remove(
            path_10_data,
            "filtered_feature_bc_matrix/"
        )

        print("STATUS: Load library information")
        libraries <- read.csv(paste0(path_to_library, "aggregation.csv"))
        ids <- rownames(libraries)
        samplenames <- as.character(libraries$library_id)
        print("STATUS: Create Seurat object")
        combined <- CreateSeuratObject(
            counts = pbmc_data, names.field = 2,
            names.delim = "-",
            min.cells = 1, min.features = lower_gene_filter_num
        )
        print("STATUS: Assign cell identities")
        combined$sample <- plyr::mapvalues(
            x = combined$orig.ident,
            from = ids,
            to = samplenames
        )
        combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
        VlnPlot(combined,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3, pt.size = 0.8, split.by = "sample"
        )
        ggsave(paste0(results_path, "/VnPlotMt.png"),
            width = 8, height = 6, dpi = 300
        )
        print("STATUS: Filtering out cells")
        combined <- subset(combined,
            subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10
        )

        VlnPlot(combined,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            pt.size = 0.8, split.by = "sample",
            ncol = 3
        )
        ggsave(paste0(results_path, "/VnPlotMt_Filtered.png"),
            width = 8, height = 6, dpi = 300
        )

        ##LOGNORMALIZE###
        combined <- NormalizeData(combined,
            normalization.method = "LogNormalize",
            scale.factor = 10000
        )

        ### GET MOST VARIABLE GENES###
        combined <- FindVariableFeatures(combined,
            selection.method = "vst",
            nfeatures = 2000
        )
        top <- head(VariableFeatures(combined), top_variable_genes)
        plot1 <- VariableFeaturePlot(combined)
        LabelPoints(plot = plot1, points = top, repel = TRUE)
        ggsave(paste0(
            results_path, "/VariableGenes.png"),
            dpi = 300
        )
    } else {
        counter <- 1
        pbmc_list <- list()
        for (path in path_10_data) {
            print(paste0("STATUS: Loading ", path, " ", sample_names[counter]))
            pbmc_data <- Read10X(data.dir = path)
            pbmc <- CreateSeuratObject(counts = pbmc_data,
                project = sample_names[counter],
                min.cells = 1, names.delim = "-",
                min.features = lower_gene_filter_num
            )
            pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
            VlnPlot(pbmc,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.8
            )
            ggsave(paste0(
                results_path, "/",
                sample_names[counter], "_VnPlotMt.png"
            ),
            width = 8, height = 6, dpi = 300
            )
            print("STATUS: Filtering out cells")
            pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10
            )

            VlnPlot(pbmc,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                pt.size = 0.8,
                ncol = 3
            )

            ggsave(paste0(
                results_path, "/",
                sample_names[counter], "_VnPlotMt_Filtered.png"
            ),
            width = 8, height = 6, dpi = 300
            )

            ###LOGNORMALIZE###
            pbmc <- NormalizeData(pbmc,
                normalization.method = "LogNormalize",
                scale.factor = 10000
            )

            ### GET MOST VARIABLE GENES###
            pbmc <- FindVariableFeatures(pbmc,
                selection.method = "vst",
                nfeatures = 2000
            )
            top <- head(VariableFeatures(pbmc), top_variable_genes)
            plot1 <- VariableFeaturePlot(pbmc)
            LabelPoints(plot = plot1, points = top, repel = TRUE)
            ggsave(paste0(
                results_path, "/",
                sample_names[counter], "_VariableGenes.png"
                ),
                dpi = 300
            )
            pbmc_list[[sample_names[counter]]] <- pbmc
            counter <- counter + 1
        }

        anchors <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:20)
        combined <- IntegrateData(
            anchorset = anchors, dims = 1:20
        )
        combined$sample <- combined$orig.ident
    }

    print("STATUS: scale data")
    all.genes <- rownames(combined)
    combined <- ScaleData(combined,
        vars.to.regress = "percent.mt",
        features = all.genes, verbose = TRUE
    )

    print("STATUS: feature reduction")
    combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
    if (isTRUE(run_umap)) {
        combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
    }
    if (isTRUE(run_tsne)) {
        combined <- RunTSNE(combined, reduction = "pca", dims = 1:20)
    }
    DimPlot(combined, reduction = "umap", group.by = "sample")
    ggsave(paste0(results_path, "/UMAP_sample.png"), dpi = 300)

    if (isTRUE(save_object)) {
        saveRDS(combined, paste0(results_path, "/SC_pbmc_norm_object.rds"))
    }
    return(combined)
}