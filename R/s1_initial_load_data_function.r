#' Initially load data into Seurat
#'
#' This generates new folder (internal function only)
#' @param path_10x_data path to 10x folder
#' @param lower_gene_filter_num minimum number of genes that has to be detected in a cell 
#' @param upper_gene_filter_num maximum number of genes that can be detected in a cell
#' @param percent_mito allowable percenatage of mitochondria DNA allowed
#' @param downsampled if downsampled feature was used when cell ranger was run
#' @param run_umap run umap for feature reduction
#' @param run_tsne run tsne for feature reduction
#' @param results_folder folder that results appear in (default: s1_quality_control_results)
#' @import Seurat
#' @import stringr
#' @import rlist
#' @import ggplot2
#' @keywords load data
#' @export

s1_load_data <- function(path_10_data, sample_names,
                         lower_gene_filter_num=200,
                         upper_gene_filter_num=2500,
                         percent_mito=10, downsampled=TRUE,
                         run_umap=TRUE, run_tsne=FALSE,
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
        pbmc <- CreateSeuratObject(
            counts = pbmc_data, names.field = 2,
            names.delim = "-",
            min.cells = 1, min.features = lower_gene_filter_num
        )
        print("STATUS: Assign cell identities")
        pbmc$sample <- plyr::mapvalues(
            x = pbmc$orig.ident,
            from = ids,
            to = samplenames
        )
    } else {
        counter <- 1
        pbmc_list <- list()
        for (path in path_10_data) {
            print(path)
            pbmc_data <- Read10X(data.dir = path)
            pbmc <- CreateSeuratObject(counts = pbmc_data,
                project = sampleID,
                min.cells = 1,
                min.features = lower_gene_filter_num
            )
            list.append(pbmc_list, sample_names[counter] <- pbmc)
            counter <- counter + 1
        }
        anchors <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:20)
        pbmc <- IntegrateData(anchorset = anchors, dims = 1:20)
        pbmc$sample <- pbmc$orig.ident
    }

    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

    VlnPlot(pbmc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, group.by = "sample", pt.size = 0.8
    )
    ggsave(paste0(results_path, "VnPlotMt.png"),
        width = 8, height = 6, dpi = 300
    )

    print("STATUS: Filtering out cells")
    pbmc <- subset(pbmc,
       subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10
    )

    VlnPlot(pbmc,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample", pt.size = 0.8,
        ncol = 3
    )
    ggsave(paste0(results_path, "VnPlotMt_Filtered.png"),
        width = 8, height = 6, dpi = 300
    )

    print("STATUS: Normalize data")
    all.genes <- rownames(pbmc)

    pbmc <- NormalizeData(pbmc,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )

    pbmc <- FindVariableFeatures(pbmc,
        selection.method = "vst",
        nfeatures = 2000
    )

    pbmc <- ScaleData(pbmc,
        vars.to.regress = "percent.mt",
        features = all.genes, verbose = FALSE
    )

    print("STATUS: feature reduction")
    pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
    if (isTRUE(run_umap)) {
        pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
    }
    if (isTRUE(run_tsne)) {
        pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:20)
    }
    DimPlot(pbmc, reduction = "umap", group.by = "sample")
    ggsave(paste0(results_path, "UMAP_sample.png"), dpi = 300)

    return(pbmc)
}