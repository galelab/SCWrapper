#' Classify PBMC SC data
#'
#' This generates new folder (internal function only)
#' @param seurat_data filtered data from seurat
#' @param path_10x_data path to 10x folder
#' @param marker_file_path path to the marker file needed for garnett classification
#' @param pbmc_classifier path to pre built classifier if one is being used (defualt is to build its own classifier)
#' @param cds_gene_type type of gene ID (defualt is Ensembl)
#' @param results_folder path to results of analysis (defualt is s2_garnett_classifier_results)
#' @import Seurat
#' @import monocle3
#' @import garnett
#' @keywords cell type classification using garnett classifier 
#' @export

final_cds_garnett_classification <- function(cds_garnett) {

    cds_garnett$cluster_ext_type_final <- cds_garnett$cluster_ext_type
    table <- cbind(
        as.character(cds_garnett$barcode),
        as.character(cds_garnett$cluster_ext_type),
        as.character(cds_garnett$garnett_cluster)
    )
    df <- as.data.frame(table)
    dfunknown <- df[df$V2 == "Unknown", ]
    clusters_unknown <- unique(dfunknown$V3)
    for (cluster in clusters_unknown) {
        temp <- df[df$V3 == cluster, ]
        cell_types <- unique(temp$V2)
        count_values <- c()
        for (cell in cell_types) {
            count_value <- length(which(temp$V2 == cell))
            count_values <- c(count_values, count_value)
        }
        max_value <- max(count_values)
        indexvalue <- match(max_value, count_values)
        final_cluster <- cell_types[indexvalue]
        if (final_cluster != "Unknown") {
            print(paste0("Assigning cell type ", final_cluster, " to cluster ", cluster))
            tempunknown <- temp[temp$V2 == "Unknown", ]
            for (row in rownames(tempunknown)) {
                cds_garnett$cluster_ext_type_final[as.double(row)] <- as.character(final_cluster)
            }
        } else {
            print(
                paste0(
                    "WARNING: Unknown is max classification for this cluster ",
                    cluster, " therefore getting next most popular classification"
                )
            )
            count_values <- count_values[-indexvalue]
            cell_types <- cell_types[-indexvalue]
            max_value <- max(count_values)
            indexvalue <- match(max_value, count_values)
            final_cluster <- cell_types[indexvalue]
            print(paste0("Assigning cell type ", final_cluster, " to cluster ", cluster))
            tempunknown <- temp[temp$V2 == "Unknown", ]
            for (row in rownames(tempunknown)) {
                cds_garnett$cluster_ext_type_final[as.double(row)] <- as.character(final_cluster)
            }
        }
    }
    return(cds_garnett)
}

classify_pbmc_cells <- function(seurat_data, path_10x_data,
                                marker_file_path,
                                pbmc_classifier = FALSE,
                                cds_gene_type = "ENSEMBL",
                                results_folder = "s2_garnett_classifier_results") {

    results_path <- generate_folder(results_folder)
    unlink(paste0(results_folder, "/*"))

    print("STATUS: load cell ranger data")
    SAMPLE <- load_cellranger_data(sample_path, umi_cutoff = 200)

    SAMPLE.counts <- SAMPLE@assays@data@listData$counts
    gene_meta_data <- rowData(SAMPLE)
    cell_meta_data <- colData(SAMPLE)

    cds <- new_cell_data_set(SAMPLE.counts,
        cell_metadata = cell_meta_data,
        gene_metadata = gene_meta_data
    )

    barcodes <- names(seurat_data$percent.mt)
    alteredbarcodes <- c()

    for (i in 1:(length(names(seurat_data$percent.mt)))) {
        alteredbarcodes[i] <- barcodes[i]
    }

    print(paste0(
        "STATUS: Seurat determined that ",
        length(barcodes),
        " should be kept based on Mitochondrial filtering"
    ))

    cds_subset <- cds[, alteredbarcodes]

    print("STATUS: Preprocessing step")
    pbmc_cds <- preprocess_cds(cds_subset,
        num_dim = 100, preprocess_method = "PCA"
    )

    print("STATUS: reduce dimensions")
    rpbmc_cds <- reduce_dimension(pbmc_cds,
        preprocess_method = "PCA",
        reduction_method = "UMAP"
    )
    if (typeof(pbmc_classifier) == "character") {
        pbmc_classifier <- readRDS(pbmc_classifier)
    } else {
        pbmc_classifier <- train_cell_classifier(
            cds = cds,
            marker_file = marker_file_path,
            db = org.Hs.eg.db,
            cds_gene_id_type = cds_gene_type,
            cluster_extend_max_frac_incorrect = 0,
            cluster_extend_max_frac_unknown = 1,
            num_unknown = 500,
            verbose = TRUE,
            marker_file_gene_id_type = "SYMBOL"
        )
    }
    ## Classify cells with classifer built previously
    rpbmc_cds <- classify_cells(rpbmc_cds, pbmc_classifier,
        db = org.Hs.eg.db,
        cluster_extend = TRUE,
        verbose = TRUE,
        cluster_extend_max_frac_incorrect = 0,
        cluster_extend_max_frac_unknown = 1,
        cds_gene_id_type = cds_gene_type
    )
    rpbmc_cds <- final_cds_garnett_classification(rpbmc_cds)


    marker_check <- check_markers(rpbmc_cds, marker_file_path,
        db = org.Hs.eg.db,
        cds_gene_id_type = "ENSEMBL",
        marker_file_gene_id_type = "SYMBOL"
    )
    # Generate marker figure
    plot_markers(marker_check)
    ggsave(paste0(results_path, "Marker_Fig.png"), dpi = 500)


    plot_cells(rpbmc_cds,
        color_cells_by = "cluster_ext_type",
        label_cell_groups = FALSE
    )

    ggsave(paste0(results_path, "Umap_cluster_ext_type.png"),
        width = 5, height = 5, dpi = 500
    )

    plot_cells(rpbmc_cds,
        group_cells_by = "cluster",
        color_cells_by = "cell_type",
        label_cell_groups = FALSE
    )
    ggsave(paste0(results_path, "Umap_cell_type.png"),
        width = 5, height = 5, dpi = 500
    )

    plot_cells(rpbmc_cds,
        color_cells_by = "cluster_ext_type_final",
        label_cell_groups = FALSE
    )
    ggsave(paste0(results_path, "Umap_cluster_ext_type_final.png"),
        width = 5, height = 5,
        dpi = 500
    )

    seurat_data[["garnett_cluster_extend"]] <- rpbmc_cds$cluster_ext_type_final

    return(seurat_data)
}