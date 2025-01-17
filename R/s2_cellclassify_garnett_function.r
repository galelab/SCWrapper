#' Classify PBMC SC data
#'
#' This generates new folder (internal function only)
#' @param seurat_data filtered data from seurat
#' @param path_10x_data path to 10x folder
#' @param lower_gene_filter_num minimum number of genes that has to be detected in a cell
#' @param marker_file_path path to the marker file needed for garnett classification
#' @param classifier path to pre built classifier if one is being used (defualt is to build its own classifier)
#' @param cds_gene_type type of gene ID (defualt is Ensembl)
#' @param marker_file_gene_type type of gene ID in marker file (marker_file_path) (defualt is SYMBOL)
#' @param save_object save seurat normalized object in RDS file so this step doesn't have to be redone everytime (defualt is TRUE)
#' @param results_folder path to results of analysis (defualt is s2_garnett_classifier_results)
#' @import Seurat
#' @import monocle3
#' @import garnett
#' @import dplyr
#' @import stringr
#' @import SingleCellExperiment
#' @import data.table
#' @keywords cell type classification using garnett classifier
#' @export

s2_classify_cells <- function(seurat_data, path_10x_data, sample_names = FALSE,
                                marker_file_path, lower_gene_filter_num = 200,
                                classifier = FALSE,
                                cds_gene_type = "ENSEMBL",
                                marker_file_gene_type="SYMBOL",
                                downsampled = TRUE,
                                save_object = TRUE,
                                results_folder = "s2_garnett_classifier_results") {
    results_path <- generate_folder(results_folder)
    unlink(paste0(results_folder, "/*"))

    print("STATUS: load cell ranger data")
    if (isTRUE(downsampled)) {
        SAMPLE <- load_cellranger_data(path_10x_data,
            umi_cutoff = lower_gene_filter_num
        )

        SAMPLE.counts <- SAMPLE@assays@data@listData$counts
        gene_meta_data <- rowData(SAMPLE)
        cell_meta_data <- colData(SAMPLE)

        cds <- new_cell_data_set(SAMPLE.counts,
            cell_metadata = cell_meta_data,
            gene_metadata = gene_meta_data
        )

    } else {
        counter <- 1
        pbmc_list <- list()
        barcode_list <- list()
        for (path in path_10x_data) {
            print(paste0("STATUS: Loading ", path, " ", sample_names[counter]))
            SAMPLE <- load_cellranger_data(path,
                umi_cutoff = lower_gene_filter_num
            )

            SAMPLE.counts <- SAMPLE@assays@data@listData$counts
            gene_meta_data <- rowData(SAMPLE)
            cell_meta_data <- colData(SAMPLE)

            cds <- new_cell_data_set(SAMPLE.counts,
                cell_metadata = cell_meta_data,
                gene_metadata = gene_meta_data
            )
            alteredbarcodes <- c()
            barcodes <- colnames(cds)
            for (i in 1:(length(barcodes))) {
                alteredbarcodes[i] <- str_replace(barcodes[i], "1", as.character(counter))
            }
            colnames(cds) <- alteredbarcodes
            pbmc_list[[sample_names[counter]]] <- cds
            barcode_list[[sample_names[counter]]] <- alteredbarcodes
            counter <- counter + 1
        }
        cds <- combine_cds(pbmc_list)
    }
    ###Adjust cell barcodes to match Seurat formating
    cds_new_barcodes <- c()
    cds_barcodes <- colnames(cds)
    for (i in 1:(length(cds_barcodes))) {
        cds_new_barcodes[i] <- str_replace(cds_barcodes[i], "_\\S+$", "")
    }
    colnames(cds) <- cds_new_barcodes

    ##Adjust Seurat barcodes
    barcodes <- names(seurat_data$percent.mt)
    alteredbarcodes <- c()

    for (i in 1:(length(names(seurat_data$percent.mt)))) {
        alteredbarcodes[i] <- str_replace(barcodes[i], "_", "-")
    }

    print(paste0(
        "STATUS: Seurat determined that ",
        length(barcodes),
        " should be kept based on Mitochondrial filtering"
    ))
    # print(colnames(cds))
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
    if (typeof(classifier) == "character") {
        pbmc_classifier <- readRDS(classifier)
    } else {
        pbmc_classifier <- train_cell_classifier(
            cds = cds,
            marker_file = marker_file_path,
            db = org.Hs.eg.db,
            cds_gene_id_type = cds_gene_type,
            num_unknown = 500,
            marker_file_gene_id_type = marker_file_gene_type
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
    print(unique(rpbmc_cds$cell_type))
    rpbmc_cds <- final_cds_garnett_classification(rpbmc_cds)

    marker_check <- check_markers(rpbmc_cds, marker_file_path,
        db = org.Hs.eg.db,
        cds_gene_id_type = cds_gene_type,
        marker_file_gene_id_type = marker_file_gene_type
    )

    # Generate marker figure
    plot_markers(marker_check)
    ggsave(paste0(results_path, "/Marker_Fig.png"), dpi = 500)

    visualize_data(rpbmc_cds, results_path, "all")

    seurat_data[["garnett_cell_type"]] <- rpbmc_cds$cell_type
    seurat_data[["garnett_cluster_extend"]] <- rpbmc_cds$cluster_ext_type_final

    if (isTRUE(downsampled)) {
        samps <- as.character(unique(seurat_data$sample))
        barcodes <- as.character(rpbmc_cds$barcode)
        for (counter in 1:(length(samps))) {
            sample_barcodes <- c()
            for (b in 1:(length(barcodes))) {
                if (endsWith(as.character(barcodes[b]), as.character(counter))) {
                    sample_barcodes[b] <- barcodes[b]
                }
            }
            overlap_vector <- intersect(sample_barcodes, colnames(rpbmc_cds))
            rpbmc_cds_sub <- rpbmc_cds[, as.character(overlap_vector)]
            visualize_data(rpbmc_cds_sub, results_path, samps[counter])
            counter <- counter + 1
        }
    } else {
        counter <- 1
        for (barcodes in barcode_list) {
            overlap_vector <- intersect(as.character(barcodes), colnames(rpbmc_cds))
            rpbmc_cds_sub <- rpbmc_cds[, as.character(overlap_vector)]
            visualize_data(rpbmc_cds_sub, results_path, sample_names[counter])
            counter <- counter + 1
        }
    }

    if (isTRUE(save_object)) {
        saveRDS(seurat_data, paste0(results_path, "/SC_pbmc_norm_object.rds"))
    }

    return(seurat_data)
}

visualize_data <- function(rpbmc_cds, results_path, sample) {
    plot_cells(rpbmc_cds,
        color_cells_by = "cluster_ext_type",
        label_cell_groups = FALSE
    )

    ggsave(paste0(results_path, "/Umap_cluster_ext_type_", sample, ".png"),
        width = 6, height = 5, dpi = 500
    )

    plot_cells(rpbmc_cds,
        group_cells_by = "cluster",
        color_cells_by = "cell_type",
        label_cell_groups = FALSE
    )
    ggsave(paste0(results_path, "/Umap_cell_type_", sample, ".png"),
        width = 6, height = 5, dpi = 500
    )

    plot_cells(rpbmc_cds,
        color_cells_by = "cluster_ext_type_final",
        label_cell_groups = FALSE
    )
    ggsave(paste0(results_path,
             "/Umap_cluster_ext_type_final_",
              sample,
              ".png"),
        width = 6, height = 5,
        dpi = 500
    )
}

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
            print(paste0("STATUS: Assigning cell type ", final_cluster, " to cluster ", cluster))
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
            print(paste0("STATUS: Assigning cell type ", final_cluster, " to cluster ", cluster))
            tempunknown <- temp[temp$V2 == "Unknown", ]
            for (row in rownames(tempunknown)) {
                cds_garnett$cluster_ext_type_final[as.double(row)] <- as.character(final_cluster)
            }
        }
    }
    return(cds_garnett)
}