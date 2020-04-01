#' perform differential gene expression
#'
#' This generates new folder (internal function only)
#' @param pbmc filtered data from seurat
#' @param ident_1 cell condition 1
#' @param ident_2 cell condition 2
#' @param LFC log fold cutoff (defualt 0.26)
#' @param pvalue pvalue (FDR cutoff (default 0.05)
#' @param clusters_DE perform DE on individual clusters in ident_1 & 2 (defualt is TRUE, needs s2_cellclasify_garnett to have been run)
#' @param individual_gene_plots generate individual violin plots for every DE gene (defualt is FALSE)
#' @param results_folder path to results of analysis (defualt is s2_garnett_classifier_results)
#' @import Seurat
#' @import org.Hs.eg.db
#' @import dplyr
#' @keywords perform differential gene analysis
#' @export
#'

s3_DE_analysis <- function(pbmc, ident_1, ident_2, LFC = 0.26,
                           pvalue = 0.05, clusters_DE = FALSE,
                           individual_gene_plots = FALSE,
                           results_folder = "s3_DE_results") {
    print("STATUS: performing DE analysis")

    results_path <- generate_folder(results_folder)
    results_path <- generate_folder(paste0(results_folder, ident_1, "_", ident_2))
    unlink(paste0(results_path, "/*"))

    Idents(pbmc) <- pbmc$sample
    print (Idents(pbmc))
    print (ident_1)
    print (ident_2)
    im.markers <- get_DE_between_conditions(
        pbmc,
        ident_1, ident_2,
        results_path,
        LFC = LFC,
        pvalue = pvalue
    )

    if (isTRUE(clusters_DE)) {
        pbmc$celltype.stim <- paste(Idents(pbmc),
            pbmc$garnett_cluster_extend,
            sep = "-"
        )
        Idents(pbmc) <- "celltype.stim"
        allmarkers <- c()
        unique_clusters <- unique(pbmc$garnett_cluster_extend)
        for (cluster in unique_clusters) {
            count_cells_1 <- length(which(Idents(pbmc) == paste0(ident_1, "-", cluster)))
            count_cells_2 <- length(which(Idents(pbmc) == paste0(ident_2, "-", cluster)))

            if ((count_cells_1 > 3) && (count_cells_2 > 3)) {
                im.markers <- get_DE_between_conditions(
                    pbmc,
                    paste0(ident_1, "-", cluster),
                    paste0(ident_2, "-", cluster),
                    results_path,
                    LFC = LFC,
                    pvalue = pvalue
                )
                if (dim(im.markers)[1] != 0) {
                    allmarkers <- append(allmarkers,
                        rownames(im.markers),
                        after = length(allmarkers)
                    )
                    immune.small <- subset(pbmc,
                        idents = c(
                            paste0(ident_1, "-", cluster),
                            paste0(ident_2, "-", cluster)
                        )
                    )

                    g <- DoHeatmap(immune.small,
                        size = 2.5,
                        features = rownames(im.markers)
                    )
                    ggsave(paste0(results_path, "heatmap_", cluster, ".png"),
                        width = 11, dpi = 300
                    )

                    VlnPlot(immune.small, features = head(rownames(im.markers)))

                    ggsave(paste0(
                        results_path, "Vnplot_",
                        cluster, ".png"
                    ),
                    width = 11, dpi = 300
                    )
                }
            } else {
                print(paste0(
                    "STATUS: not enough cells in cluster ",
                    cluster,
                    " ident 1 or 2 to do DE analysis"
                ))
            }
        }
        if (individual_gene_plots == TRUE) {
            results_path <- generate_folder(paste0(result_path, "IndividualGenePlots/"))
            # unlink(paste0(results_folder, "/*"))
            for (gene in unique(allmarkers)) {
                VlnPlot(pbmc,
                    features = c(gene),
                    split.by = "orig.ident",
                    group.by = "garnett_cluster_extend",
                    pt.size = 0, combine = FALSE
                )
                ggsave(paste0(results_path, gene, "_Vnplot.png"), dpi = 300)
            }
        }
    }
}

get_DE_between_conditions <- function(pbmc, ident_1, ident_2,
                                      results_folder, LFC = 0.26,
                                      pvalue = 0.05) {
    print("STATUS: getting DEs...")
    ### Parameters for FindMarkers
    #### test.use: mast (default is wilcox)
    #### min.pct: 0.1(default: 0.1) filter out genes (features) that are detected at less than 10 percent frequency in cells in ident_1 or ident_2
    #### logfc: 0 (default is .25) logfc must be higher than 0 (I set this at 0 because I filter this out at a later stage - this was primarily done to
    ####                          understand how the tool (Seurat) works
    #### min.cells.feature: 3 (default is 3) minimum number of cells expressing the feature in at least one of the two groups (similar to min.pct)
    #### min.cells.group: 3 (defualt is 3) minimum number of cells in the group
    #### max.cells.per.ident: Inf (default Inf-means no down sampling) Down sample each identity class (cluster of cells) to a max number
    #### min.dif.pct: Inf (default Inf) Only test genes that show minimum difference in the fraction of detection between the two identities (cluster)

    DEgenes <- FindMarkers(pbmc,
        ident.1 = ident_1,
        ident.2 = ident_2,
        test.use = "MAST",
        logfc.threshold = 0
    )
    write.table(data.frame(DEgenes),
        paste0(
            results_path,
            "DEgenes_full_",
            ident_1,
            "_", ident_2,
            ".txt"
        ),
        sep = "\t", quote = FALSE
    )
    print(dim(DEgenes))
    print(length(DEgenes))
    DEgenes[["gene_name"]] <- rownames(DEgenes)

    DE_sig_final <- DEgenes %>%
        filter(avg_logFC >= LFC | avg_logFC <= -LFC) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
    DE_sig_final <- DE_sig_final %>%
        filter(p_val_adj <= pvalue) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    if (length(DE_sig_final$gene_name) != 0) {
        rownames(DE_sig_final) <- DE_sig_final$gene_name
        DE_sig_final$gene_name <- NULL
    } else {
        print("No significant DE genes")
    }
    print(dim(DE_sig_final))
    write.table(data.frame(DE_sig_final),
        paste0(
            results_path,
            "DEgenes_sig_",
            ident_1, "_", ident_2,
            ".txt"
        ),
        sep = "\t", quote = FALSE
    )
    return(DE_sig_final)
}
