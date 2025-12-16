library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(xlsx)
library(enrichR)
library(readxl)
library(ggrepel)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(sctransform)
library(scCustomize)
library(pheatmap)
library(openxlsx)
library(tidyr)
library(Matrix)
library(plyr)
library(powsimR)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(enrichplot)
library(ggalluvial)
library(CellChat)



setwd("C:/Users/shenykk/OneDrive - A STAR/Desktop/R projects/sarah")
options(future.globals.maxSize = 2 * 1024^3)

##############################################
#
# Load data and QC
#
##############################################

t1 <- readRDS("trial1UPBliver.rds") #2706 cells

t1[["percent.mt"]] <- PercentageFeatureSet(t1, pattern = "^mt-")
VlnPlot(t1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
t1 <- subset(t1, subset = nFeature_RNA >= 500 & nFeature_RNA <= 12000 & percent.mt <= 10) #2281 cells


##############################################
#
# Set UPB idents for t1 and t2, then merge (merge before normalization)
#
##############################################

UPB_positive <- FetchData(t1, vars = "UPB")[,1] > 0
t1$UPB_label <- ifelse(UPB_positive, "UPB+", "UPB-")

gc()
table(t1$UPB_label)

##############################################
#
# Clustering and cell annotation
#
##############################################

t1 <- SCTransform(t1, vars.to.regress = "percent.mt")
t1 <- RunPCA(t1)
ElbowPlot(t1, ndims = 40) # elbows are at 12 and 18 and 25
t1 <- RunUMAP(t1, dims = 1:18)
t1 <- FindNeighbors(t1, dims = 1:18)
t1 <- FindClusters(t1)


x <- DimPlot_scCustom(t1, label = FALSE, aspect_ratio = 1)
ggsave("ng_umap_clusternumbers1.jpeg", plot = x, width = 5, height = 5, units = "in", dpi = 600)

#t1$UPB_label <- factor(t1$UPB_label, levels = c("UPB-", "UPB+"))
#x <- DimPlot_scCustom(t1, order = TRUE, label = FALSE, aspect_ratio = 1, group.by = "UPB_label", colors_use = c("lightgrey", "black"))
#ggsave("ng_umap_UPBannotation1.jpeg", plot = x, width = 5.5, height = 5.5, units = "in", dpi = 600)

t1.markers <- FindAllMarkers(t1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t <- t1.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t)

x <- FeaturePlot_scCustom(t1, aspect_ratio = 1, 
                          features = c("Pecam1","Emcn", "Adgre1","Cd68", "Ptprc", "Ccr2", "Cd86","Ccr5", "Cd3g", "Il2rb",
                                       "Col1a1","Fap", "Lum", "Fbln1","Kcnj10", "Fabp7", "Bcan", "Ptprz1", "Mpz", "Sox10", 
                                       "Mbp", "Mag", "Slc17a6", "Phox2b", "Isl1", "Snap25"), 
                          num_columns = 5)
ggsave("ng_celltypemarkers_featureplots1.jpeg", plot = x, width = 20, height = 24, dpi = 600, units = "in")

Idents(t1) <- t1$seurat_clusters


new.cluster.ids <- c("Fibroblast","Other neurons","Schwann cells", "Satellite glial cells","Schwann cells",
                     "Monocytes and macrophages", "Monocytes and macrophages", "Schwann cells", "Other neurons",
                     "Fibroblast", "Endothelial cells", 
                     "Monocytes and macrophages", "Other neurons", "Other neurons", "Other neurons", "Other neurons",
                     "Other neurons", "Vascular cells")
names(new.cluster.ids) <- levels(t1)
t1 <- RenameIdents(t1, new.cluster.ids)
t1$celltype <- Idents(object = t1)
DimPlot_scCustom(t1, label = FALSE, aspect_ratio = 1)

#t1$celltype2 <- ifelse(t1$UPB_label == "UPB+", 
#                           "Liver-specific NG neurons", 
#                           as.character(t1$celltype))

t1$celltype2 <- ifelse(t1$UPB_label == "UPB+" & t1$celltype == "Other neurons", 
                       "Liver-specific NG neurons", 
                       as.character(t1$celltype))

Cd3g_positive <- FetchData(t1, vars = "Cd3g")[,1] > 0
t1$Cd3g_label <- ifelse(Cd3g_positive, "Cd3g+", "Cd3g-")
table(t1$Cd3g_label)

t1$celltype2 <- ifelse(t1$Cd3g_label == "Cd3g+", 
                       "T cells", 
                       as.character(t1$celltype2))

ordered_levels <- c(
  "Other neurons",                    # 2573
  "Satellite glial cells",            # 6185
  "Schwann cells",                    # 1231
  "Fibroblast",                       # 1076
  "Monocytes and macrophages",       # 746
  "Endothelial cells",               # 349
  "Vascular cells",                  # 157
  "T cells",                         # 75
  "Liver-specific NG neurons"
)

# Set factor levels
t1$celltype2 <- factor(t1$celltype2, levels = ordered_levels)

Idents(t1) <- t1$celltype2

x <- DimPlot_scCustom(t1, label = FALSE, order = TRUE, aspect_ratio = 1, colors_use = c(
  "#E6E6E6",  # light grey
  "#F13030",  # red
  "#FF00FF",  # magenta
  "#00FF00",  # neon green
  "#3399FF",  # sky blue
  "#FFB733",  # orange
  "#A02C6D",  # dark magenta
  "#33FFF3",  # cyan
  "#4B424A"  # dark grey
))
ggsave("ng_UMAP_celltype1.jpeg", plot = x, width = 6.5, height =6.5, dpi = 600, units = "in")


# dotplot of cell type markers with LSN removed
t12 <- subset(x = t1, idents = c("Liver-specific NG neurons"), invert = TRUE)
celltype.markers <- FindAllMarkers(t12, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
t <- celltype.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
x <- DotPlot_scCustom(seurat_object = t12, features = unique(t$gene), flip_axes = T,
                      x_lab_rotate = TRUE)
ggsave("ng_dotplot_celltype1.jpeg", plot = x, width = 6, height =6, dpi = 600, units = "in")

##############################################
#
# Zhang Na's request for violin and umap plots
#
##############################################


x <- VlnPlot(t1, features = c(
  "Ptger3", "Ptger4", "Gpr65", "Glp1r", "Sst", "Tac1", "Trpv1", "Trpa1", "Oxtr", "Cysltr2",
  "Ahr", "Tlr4", "Il1r1", "Tnfrsf1a", "Tnfrsf1b", "Il6r", "Il6st", "Ghsr", "Th", "Mrgpre",
  "Mrgprf", "Slc2a4", "Insr", "Ackr1", "C1ql2", "Chrnb4", "Igf1r", "Lepr", "Cckbr", "Grp", 
  "Phox2b", "Slc17a6", "Prph", "Nefh", "UPB", "tdTomato"
))

ggsave("zhangna_request_RNA2.jpeg", plot = x, width = 16, height =35, dpi = 300, units = "in")

x <- FeaturePlot_scCustom(t1, aspect_ratio = 1, features = c(
  "Ptger3", "Ptger4", "Gpr65", "Glp1r", "Sst", "Tac1", "Trpv1", "Trpa1", "Oxtr", "Cysltr2",
  "Ahr", "Tlr4", "Il1r1", "Tnfrsf1a", "Tnfrsf1b", "Il6r", "Il6st", "Ghsr", "Th", "Mrgpre",
  "Mrgprf", "Slc2a4", "Insr", "Ackr1", "C1ql2", "Chrnb4", "Igf1r", "Lepr", "Cckbr", "Grp",
  "Phox2b", "Slc17a6", "Prph", "Nefh", "UPB", "tdTomato"
))

ggsave("zhangna_request_SCT_feature2.jpeg", plot = x, width = 16, height =35, dpi = 300, units = "in")



##############################################
#
# DEGs 
#
##############################################
saveRDS(t1, "t1.RDS")

Idents(t1) <- t1$celltype2
markers <- FindMarkers(t1,
                       ident.1 = "Liver-specific NG neurons",
                       ident.2 = "Other neurons",
                       min.pct = 0.1,
                       logfc.threshold = 0)
markers$gene <- rownames(markers)
markers <- markers[!markers$gene %in% c("UPB1", "UPB2", "UPB", "tdTomato"), ]
markers <- markers[!grepl("^mt-", markers$gene, ignore.case = TRUE), ]
markers$fdr <- p.adjust(markers$p_val, method = "BH")

write.xlsx2(markers, "DEGs_trial1.xlsx")

markers$color <- "NS."
markers$color[markers$avg_log2FC > 1 & markers$p_val < 0.05] <- "Up"
markers$color[markers$avg_log2FC < -1 & markers$p_val < 0.05] <- "Down"
markers$neglog10pval <- -log10(markers$p_val)

ptger3_row <- subset(markers, gene == "Ptger3")
color_map <- c("Up" = "red", "Down" = "blue", "NS." = "gray")

x <- ggplot(markers, aes(x = avg_log2FC, y = neglog10pval)) +
  geom_point(aes(color = color), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = color_map) +
  theme_classic(base_size = 14) +
  labs(x = "Log2 Fold Change", y = "-Log10(p_val)", color = "") +
  
  geom_point(data = ptger3_row,
             mapping = aes(x = avg_log2FC, y = neglog10pval),
             shape = 21, fill = "red", color = "black", size = 3, stroke = 1) +
  
  # Add Ptger3 label
  geom_text(data = ptger3_row,
            aes(x = avg_log2FC, y = neglog10pval, label = gene),
            vjust = -1.2, hjust = 0.5,
            size = 4, fontface = "bold") 
ggsave("ng_DEG_volcano1.jpeg", plot = x, width = 5.5, height = 5, dpi = 600, units = "in")




##############################################
#
# To get GO and KEGG terms
#
##############################################

up_genes <- markers %>%
  dplyr::filter(color == "Up", !is.na(avg_log2FC)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Convert to ENTREZ ID
gene_entrez <- clusterProfiler::bitr(up_genes$gene,
                                     fromType = "SYMBOL",
                                     toType = "ENTREZID",
                                     OrgDb = org.Mm.eg.db)

# Merge logFC and ENTREZ
ranked_up <- dplyr::inner_join(up_genes, gene_entrez, by = c("gene" = "SYMBOL"))

# Create named vector for GSEA
gene_list <- ranked_up$avg_log2FC
names(gene_list) <- ranked_up$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_kegg <- clusterProfiler::gseKEGG(geneList = gene_list,
                                      organism = "mmu",
                                      pvalueCutoff = 0.1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      verbose = FALSE, 
                                      pAdjustMethod = "none")

x <- enrichplot::dotplot(gsea_kegg, showCategory = 15, title = "GSEA: KEGG Pathways") +
  ggplot2::theme_minimal()
ggsave("ng_enrichplot_KEGG1.jpeg", plot = x, width = 5, height = 7.5, dpi = 600, units = "in")
write.xlsx2(as.data.frame(gsea_kegg), file = "ng_KEGG_terms1.xlsx")


gsea_go <- clusterProfiler::gseGO(geneList = gene_list,
                                  OrgDb = org.Mm.eg.db,
                                  ont = "BP",
                                  keyType = "ENTREZID",
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = 0.1,
                                  pAdjustMethod = "none",
                                  verbose = FALSE)
x <- enrichplot::dotplot(gsea_go, showCategory = 15, title = "GSEA: GO Biological Process") +
  ggplot2::theme_minimal()
ggsave("ng_enrichplot_GO_BP1.jpeg", plot = x, width = 5.5, height = 6, dpi = 600, units = "in")
write.xlsx2(as.data.frame(gsea_go), file = "ng_GO_BP_terms1.xlsx")



##############################################
#
# Power analysis 
#
##############################################

.roc_plot <- function(ROCData) {
  roc.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ROCData,
                       ggplot2::aes_(x = quote(FPR_Mean),
                                     y = quote(TPR_Mean),
                                     group = quote(Samples),
                                     color = quote(Samples)),
                       size = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         size = 0.75,
                         linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(x = "1 - Specificity (FPR)",
                  y = "Sensitivity (TPR)") +
    .theme_eval_roc()
  
  return(roc.plot)
}
.theme_eval_roc <- function() {
  ggplot2::theme_linedraw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size=10, color='black'),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   legend.key.size = grid::unit(1, "lines"),
                   axis.text = ggplot2::element_text(size=10, color='black'),
                   axis.title = ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title =ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white",
                                                            colour = NA))
}



plotEvalROC <- function(evalRes, cutoff=c('liberal', 'conservative'), Annot=TRUE) {
  
  cutoff = match.arg(cutoff)
  cutoff = ifelse(cutoff == "liberal", "lib", "conv")
  
  # ROC curve
  roc.plot <- .roc_plot(ROCData = evalRes$Performances$`ROC-Curve`)
  
  # PR curve
  pr.plot <- .pr_plot(ROCData = evalRes$Performances$`PR-Curve`)
  
  # TPR vs FDR curve
  alpha.nominal <- as.numeric(evalRes$Settings["alpha.nominal"])
  tprvsfdr.plot <-  .tprvsfdr_plot(ROCData = evalRes$TPRvsFDR,
                                   alpha.nominal = alpha.nominal)
  
  # table with summary statistics
  summary.tbl <- .summary_table_print(TblData = evalRes$Scores,
                                      cutoff = cutoff)
  
  
  curve.legend <- cowplot::get_legend(tprvsfdr.plot)
  # combine the plots
  top_row <- suppressWarnings(cowplot::plot_grid(roc.plot + ggplot2::theme(legend.position = "none"),
                                                 pr.plot + ggplot2::theme(legend.position = "none"),
                                                 labels= LETTERS[1:2],
                                                 ncol=2, nrow=1))
  bottom_row <- suppressWarnings(cowplot::plot_grid(tprvsfdr.plot + ggplot2::theme(legend.position = "none"),
                                                    summary.tbl,
                                                    labels= LETTERS[3:4],
                                                    ncol=2, nrow=1,
                                                    rel_widths = c(0.4, 0.6)))
  
  p.combined <- suppressWarnings(cowplot::plot_grid(top_row,
                                                    bottom_row,
                                                    ncol=1, nrow=2))
  
  p.final <- suppressWarnings(cowplot::plot_grid(p.combined,
                                                 curve.legend, rel_heights = c(1,0.1),
                                                 ncol=1, nrow=2))
  
  # annotation under plot
  if(Annot) {
    delta.text <- ifelse(evalRes$Settings$delta != 0,
                         paste0("considering genes with at least ", evalRes$Settings$delta, " as biologically meaningful DE genes.\n"),
                         " considering all DE genes.\n")
    settings.text <- paste0(evalRes$Settings$alpha.type, " p-values with nominal level equal to ", evalRes$Settings$alpha.nominal, delta.text)
    annot.text <- c("A) Receiver-Operator-Characteristics (ROC) Curve per sample size setup. \nB) Precision-Recall (PR) Curve per sample size setup. \nC) TPR versus observed FDR per sample size setup. The filling of the point indicates whether FDR is controlled at the chosen nominal level. \nD) Summary Statistics per sample size setup rounded to two digits.")
    p.final <- cowplot::add_sub(plot = p.final,
                                label = paste0(settings.text, annot.text),
                                x = 0, hjust = 0, size=8)
  }
  
  # draw final plot
  cowplot::ggdraw(p.final)
}

.pr_plot <- function(ROCData) {
  pr.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ROCData,
                       ggplot2::aes_(x = quote(TPR_Mean),
                                     y = quote(PPV_Mean),
                                     group = quote(Samples),
                                     color = quote(Samples)),
                       size = 1) +
    ggplot2::scale_x_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0,1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(x = "Recall (TPR)", y = "Precision (PPV)") +
    .theme_eval_roc()
  
  return(pr.plot)
}
# tpr vs fdr curve
.tprvsfdr_plot <- function(ROCData, alpha.nominal) {
  tprvsfdr.nominal <- ROCData %>%
    dplyr::filter(.data$Threshold == alpha.nominal) %>%
    dplyr::mutate(`FDR Control` = ifelse(.data$FDR_Mean <= .data$Threshold, "yes", "no"),
                  FDRLower = .data$FDR_Mean - .data$FDR_SE,
                  FDRUpper = .data$FDR_Mean + .data$FDR_SE,
                  TPRLower = .data$TPR_Mean - .data$TPR_SE,
                  TPRUpper = .data$TPR_Mean + .data$TPR_SE)
  
  upperlimit <- ifelse(max(tprvsfdr.nominal$FDR_Mean) <= alpha.nominal,
                       max(tprvsfdr.nominal$FDR_Mean)*5,
                       max(tprvsfdr.nominal$FDR_Mean)*2)
  
  tprvsfdr.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ROCData,
                       ggplot2::aes_(x = quote(FDR_Mean),
                                     y = quote(TPR_Mean),
                                     group = quote(Samples),
                                     color = quote(Samples)),
                       size = 1) +
    ggplot2::geom_point(data = tprvsfdr.nominal,
                        ggplot2::aes_(x = quote(FDR_Mean),
                                      y = quote(TPR_Mean),
                                      group = quote(Samples),
                                      shape = quote(`FDR Control`)),
                        size = 3) +
    ggplot2::geom_linerange(data = tprvsfdr.nominal,
                            ggplot2::aes_(ymin=quote(TPRLower),
                                          ymax=quote(TPRUpper),
                                          x = quote(FDR_Mean),
                                          group = quote(Samples),
                                          color = quote(Samples)),
    ) +
    ggstance::geom_linerangeh(data = tprvsfdr.nominal,
                              ggplot2::aes_(xmin=quote(FDRLower),
                                            xmax=quote(FDRUpper),
                                            y = quote(TPR_Mean),
                                            group = quote(Samples),
                                            color = quote(Samples))) +
    ggplot2::scale_fill_manual(values = ) +
    ggplot2::scale_x_continuous(limits = c(0,upperlimit),
                                expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(y = "TPR",
                  x = "observed FDR") +
    ggplot2::guides() +
    .theme_eval_roc()
  
  return(tprvsfdr.plot)
}


.summary_table_print <- function(TblData, cutoff) {
  
  fixed.scores <- c("TPRvsFPR_AUC", "TPRvsPPV_AUC")
  select.scores <- paste(c("ACC", "MCC", "F1score", "TPRvsFDR_pAUC"), cutoff, sep="_")
  scores <- c(select.scores, fixed.scores)
  
  tbl.dat <- TblData %>%
    dplyr::filter(.data$Score %in% scores) %>%
    dplyr::mutate(`Summary Statistic` = gsub(pattern = paste0("_", cutoff), replacement = "", .data$Score)) %>%
    dplyr::mutate(`Summary Statistic` = gsub(pattern = "_", replacement = " ", .data$`Summary Statistic`)) %>%
    dplyr::mutate(`Summary Statistic` = case_when(.data$`Summary Statistic` == "TPRvsFPR AUC" ~ "ROC AUC",
                                                  .data$`Summary Statistic` == "TPRvsPPV AUC" ~ "PRC AUC",
                                                  TRUE ~ as.character(.data$`Summary Statistic`))) %>%
    tidyr::separate(.data$Samples, c('n1', 'n2'), " vs ", remove=FALSE) %>%
    dplyr::mutate(SumN = as.numeric(.data$n1)+as.numeric(.data$n2),
                  Mean = round(.data$Mean, digits = 2),
                  SE = round(.data$SE, digits = 2)) %>%
    dplyr::arrange(.data$`Summary Statistic`, .data$SumN) %>%
    tidyr::unite("Value", c("Mean", "SE"), sep = "\u00B1") %>%
    dplyr::select(.data$`Summary Statistic`, .data$Samples, .data$Value) %>%
    tidyr::pivot_wider(id_cols = "Samples", names_from = "Summary Statistic", values_from = "Value")
  
  summary.tbl <- ggpubr::ggtexttable(tbl.dat,
                                     rows = NULL,
                                     theme = ggpubr::ttheme("lBlackWhite"))
  
  return(summary.tbl)
}

countData <- as.data.frame(GetAssayData(t1, slot = "counts"))
countData <- as.matrix(countData)
gene_vars <- apply(countData, 1, var)
countData <- countData[gene_vars > 0, ]
countData <- countData[complete.cases(countData), ]
# This takes a long time (~2h)
estparam_gene <- estimateParam(countData,
                               readData = NULL,
                               batchData = NULL,
                               spikeData = NULL,
                               spikeInfo = NULL,
                               Lengths = NULL,
                               MeanFragLengths = NULL,
                               RNAseq = 'singlecell',
                               Protocol = 'UMI',
                               Distribution = 'NB',
                               Normalisation = c("TMM"),
                               GeneFilter = 0.05,
                               SampleFilter = 5,
                               sigma = 3.96,
                               NCores = NULL,
                               verbose=TRUE)

plotParam(estParamRes = estparam_gene, Annot = T)
x<- plotParam(estParamRes = estparam_gene, Annot = T)
ggsave("powsimr_t1.jpeg", plot = x, width =16 , height = 8, dpi = 600, units = "in")


# define log fold change
p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 2, rate = 0.5)

# set up simulations
setupres <- Setup(ngenes = NULL, nsims = 5, #nsims is the number of simulations
                  p.DE = 0.1, pLFC = 1, p.G = 1, # change these?
                  #n1 = c(10, 20, 50, 100, 200), n2 = c(10, 20, 50, 100, 200), #change cell number for real thing
                  Thinning = NULL, LibSize = 'equal',
                  estParamRes = estparam_gene,
                  estSpikeRes = NULL,
                  DropGenes = TRUE,
                  setup.seed = 5299, verbose = TRUE)

simres <- simulateDE(SetupRes = setupres,
                     Prefilter = NULL, Imputation = NULL,
                     Normalisation = 'scran',
                     DEmethod = "limma-trend", DEFilter = FALSE,
                     NCores = NULL, verbose = TRUE)

evalderes = evaluateDE(simRes = simres,
                       alpha.type = 'adjusted',
                       MTC = 'BH',
                       alpha.nominal = 0.1,
                       stratify.by = 'dispersion', #mean, dispersion, dropout, lfc
                       filter.by = 'dispersion',
                       strata.filtered = 1,
                       target.by = 'lfc',
                       delta = 0)

plotEvalDE(evalRes = evalderes, rate = 'marginal', quick = TRUE, Annot = TRUE)
plotEvalDE(evalRes = evalderes, rate = 'conditional', quick = TRUE, Annot = TRUE)

x<- plotEvalDE(evalRes = evalderes, rate = 'marginal', quick = TRUE, Annot = TRUE)
ggsave("powsimr_t1_2.jpeg", plot = x, width =8 , height = 8, dpi = 600, units = "in")

x<- plotEvalDE(evalRes = evalderes, rate = 'conditional', quick = TRUE, Annot = TRUE)
ggsave("powsimr_t1_3.jpeg", plot = x, width =8 , height = 8, dpi = 600, units = "in")


##############################################
#
# Liver cell annotation
#
##############################################

liver <- readRDS("Mouseliverstst.rds")
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^mt-")
liver <- subset(liver, subset = nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt <= 10) #282453 cells
liver <- liver[, sample(colnames(liver), size =20000, replace=F)] #downsample to 20k cells

liver <- SCTransform(liver, vars.to.regress = "percent.mt")
liver <- RunPCA(liver)
ElbowPlot(liver, ndims = 40) # elbows are at 14 and 22
liver <- RunUMAP(liver, dims = 1:22)
liver <- FindNeighbors(liver, dims = 1:22)
liver <- FindClusters(liver, resolution = 0.4)

x <- DimPlot_scCustom(liver, label = FALSE, aspect_ratio = 1)
ggsave("liver_umap_clusternumbers.jpeg", plot = x, width = 6, height = 6, units = "in", dpi = 600)

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t <- liver.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t)

x <- FeaturePlot_scCustom(liver, aspect_ratio = 1, 
                          features = c("Alb","Apoa1", "Ttr","Cyp3a11", "Clu", "Aldob", "Pck1","Uox", "Ass1", "Gls2", "Cps1",
                                       "Dcn","Ecm1", "Lrat", "Reln","Rbp1", "Clec4f", "Cd5l", "Timd4", "Plpp3", "Igfbp7", 
                                       "Clec4g", "Pecam1", "Stab2", "Lyve1", "Fcgr2b", "Krt19", "Sox9", "Epcam", "Tspan8"), 
                          num_columns = 5)

x <- FeaturePlot_scCustom(liver, aspect_ratio = 1, 
                          features = c( 
                            "Ccr2", "Lyz1", "Cx3cr1", "Ace", "Fn1", "Cd9", "Col1a1", "Msln", "Pdgfra", "Pdgfrb", "Lum", 
                            "Tcf7", "Il7r", "Cd3e", "Cd3g", "Ms4a1", "Igkc", "Ighm", "Retnlg", "Ly6g", "Cxcr2", "S100a9", "Siglech", 
                            "Flt3", "Zbtb46", "Itgax", "Ncr1", "Nkg7", "Gzma", "Ccl5"), 
                          num_columns = 5)
ggsave("liver_celltypemarkers_featureplots2.jpeg", plot = x, width = 20, height = 20, dpi = 600, units = "in")

liver <- subset(x = liver, idents = c("21", "22", "16"), invert = TRUE) # removing unknown and nerve

# remove those clusters and repeat. Cluster 16 is nerves
liver <- SCTransform(liver, vars.to.regress = "percent.mt")
liver <- RunPCA(liver)
ElbowPlot(liver, ndims = 40) # elbows are at 14 and 22
liver <- RunUMAP(liver, dims = 1:21)
liver <- FindNeighbors(liver, dims = 1:21)
liver <- FindClusters(liver, resolution = 0.5)

x <- DimPlot_scCustom(liver, label = FALSE, aspect_ratio = 1)
ggsave("liver_umap_clusternumbers.jpeg", plot = x, width = 6, height = 6, units = "in", dpi = 600)

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t <- liver.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t)

x <- FeaturePlot_scCustom(liver, aspect_ratio = 1, 
                          features = c("Alb","Apoa1", "Ttr","Cyp3a11", "Clu", "Aldob", "Pck1","Uox", "Ass1", "Gls2", "Cps1",
                                       "Dcn","Ecm1", "Lrat", "Reln","Rbp1", "Clec4f", "Cd5l", "Timd4", "Plpp3", "Igfbp7", 
                                       "Clec4g", "Pecam1", "Stab2", "Lyve1", "Fcgr2b", "Krt19", "Sox9", "Epcam", "Tspan8"), 
                          num_columns = 5)

x <- FeaturePlot_scCustom(liver, aspect_ratio = 1, 
                          features = c( 
                            "Ccr2", "Lyz1", "Cx3cr1", "Ace", "Fn1", "Cd9", "Col1a1", "Msln", "Pdgfra", "Pdgfrb", "Lum", 
                            "Tcf7", "Il7r", "Cd3e", "Cd3g", "Ms4a1", "Igkc", "Ighm", "Retnlg", "Ly6g", "Cxcr2", "S100a9", "Siglech", 
                            "Flt3", "Zbtb46", "Itgax", "Ncr1", "Nkg7", "Gzma", "Ccl5"), 
                          num_columns = 5)
ggsave("liver_celltypemarkers_featureplots2.jpeg", plot = x, width = 20, height = 20, dpi = 600, units = "in")

liver <- subset(x = liver, idents = c("22", "23"), invert = TRUE)

Idents(liver) <- liver$seurat_clusters
new.cluster.ids <- c("Liver sinusoid epithelial","Kupffer", "Liver sinusoid epithelial","Monocyte and macrophages",
                     "Kupffer", "Periportal hepatocytes", "Monocyte and macrophages", "Dendritic cells", "T cells",
                     "Periportal hepatocytes", "Hepatocytes",  "NK cells", "B cells", "Hepatocytes","Liver sinusoid epithelial",
                     "Liver sinusoid epithelial", "Kupffer", "Neutrophils", "Liver sinusoid epithelial", "Dividing cells", 
                     "Hepatic stellate cells", "Cholangiocytes")

names(new.cluster.ids) <- levels(liver)
liver <- RenameIdents(liver, new.cluster.ids)
liver$celltype <- Idents(object = liver)
DimPlot_scCustom(liver, label = FALSE, aspect_ratio = 1)

t1$celltype2 <- ifelse(t1$UPB_label == "UPB+", 
                           "Liver-specific NG neurons", 
                           as.character(t1$celltype))

ordered_levels <- c(
  "Hepatocytes",                    # 2573
  "Periportal hepatocytes",            # 6185
  "Hepatic stellate cells",                    # 1231
  "Kupffer",
  "Liver sinusoid epithelial",
  "Cholangiocytes", # 1076
  "Monocyte and macrophages",       # 746
  "T cells",               # 349
  "B cells",                  # 157
  "Neutrophils",                        # 75
  "Dendritic cells",
  "NK cells", 
  "Dividing cells"
)

# Set factor levels
liver$celltype <- factor(liver$celltype, levels = ordered_levels)

Idents(liver) <- liver$celltype

x <- DimPlot_scCustom(liver, label = FALSE, order = TRUE, aspect_ratio = 1)
ggsave("liver_UMAP_celltype.jpeg", plot = x, width = 8, height =8, dpi = 600, units = "in")


# dotplot of cell type markers

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.3)
liver.markers <- dplyr::filter(liver.markers, !grepl("^Gm", gene))
t <- liver.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
x <- DotPlot_scCustom(seurat_object = liver, features = unique(t$gene), flip_axes = T,
                      x_lab_rotate = TRUE)
ggsave("liver_dotplot_celltype.jpeg", plot = x, width = 6, height =6, dpi = 600, units = "in")

# ask Zhang Na if need to differentiate between periportal and pericentral hepatocytes

Idents(liver) <- plyr::mapvalues(
  Idents(liver),
  from = "Periportal hepatocytes",
  to = "Hepatocytes"
)

saveRDS(liver, "liver.RDS")
liver <- readRDS("liver.RDS")

##############################################
#
# Cell type proportion graph
#
##############################################

cell_counts <- table(Idents(liver))
cell_df <- as.data.frame(cell_counts)
colnames(cell_df) <- c("celltype", "count")

cell_df$proportion <- cell_df$count / sum(cell_df$count)

color_map <- c(
  "Hepatocytes" = "#4D4D4D",        # dark grey
  "Hepatic stellate cells" = "#E6E6E6",  # light grey
  "Kupffer" = "red",            # hot pink
  "Liver sinusoid epithelial" = "magenta",  # bright green
  "Cholangiocytes" = "#00FF00",     # cyan
  "Monocyte and macrophages" = "#3399FF", # light blue
  "T cells" = "#FFB300",            # orange
  "B cells" = "maroon",            # red
  "Neutrophils" = "#66F7F1",        # teal
  "Dendritic cells" = "#808000",    # darker teal
  "NK cells" = "#66CCFF",           # turquoise
  "Dividing cells" = "#FFB3FF"      # light pink
)



cell_counts <- table(Idents(t1))
cell_df <- as.data.frame(cell_counts)
colnames(cell_df) <- c("celltype", "count")

cell_df$proportion <- cell_df$count / sum(cell_df$count)

color_map <- c(
  "Other neurons" = "#E6E6E6",  # light grey
  "Satellite glial cells" = "red",            # hot pink
  "Schwann cells" = "magenta",  # bright green
  "Fibroblast" = "#00FF00",     # cyan
  "Monocytes and macrophages" = "#3399FF", # light blue
  "Endothelial cells" = "#FFB300",            # orange
  "Vascular cells" = "maroon",            # red
  "T cells" = "#66F7F1",        # teal
  "Liver-specific NG neurons" = "#4D4D4D"      # light pink
)

ggplot(cell_df, aes(x = "All Cells", y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = color_map) +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))



##############################################
#
# Cellchat
#
##############################################

# For this analysis include only: Hepatocytes, 
# Hepatic stellate cells, Kupffer, Liver sinusoid epithelial, 
# Cholangiocytes, Monocyte and macrophages, LSN
liver <- readRDS("liver.RDS")
t1 <- readRDS("t1.RDS")

Idents(t1) <- t1$celltype2
t12 <- subset(x = t1, idents = c("Liver-specific NG neurons"))
table(Idents(t12))
x <- t12
y <- subset(liver, idents = c("Hepatocytes",
                              "Hepatic stellate cells", 
                              "Kupffer",
                              "Liver sinusoid epithelial",
                              "Cholangiocytes",
                              "Monocyte and macrophages"))
cc <- merge(x, y = y)
table(Idents(cc))
cc <- SCTransform(cc, vars.to.regress = "percent.mt")




# https://rdrr.io/github/sqjin/CellChat/f/tutorial/CellChat-vignette.Rmd
# Convert Seurat objects to CellChat objects
cellchat1 <- createCellChat(object = cc, group.by = "ident")
CellChatDB <- CellChatDB.mouse
cellchat1@DB <- CellChatDB

cellchat1 <- subsetData(cellchat1)
cellchat1 <- identifyOverExpressedGenes(cellchat1) 
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- computeCommunProb(cellchat1) # takes 1h 
cellchat1 <- filterCommunication(cellchat1) 
cellchat1 <- computeCommunProbPathway(cellchat1) 

cellchat1 <- updateCellChat(cellchat1)
cellchat1 <- aggregateNet(cellchat1)

saveRDS(cellchat1, "cellchat_new.RDS") # use this cellchat_LSN82_ON2573.RDS

cellchat1 <- readRDS("cellchat_new.RDS")

groupSize <- as.numeric(table(cellchat1@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# outgoing
mat <- cellchat1@net$count
par(mfrow = c(2, 4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#incoming
mat <- cellchat1@net$weight
par(mfrow = c(2, 4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[ ,i] <- mat[ ,i]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_heatmap(cellchat1)
netVisual_heatmap(cellchat1, measure = "weight")


netVisual_bubble(cellchat1)
netVisual_bubble(cellchat1, sources.use = c("Liver-specific NG neurons"), 
                 remove.isolate = FALSE)
netVisual_bubble(cellchat1, targets.use = c("Liver-specific NG neurons"), 
                 remove.isolate = FALSE)

cellchat1@netP$pathways #all

# neurons only
neuron_sending <- subsetCommunication(cellchat1, sources.use =c("Liver-specific NG neurons"))
neuron_sending_pathways <- unique(neuron_sending$pathway_name)
neuron_receiving <- subsetCommunication(cellchat1, targets.use = c("Liver-specific NG neurons"))
neuron_receiving_pathways <- unique(neuron_receiving$pathway_name)
neuron_combined <- unique(c(neuron_sending_pathways, neuron_receiving_pathways))
print(neuron_combined)

for (pathway in neuron_combined) {
  message("Plotting pathway: ", pathway)
  filename <- paste0("cellchat_", gsub(" ", "_", pathway), ".jpeg")
  jpeg(filename, width = 1500, height = 1550, res = 300)
  #par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))  # 1 row, 2 columns; top margin allows for title
  netVisual_aggregate(cellchat1, signaling = pathway, layout = "circle", arrow.size = 1)
  title(main = pathway, line = 2, font.main = 2)
  dev.off()
}


cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")

for (pathway in neuron_combined) {
  message("Plotting pathway: ", pathway)
  output_file <- paste0("cellchat_role_", gsub(" ", "_", pathway), ".pdf")
  pdf(output_file, width = 8, height = 8)
  netAnalysis_signalingRole_network(cellchat1, signaling = pathway)
  dev.off()
}

netAnalysis_signalingRole_scatter(cellchat1)

