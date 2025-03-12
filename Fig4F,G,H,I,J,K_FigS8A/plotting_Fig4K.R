#Plot figure panel Fig4K + additional plots when DEGs are selected on Med12 mutant
# adjust maximal RAM and load libraries
options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggh4x)
library(matrixStats)
library(ggpubr)
library(rstatix)

# Define Path with all Samples
path_list <- list.dirs(path = "./CellRanger_Matrix", recursive = FALSE)

# Load all samples
for (i in 1:length(path_list)){
  # Read in CellRanger Output
  temp <- Read10X(path_list[i])
  temp <- CreateSeuratObject(counts=temp$`Gene Expression`, project=sub(".*/", "", path_list[i]))
  # Filter based on mt percent
  temp$percent.MT <- PercentageFeatureSet(temp, pattern="^mt-" )
  temp <- subset(temp, subset = nFeature_RNA > 2500 & percent.MT < 15)
  # Assign Genotype
  if (grepl("M12KO", path_list[i], fixed=TRUE)==TRUE){
    temp$Med12_Genotype <- "Med12 KO"} else {
      temp$Med12_Genotype <- "wt"
    }
  # Assign Media condition
  if (grepl("Dif", path_list[i], fixed=TRUE)==TRUE){
    temp$Medium <- "Dif"} else {
      temp$Medium <- "2iL"
    }
  # Return as individual seurat obejcts
  assign(sub(".*/", "", path_list[i]), temp)
}


# Integrate Differentiated Samples for cluster assignment
data.list <- list(iG4_Dif,iG6_Dif, iG6_M12KO_A6_Dif, iG6_M12KO_A11_Dif, iG6_M12KO_A12_Dif)

data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, 
                                  reduction = "rpca")
# this command creates an 'integrated' data assay
Dif_Integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(Dif_Integrated) <- "integrated"

Dif_Integrated <- ScaleData(Dif_Integrated, verbose = FALSE)
Dif_Integrated <- RunPCA(Dif_Integrated, features = VariableFeatures(object = Dif_Integrated))
Dif_Integrated <- FindNeighbors(Dif_Integrated, dims = 1:12)
Dif_Integrated <- FindClusters(Dif_Integrated, resolution = 0.05, algorithm = 1)
Dif_Integrated <- RenameIdents(Dif_Integrated, "0" = "PrE", "1" = "Epi", "2" = "undefined")
Dif_Integrated$Celltype <- Idents(Dif_Integrated)


# Go back to all samples (Dif and Pluri)
All_samples <- merge(x = iG4_2iL,
                     y = c(iG4_Dif, iG6_2iL, iG6_Dif, iG6_M12KO_A11_2iL, iG6_M12KO_A11_Dif, iG6_M12KO_A12_2iL, iG6_M12KO_A12_Dif, iG6_M12KO_A6_2iL, iG6_M12KO_A6_Dif))

# Transfer clustering information from integrated Dif analysis to All_samples
# all cells not in Differentiated samples called Pluripotent
Pluri_cells <- colnames(All_samples)[which(!(colnames(All_samples) %in% colnames(Dif_Integrated)))]
Celltype_Meta <- rbind(data.frame(Dif_Integrated$Celltype), 
                       data.frame(row.names = Pluri_cells, 
                                  Dif_Integrated.Celltype = rep("Pluripotent", length(Pluri_cells))))
All_samples <- AddMetaData(All_samples, Celltype_Meta, col.name = "Celltype")


All_samples <- NormalizeData(All_samples, normalization.method = "LogNormalize", scale.factor = 10000)
All_samples <- FindVariableFeatures(All_samples, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(All_samples)
All_samples <- ScaleData(All_samples, features = all.genes)

# Find lists of DEGs between wt and Med12 KO
# based on wt (used for Figure 4K)
DefaultAssay(All_samples) <- "RNA"
All_samples$Celltype_Geno_Medium <- paste(All_samples$Celltype, All_samples$Med12_Genotype, All_samples$Medium)
Idents(All_samples) <- All_samples$Celltype_Geno_Medium

marker.Epi.wt <- FindMarkers(All_samples, ident.1 = "Epi wt Dif", ident.2 = "Pluripotent wt 2iL",logfc.threshold = 0.5)
marker.PrE.wt <- FindMarkers(All_samples, ident.1 = "PrE wt Dif", ident.2 = "Pluripotent wt 2iL",logfc.threshold = 0.5)

Epi_up <- rownames(subset(marker.Epi.wt, avg_log2FC > 0))
Epi_down <- rownames(subset(marker.Epi.wt, avg_log2FC < 0))
PrE_up <- rownames(subset(marker.PrE.wt, avg_log2FC > 0))
PrE_down <- rownames(subset(marker.PrE.wt, avg_log2FC < 0))


marker_expression <- data.frame(log2(AverageExpression(All_samples, features = c(PrE_up, PrE_down), group.by = "Celltype_Geno_Medium", slot = "data")$RNA+1))

marker_expression <- rbind(marker_expression, data.frame(log2(AverageExpression(All_samples, features = c(Epi_up, Epi_down), group.by = "Celltype_Geno_Medium", slot = "data")$RNA+1)))

marker_expression$Up_Down <- c(rep("up", length(PrE_up)), rep("down", length(PrE_down)), rep("up", length(Epi_up)), rep("down", length(Epi_down)))
marker_expression$Celltype <- c(rep("PrE", length(PrE_up)+length(PrE_down)), rep("Epi", length(Epi_up)+length(Epi_down)))

marker_expression$FC_wt_Epi <- marker_expression$Epi.wt.Dif - marker_expression$Pluripotent.wt.2iL
marker_expression$FC_mutant_Epi <- marker_expression$Epi.Med12.KO.Dif - marker_expression$Pluripotent.Med12.KO.2iL
marker_expression$FC_wt_PrE <- marker_expression$PrE.wt.Dif - marker_expression$Pluripotent.wt.2iL
marker_expression$FC_mutant_PrE <- marker_expression$PrE.Med12.KO.Dif - marker_expression$Pluripotent.Med12.KO.2iL

marker_expression$FC_wt <- ifelse(marker_expression$Celltype == "Epi", yes = marker_expression$FC_wt_Epi, no = marker_expression$FC_wt_PrE)
marker_expression$FC_mutant <- ifelse(marker_expression$Celltype == "Epi", yes = marker_expression$FC_mutant_Epi, no = marker_expression$FC_mutant_PrE)

marker_expression$Gene <- rownames(marker_expression)

marker_expression_long <- melt(marker_expression[,c(9:10,15:17)], )

stat.test <- marker_expression_long %>%
  group_by(Up_Down, Celltype) %>%
  wilcox_test(value ~ variable, paired = TRUE) %>%
  # t_test(FC_wt, FC_mutant) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

stat.test$y.position = c(0.5,0.5,3,3)

ggplot(marker_expression_long, aes(x=variable, y=value)) +
  geom_violin(aes(fill = variable), linewidth = 0.2) +
  geom_jitter(size = 0.2, width = 0.3, alpha = 0.7, pch = 16)+
  facet_grid(Up_Down~Celltype, scales = "free_y") +
  stat_pvalue_manual(stat.test, bracket.size = 0.46, tip.length = 0)+
  xlab("Genotype") +
  ylab("Expression Foldchange (log2)") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("#2E6A67", "#E98624")) +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) + 
  force_panelsizes(rows = unit(2.2, "cm"),cols = unit(2.2, "cm"))
ggsave("./Plots/Markers_during_Dif/20231124_wt_Differentiation_markers_Foldchanges_in_mutant_and_wt_logfc.threshold = 0.5.pdf", width = 4, height = 4)




# Compare DEGs based on mutant
DefaultAssay(All_samples) <- "RNA"
All_samples$Celltype_Geno_Medium <- paste(All_samples$Celltype, All_samples$Med12_Genotype, All_samples$Medium)
Idents(All_samples) <- All_samples$Celltype_Geno_Medium

marker.Epi.mutant <- FindMarkers(All_samples, ident.1 = "Epi Med12 KO Dif", ident.2 = "Pluripotent Med12 KO 2iL", logfc.threshold = 0.5)
marker.PrE.mutant <- FindMarkers(All_samples, ident.1 = "PrE Med12 KO Dif", ident.2 = "Pluripotent Med12 KO 2iL", logfc.threshold = 0.5)

PrE_up <- rownames(subset(marker.PrE.mutant, avg_log2FC > 0))
PrE_down <- rownames(subset(marker.PrE.mutant, avg_log2FC < 0))
Epi_up <- rownames(subset(marker.Epi.mutant, avg_log2FC > 0))
Epi_down <- rownames(subset(marker.Epi.mutant, avg_log2FC < 0))


marker_expression <- data.frame(log2(AverageExpression(All_samples, features = c(PrE_up, PrE_down), group.by = "Celltype_Geno_Medium", slot = "data")$RNA+1))

marker_expression <- rbind(marker_expression, data.frame(log2(AverageExpression(All_samples, features = c(Epi_up, Epi_down), group.by = "Celltype_Geno_Medium", slot = "data")$RNA+1)))

marker_expression$Up_Down <- c(rep("up", length(PrE_up)), rep("down", length(PrE_down)), rep("up", length(Epi_up)), rep("down", length(Epi_down)))
marker_expression$Celltype <- c(rep("PrE", length(PrE_up)+length(PrE_down)), rep("Epi", length(Epi_up)+length(Epi_down)))

marker_expression$FC_wt_Epi <- marker_expression$Epi.wt.Dif - marker_expression$Pluripotent.wt.2iL
marker_expression$FC_mutant_Epi <- marker_expression$Epi.Med12.KO.Dif - marker_expression$Pluripotent.Med12.KO.2iL
marker_expression$FC_wt_PrE <- marker_expression$PrE.wt.Dif - marker_expression$Pluripotent.wt.2iL
marker_expression$FC_mutant_PrE <- marker_expression$PrE.Med12.KO.Dif - marker_expression$Pluripotent.Med12.KO.2iL

marker_expression$FC_wt <- ifelse(marker_expression$Celltype == "Epi", yes = marker_expression$FC_wt_Epi, no = marker_expression$FC_wt_PrE)
marker_expression$FC_mutant <- ifelse(marker_expression$Celltype == "Epi", yes = marker_expression$FC_mutant_Epi, no = marker_expression$FC_mutant_PrE)

marker_expression$Gene <- rownames(marker_expression)

marker_expression_long <- melt(marker_expression[,c(9:10,15:17)], )

stat.test <- marker_expression_long %>%
  group_by(Up_Down, Celltype) %>%
  wilcox_test(value ~ variable, paired = TRUE) %>%
  # t_test(FC_wt, FC_mutant) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

stat.test$y.position = c(0.5,0.5,3,3)

ggplot(marker_expression_long, aes(x=variable, y=value)) +
  geom_violin(aes(fill = variable), linewidth = 0.2) +
  geom_jitter(size = 0.2, width = 0.3, alpha = 0.7, pch = 16)+
  facet_grid(Up_Down~Celltype, scales = "free_y") +
  stat_pvalue_manual(stat.test, bracket.size = 0.46, tip.length = 0)+
  xlab("Genotype") +
  ylab("Expression Foldchange (log2)") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("#2E6A67", "#E98624")) +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) + 
  force_panelsizes(rows = unit(2.2, "cm"),cols = unit(2.2, "cm"))
ggsave("./Plots/Markers_during_Dif/20231124_Med12Mutant_Differentiation_markers_Foldchanges_in_mutant_and_wt_logfc.threshold = 0.5.pdf", width = 4, height = 4)

