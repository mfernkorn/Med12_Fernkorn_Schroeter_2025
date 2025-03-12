# Plotting Figure Panels 4F and 4G
# load libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(stats)
library(ggpubr)
library(reshape2)
library(dplyr)

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

# Merge Samples to one Seurat Object
All_Samples <- merge(x = iG4_2iL, 
                 y = c(iG4_Dif, iG6_2iL, iG6_Dif, iG6_M12KO_A11_2iL, iG6_M12KO_A11_Dif, iG6_M12KO_A12_2iL, iG6_M12KO_A12_Dif, iG6_M12KO_A6_2iL, iG6_M12KO_A6_Dif), 
                 add.cell.ids = sub(".*/", "", path_list), 
                 project = "PrE_Epi")
All_Samples$Med12_Genotype <- factor(All_Samples$Med12_Genotype, levels = c("wt","Med12 KO"), ordered = TRUE)

All_Samples <- NormalizeData(All_Samples, normalization.method = "LogNormalize", scale.factor = 10000)
All_Samples <- FindVariableFeatures(All_Samples, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(All_Samples)
All_Samples <- ScaleData(All_Samples, features = all.genes)
All_Samples <- RunPCA(All_Samples)
All_Samples <- RunUMAP(All_Samples, dims = 1:12)
DimPlot(All_Samples)

# Access UMAP data in Seurat object, including. orig.ident for sample info (Figure 4F)
UMAP_data <- data.frame(All_Samples@reductions[["umap"]]@cell.embeddings[,c(1,2)], Sample = All_Samples@meta.data[["orig.ident"]], Medium = All_Samples@meta.data[["Medium"]])
UMAP_data <- UMAP_data[sample(nrow(UMAP_data)), , drop = FALSE]
UMAP_data$Sample <- factor(UMAP_data$Sample, ordered = TRUE)

ggplot(UMAP_data, aes(x=umap_1, y=umap_2, color = Sample)) + 
  geom_point(size = 0.3)+
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_color_manual(values = c("#35514F","#489D9A","#AED6D3","#B6FFFD","#CE8B49","#FF9228","#DCCFBF","#FFD8AA","#764A32","#C3602B"))+
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.48),
        panel.border = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill = "white", color = "white"),  # Remove the border
        strip.text = element_text(color = "black")) +
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(4, "cm"))
ggsave("./Plots/Figures/20240111_All_Samples_UMAP.pdf", width = 6, height = 3)



# FeaturePLots  (Figure 4G)
FeaturePlot(All_Samples, "Sox17", slot = "data")
# All_Samples <- JoinLayers(All_Samples)

UMAP_data <- FetchData(object = All_Samples, vars = c("umap_1","umap_2","Sp5", "Nanog","Fgf4","Dnmt3l","Sox17","Dab2","Cubn"), layer = "data")
UMAP_data <- melt(UMAP_data, id.vars = c("umap_1", "umap_2"))

# Plot UMAP colored by cluster
ggplot(UMAP_data, aes(x=umap_1, y=umap_2, color = value)) + 
  geom_point(size = 0.02)+ 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_viridis_c(rescaler = function(x, to = c(0, 1), from = NULL) {
    ifelse(x<1.8, 
           scales::rescale(x,
                           to = to,
                           from = c(min(x, na.rm = TRUE), 1.8)),
           1)}) +
  theme_bw(base_size = 10) + 
  facet_wrap(~variable, nrow = 1) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.48),
        panel.border = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill = "white", color = "white"),  # Remove the border
        strip.text = element_text(color = "black")) +
  force_panelsizes(rows = unit(1.8, "cm"),
                   cols = unit(1.8, "cm"))
ggsave("./Plots/Figures/20240111_Marker_gene_expression.pdf", width = 10, height = 3)






