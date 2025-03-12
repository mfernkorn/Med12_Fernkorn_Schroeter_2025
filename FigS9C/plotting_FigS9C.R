# Plot figure panel S9C
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

Pluri <- merge(x = iG4_2iL, 
               y = c(iG6_2iL, iG6_M12KO_A6_2iL, iG6_M12KO_A11_2iL, iG6_M12KO_A12_2iL))
Pluri$orig.ident <- factor(Pluri$orig.ident) # for heatmaps later
Pluri$Med12_Genotype <- factor(Pluri$Med12_Genotype, levels = c("wt","Med12 KO"), ordered = TRUE)

# Display numbers of cells per sample after filtering
table(Pluri@meta.data[["orig.ident"]])

# Normalize Data
Pluri <- NormalizeData(Pluri, normalization.method = "LogNormalize", scale.factor = 10000)
Pluri <- FindVariableFeatures(Pluri, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Pluri)
Pluri <- ScaleData(Pluri, features = all.genes)


Idents(Pluri) <-Pluri$Med12_Genotype
Ave_Exp <- AverageExpression(Pluri, features = c("Esrrb", "Klf4", "Nanog", "Prdm14", "Tbx3", "Tfcp2l1", "Zfp42"))$RNA # Aggregate Expression useless here, since different numbers of cells per category
ordered_genes <- rev(rownames(Ave_Exp[order(Ave_Exp[,1]/Ave_Exp[,2]),]))

Pluri <- JoinLayers(Pluri) # needed since Seurat v5
marker_exp <- melt(data.frame(t(as.matrix(Pluri@assays[["RNA"]]$data[c("Esrrb", "Klf4", "Nanog", "Prdm14", "Tbx3", "Tfcp2l1", "Zfp42"),])), 
                              Genotype = Pluri@meta.data[["Med12_Genotype"]]))
marker_exp$variable <- factor(marker_exp$variable, levels = ordered_genes, ordered = TRUE)

ggplot(marker_exp, aes(x=variable, y= value, fill = Genotype)) +
  geom_violin(position = position_dodge(width = 0.95), scale = "width", linewidth = 0.2)+
  scale_fill_manual(values = c("#2E6A67","#E98624"))+
  xlab("Gene") +
  ylab("Expression Level") +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black", angle = 45,hjust = 0.9),
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
  force_panelsizes(rows = unit(2.5, "cm"),
                   cols = unit(4, "cm"))
ggsave("./Plots/Figures/20231116_2iL_ViolinePlot_naive_markers.pdf", width = 10, height = 5)

