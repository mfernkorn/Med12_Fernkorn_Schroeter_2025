# Plot Figure Panels S8A, 4H, 4I and 4J
# load libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(stats)
library(ggpubr)
library(reshape2)
library(dplyr)

# Define Path with all Samples (from Cellranger output)
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

# Cell number per sample
table(All_Samples@meta.data[["orig.ident"]])

# Compare number of RNAs between samples before normalization (Figure S8A)
meta_data <- All_Samples@meta.data
meta_data$Geno_Med <- paste(meta_data$Med12_Genotype, meta_data$Medium, sep = "_")
meta_data$Geno_Med <- factor(meta_data$Geno_Med, levels = c("wt_2iL", "Med12 KO_2iL", "wt_Dif", "Med12 KO_Dif"), ordered = TRUE)
# load  datasets for k-s-tesing
data_list <- list(subset(meta_data, Medium == "2iL" & Med12_Genotype == "wt")$nCount_RNA, 
                  subset(meta_data, Medium == "2iL" & Med12_Genotype == "Med12 KO")$nCount_RNA,
                  subset(meta_data, Medium == "Dif" & Med12_Genotype == "wt")$nCount_RNA,
                  subset(meta_data, Medium == "Dif" & Med12_Genotype == "Med12 KO")$nCount_RNA)

num_samples <- length(data_list)
p_values_matrix <- matrix(nrow = num_samples, ncol = num_samples)

# Perform KS tests for all sample pairs
for (i in 1:(num_samples - 1)) {
  for (j in (i + 1):num_samples) {
    ks_result <- ks.test(data_list[[i]], data_list[[j]], exact = TRUE)
    p_values_matrix[i, j] <- ks_result$p.value
    p_values_matrix[j, i] <- ks_result$p.value
  }
}

# Apply Bonferroni correction
adjusted_p_values <- p.adjust(as.vector(p_values_matrix), method = "bonferroni")

# Reshape the adjusted p-values matrix
adjusted_p_values_matrix <- data.frame(matrix(adjusted_p_values, nrow = num_samples))
colnames(adjusted_p_values_matrix) <- c("wt_2iL","Med12 KO_2iL", "wt_Dif", "Med12 KO_Dif")
rownames(adjusted_p_values_matrix) <- c("wt_2iL","Med12 KO_2iL", "wt_Dif", "Med12 KO_Dif")
adjusted_p_values_matrix <- data.frame(group1 = c("wt_2iL","wt_2iL","Med12 KO_2iL", "wt_Dif"), 
                                       group2 = c("Med12 KO_2iL", "wt_Dif", "Med12 KO_Dif", "Med12 KO_Dif"),
                                       p = c(adjusted_p_values_matrix[2,1], adjusted_p_values_matrix[3,1], 
                                             adjusted_p_values_matrix[2,4],adjusted_p_values_matrix[3,4]),
                                       y.position = c(5.2,5.3,5.4,5.5),
                                       signif_level = c("n.s.","***","***","***"))
# Print the adjusted p-values matrix
print(adjusted_p_values_matrix)

#Compare medians
meta_data %>% group_by(Geno_Med) %>% summarize(median = median(nCount_RNA))
meta_data %>% group_by(orig.ident) %>% summarize(median = median(nCount_RNA))
meta_data %>% group_by(orig.ident) %>% summarize(median = median(nFeature_RNA))
meta_data %>% group_by(orig.ident) %>% summarize(median = median(percent.MT))

# Plot violin plots (Figure S8A)
ggplot(meta_data, aes(x=Geno_Med, y= nCount_RNA)) +
  geom_violin(aes(fill=Geno_Med), position = position_dodge(width = 0.9), scale = "width",width = 0.85, linewidth = 0.2) +
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.2,outlier.size = 0.3,linewidth = 0.2) + scale_y_log10() +
  scale_fill_manual(values = c("#2E6A67","#E98624","#91E2DA","#F7D1A5"))+
  xlab("Condition") +
  ylab("mRNA Counts") +
  stat_pvalue_manual(adjusted_p_values_matrix, size = 4, label = "signif_level", bracket.size =0.46,tip.length = 0) +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text.x = element_text(colour = "black", angle = 45,hjust = 0.9),
        axis.text.y = element_text(colour = "black"),
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
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(2.5, "cm"))
ggsave("./Plots/Figures/20231108_Per_Genotype_and_Medium_Counts_per_cell.pdf", width = 4, height = 5)



# For Figure 4H to K integrate data for differentiated samples
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

table(Dif_Integrated$seurat_clusters)

Dif_Integrated@reductions[["pca"]]@stdev^2 / sum(Dif_Integrated@reductions[["pca"]]@stdev^2)*100
total_variance <- sum(Dif_Integrated@reductions[["pca"]]@stdev^2)
total_variance

PCA_data <- data.frame(Dif_Integrated@reductions[["pca"]]@cell.embeddings[,c(1,2)], Cluster = Dif_Integrated@meta.data[["seurat_clusters"]],
                       Sample = Dif_Integrated@meta.data[["orig.ident"]])

# Plot PCA colored by cluster (Figure 4H)
ggplot(PCA_data, aes(x=PC_1, y=PC_2, color = Cluster)) + 
  geom_point(pch = 16,size = 0.1)+ 
  xlab("PC 1 (15.7%)") +
  ylab("PC 2 (12.4%)") +
  facet_wrap(~ Sample, ncol = 3) +
  scale_color_manual(values = c("#D86FAA","#00A76E","black"))+
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
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill = "white", color = "white"),  # Remove the border
        strip.text = element_text(color = "black")) +
  force_panelsizes(rows = unit(2, "cm"),
                   cols = unit(2, "cm"))
ggsave("./Plots/Figures/20231116_Dif_Samples_Clustering_PCA_After_Integration.pdf", width = 8, height = 3)

# Violin Plots of key marker genes for Epi and PrE in both cluster (Figure 4I)
Dif_Integrated@assays[["RNA"]] <- JoinLayers(Dif_Integrated@assays[["RNA"]])
marker_exp <- melt(data.frame(Gata6 =  Dif_Integrated@assays[["RNA"]]$data["Gata6",],
                              Nanog =  Dif_Integrated@assays[["RNA"]]$data["Nanog",],
                              Cluster = Dif_Integrated@meta.data[["seurat_clusters"]]))
marker_exp <- subset(marker_exp, Cluster != "2")
ggplot(marker_exp, aes(x=variable, y= value, fill = Cluster)) +
  geom_violin(position = position_dodge(width = 0.95), scale = "width", linewidth = 0.2)+
  scale_fill_manual(values = c("#D86FAA","#00A76E"))+
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
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(2, "cm"))
ggsave("./Plots/Figures/20231117_Dif_ViolinePlot_Gata6_Nanog_by_cluster.pdf", width = 4, height = 3)

# Bar chart of relative cluster frequencies per sample in Differentiated samples (Figure 4J)
# Obtain dataframe with relative cluster frequencies and meta data information
heatmap_prop_pre = data.frame(Dif_Integrated@meta.data[["seurat_clusters"]],
                              Dif_Integrated@meta.data[["orig.ident"]])
heatmap_prop <- prop.table(table(heatmap_prop_pre), margin = 2)
heatmap_prop_long <- melt(heatmap_prop)
colnames(heatmap_prop_long) <- c("Cluster", "Sample", "Frequency")
heatmap_prop_long$Cluster <- factor(heatmap_prop_long$Cluster)
heatmap_prop_long$Genotype <- c(rep("wt", 6), rep("mutant", 9))
heatmap_prop_long$Rep <- c(rep("1", 3), rep("2", 3),rep("1", 3), rep("2", 3),rep("3", 3))
heatmap_prop_long$Genotype <- factor(heatmap_prop_long$Genotype, levels = c("wt", "mutant"), ordered = TRUE)

ggplot(heatmap_prop_long, aes(x=Genotype, y=Frequency, fill= Cluster, group = Rep)) + 
  geom_bar(stat = "identity", width = 1.2, lwd = 0.2,color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#D86FAA","#00A76E","#686868"))+
  facet_wrap(Genotype~ Rep,scales = "free_x", nrow = 1) +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.minor = element_blank(),
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
        strip.text = element_text(color = "black"),
        panel.spacing.x = unit(c(0,4,0,0), "pt")) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(0.4, "cm"))
ggsave("./Plots/Figures/20231117_Dif_Samples_Cluster_Frequencies_Barchart_After_Integration.pdf", width = 4, height = 4)


# Save a list of marker genes between Epi and PrE Clusters
Dif_Integrated.markers <- FindMarkers(Dif_Integrated, ident.1 = 0, ident.2 = 1)
write.csv2(Dif_Integrated.markers, "./Plots/Figures/202311009_PrE_Epi_Cluster_Marker_Genes.csv")
