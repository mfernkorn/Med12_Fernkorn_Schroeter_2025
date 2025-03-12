# Plot Figure S7C
# Set RAM and load libraries

options(future.globals.maxSize = 8000 * 1024^2)

library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggh4x)
library(matrixStats)
library(viridis)

# Define Path with all Samples (output from cellranger)
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


# Integrate Differentiated Samples for proper cluster assignment
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
DefaultAssay(All_samples) <- "RNA"
All_samples$Celltype_Geno_Medium <- paste(All_samples$Celltype, All_samples$Med12_Genotype, All_samples$Medium)
Idents(All_samples) <- All_samples$Celltype_Geno_Medium

marker.pluri <- FindMarkers(All_samples, ident.1 = "Pluripotent wt 2iL", ident.2 = "Pluripotent Med12 KO 2iL")
marker.Epi <- FindMarkers(All_samples, ident.1 = "Epi wt Dif", ident.2 = "Epi Med12 KO Dif")
marker.PrE <- FindMarkers(All_samples, ident.1 = "PrE wt Dif", ident.2 = "PrE Med12 KO Dif")


# Plot bulk expression of hitlists in heatmap for Pluri, Epi and PrE
# Extract Top and bottum 10 genes by foldchange of every comparison
marker.pluri_top10 <- rownames(marker.pluri[order(marker.pluri$avg_log2FC, decreasing = TRUE)[c(1:10, (length(marker.pluri$avg_log2FC)-9):(length(marker.pluri$avg_log2FC)))],])
marker.Epi_top10 <- rownames(marker.Epi[order(marker.Epi$avg_log2FC, decreasing = TRUE)[c(1:10, (length(marker.Epi$avg_log2FC)-9):(length(marker.Epi$avg_log2FC)))],])
marker.PrE_top10 <- rownames(marker.PrE[order(marker.PrE$avg_log2FC, decreasing = TRUE)[c(1:10, (length(marker.PrE$avg_log2FC)-9):(length(marker.PrE$avg_log2FC)))],])

# Merge Gene lists of all Comparisons and information for color coding of tile border
marker <- merge(data.frame(Gene = marker.pluri_top10, Pluripotent.wt.2iL = "Pluri",Pluripotent.Med12.KO.2iL = "Pluri"), 
                data.frame(Gene = marker.Epi_top10, Epi.wt.Dif = "Epi",	Epi.Med12.KO.Dif = "Epi"), all=TRUE)
marker <- merge(marker, data.frame(Gene = marker.PrE_top10, PrE.wt.Dif = "PrE", PrE.Med12.KO.Dif = "PrE"), all=TRUE)
marker <- replace(marker, is.na(marker), "none")

# Get average expression data of define list
marker_expression <- data.frame(AverageExpression(All_samples, features =  marker$Gene)$RNA[,c(1,5,3,6,2,7)])
marker_expression <- as.data.frame(t(scale(t(marker_expression))))
marker_expression$Up_Down <- rowSums(marker_expression[,c(1,3,5)]) > rowSums(marker_expression[,c(2,4,6)])
marker_expression$Up_Down <- factor(marker_expression$Up_Down, levels = c(TRUE, FALSE), ordered = TRUE)
levels(marker_expression$Up_Down) <- c("Downregulated in Med12 mutant", "Upregulated in Med12 mutant")

marker_expression$Gene <- rownames(marker_expression)

# Restructure (melt) for ggplot
marker_expression_long <- melt(marker_expression)

# Add meta data (Up_Down) to marker list and restructure for ggplot
marker$Up_Down <- marker_expression$Up_Down
marker_long <- melt(marker, id.vars = c("Gene", "Up_Down"))

# Define meta data (Genotype and Celltype) and merge with data
meta_data <- data.frame(variable = levels(factor(marker_expression_long$variable)),
                        Celltype = c("Pluripotency", "Pluripotency","Epiblast", "Epiblast", "PrE", "PrE"),
                        Genotype = c("wt", "mutant","wt", "mutant", "wt", "mutant"))

marker_expression_long <- merge(marker_expression_long, meta_data)
marker_long <- merge(marker_long, meta_data)

# Define order of Genes
gene_order <- unique(c(rownames(marker.pluri[order(marker.pluri$avg_log2FC, decreasing = TRUE)[c(1:10)],]),
                     rownames(marker.Epi[order(marker.Epi$avg_log2FC, decreasing = TRUE)[c(1:10)],]),
                     rownames(marker.PrE[order(marker.PrE$avg_log2FC, decreasing = TRUE)[c(1:10)],]),
                     rownames(marker.pluri[order(marker.pluri$avg_log2FC, decreasing = TRUE)[(length(marker.pluri$avg_log2FC)-9):(length(marker.pluri$avg_log2FC))],]),
                     rownames(marker.Epi[order(marker.Epi$avg_log2FC, decreasing = TRUE)[(length(marker.Epi$avg_log2FC)-9):(length(marker.Epi$avg_log2FC))],]),
                     rownames(marker.PrE[order(marker.PrE$avg_log2FC, decreasing = TRUE)[(length(marker.PrE$avg_log2FC)-9):(length(marker.PrE$avg_log2FC))],])))
marker_expression_long$Gene <- factor(marker_expression_long$Gene, levels = rev(gene_order), ordered = TRUE)               
marker_long$Gene <- factor(marker_long$Gene, levels = rev(gene_order), ordered = TRUE)                
marker_expression_long$Celltype <- factor(marker_expression_long$Celltype, levels = c("Pluripotency", "Epiblast","PrE"), ordered = TRUE)            
marker_long$Celltype <- factor(marker_long$Celltype, levels = c("Pluripotency", "Epiblast","PrE"), ordered = TRUE)                
marker_expression_long$Genotype <- factor(marker_expression_long$Genotype, levels = c("wt", "mutant"), ordered = TRUE)               
marker_long$Genotype <- factor(marker_long$Genotype, levels = c("wt", "mutant"), ordered = TRUE)      

# Plotting
ggplot(marker_expression_long, aes(x=Genotype, y=Gene)) +
  geom_tile(aes(fill = value),height = 1, width = 1) +
  geom_tile(data = marker_long, aes(color = value), fill= NA, linewidth = 1,height = 1, width = 1) +
  scale_color_manual(values = c("#00A76E", "NA", "lightblue","#D86FAA")) +
  scale_fill_viridis(option = "cividis") +
  xlab("Med12 Genotype") + 
  facet_grid(cols = vars(Celltype),rows = vars(Up_Down), scales = "free_y", space = "free")+
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"),
        axis.line = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.text = element_text(color = "black"))
ggsave("./20231117_DEGs_Med12mutant_vs_wt_by_Celltype_Heatmap.pdf", width = 6, height = 6)

