# Plot Figure 5C
# load libraries
library(DESeq2)
library(reshape2)
library(ggplot2)
library(pheatmap)

# Compare timing between Lackner et al data and own RNA seq results (from Figure 3)
# Load Data from Lackner et al 2021
Lackner_raw <- read.csv("GSE145653_counts_2hrTC.csv")[,c(2,3,6:22,24)]
Lackner_raw <- Lackner_raw[!duplicated(Lackner_raw$mgi),] # Remove duplicate Genes
Lackner_raw <- Lackner_raw[which(Lackner_raw$mgi != ""),] # Remove unnamed genes
rownames(Lackner_raw) <- Lackner_raw$mgi
Lackner_raw$mgi <- NULL

# Load own raw data
Med12comp <- read.csv(file = "Annotated Probe Report for All Probes_Raw_Counts.txt",sep="\t")
Med12comp <- Med12comp[!duplicated(Med12comp$Probe),] # Remove duplicate Genes
Med12comp_raw <- Med12comp[,c(1,13:24)]
rownames(Med12comp_raw) <- Med12comp_raw$Probe
Med12comp_raw$Probe <- NULL

# Deseq analysis for both data sets on their own
# Own data
# Define column with all metadata information (FGF state and Med12 genotype)
FGF_state = rep(c(rep("twoi", 3),rep("N2", 3)),2)
Med12_state = c(rep("wt", 6), rep("Med12_KO", 6))
coldata <- data.frame(State = paste(FGF_state, Med12_state, sep = "_"),
                      row.names = colnames(Med12comp_raw))

# Combine Counts and Metadata to DESeq object
dds <- DESeqDataSetFromMatrix(countData = Med12comp_raw,
                              colData = coldata,
                              design = ~ State)

# Filter lowly expressed genes (Raw count >= 100 in at least one sample)
keep <- rowSums(counts(dds) >= 100) >= 1
dds <- dds[keep,]

# Run Deseq2 and extract results as a list of all combinations+ applying log foldchange shrinkage
dds <- DESeq(dds)
results <- list()
for (i in 1:6){
  cur_name <- paste(combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][1], 
                    combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][2], sep = "_vs_")
  results[cur_name] <- lfcShrink(dds, contrast = c("State", combn(unique(coldata$State)[c(2,4,1,3)],2)[,i]), type="ashr")
}
Comparison_own <- rbind(data.frame("Gene" = rownames(results$N2_wt_vs_N2_Med12_KO), 
                               "FC_wt" = results$N2_wt_vs_twoi_wt$log2FoldChange,
                               "FC_Med12KO" = results$N2_Med12_KO_vs_twoi_Med12_KO$log2FoldChange))


# Lackner data
coldata <- data.frame(State = c("twoi","twoi","N2","N4","N6","N8","N10","N12","N14","N16",
                                "N18","N20","N22","N24","N26","N28","N30","N32","N32"),
                      row.names = colnames(Lackner_raw[,1:19]))
dds <- DESeqDataSetFromMatrix(countData = Lackner_raw[,1:19],
                              colData = coldata,
                              design = ~ State)

# Filter lowly expressed genes (Raw count >= 100 in at least one sample)
keep <- rowSums(counts(dds) >= 100) >= 1
dds <- dds[keep,]
dds <- DESeq(dds)


results <- list()
for (i in 2:17){
  cur_name <- paste(unique(coldata$State), "twoi", sep = "_vs_")
  results[cur_name[i]] <- lfcShrink(dds, contrast = c("State",unique(coldata$State)[i],"twoi"), type="ashr")
}

# Extract FC for all timepoints
Comparison_timecourse <- rbind(data.frame("Gene" = rownames(results$N2_vs_twoi), 
                               "FC_2h" = results$N2_vs_twoi$log2FoldChange,
                               "FC_4h" = results$N4_vs_twoi$log2FoldChange,
                               "FC_6h" = results$N6_vs_twoi$log2FoldChange,
                               "FC_8h" = results$N8_vs_twoi$log2FoldChange,
                               "FC_10h" = results$N10_vs_twoi$log2FoldChange,
                               "FC_12h" = results$N12_vs_twoi$log2FoldChange,
                               "FC_14h" = results$N14_vs_twoi$log2FoldChange,
                               "FC_16h" = results$N16_vs_twoi$log2FoldChange,
                               "FC_18h" = results$N18_vs_twoi$log2FoldChange,
                               "FC_20h" = results$N20_vs_twoi$log2FoldChange,
                               "FC_22h" = results$N22_vs_twoi$log2FoldChange,
                               "FC_24h" = results$N24_vs_twoi$log2FoldChange,
                               "FC_26h" = results$N26_vs_twoi$log2FoldChange,
                               "FC_28h" = results$N28_vs_twoi$log2FoldChange,
                               "FC_30h" = results$N30_vs_twoi$log2FoldChange,
                               "FC_32h" = results$N32_vs_twoi$log2FoldChange))

# Subset both dataset to just contain pluripotency marker genes
pluripotency_markers <- c("Sox2", "Klf4", "Nanog", "Esrrb", "Tbx3", "Tfcp2l1","Prdm14", "Zfp42") # Defined based on list in Lacker et al, unclear if the correct list, added Sox2
Differentiation_DEGs <- read.csv2("./DEGs_Lackner.csv", header = FALSE)

Comparison_timecourse <- subset(Comparison_timecourse, Gene %in% pluripotency_markers)
Comparison_timecourse <- subset(Comparison_timecourse, Gene %in% Differentiation_DEGs$V1)

Comparison_timecourse$FC_0h <- 0
Comparison_timecourse <- Comparison_timecourse[order(Comparison_timecourse$Gene),]
Comparison_own <- subset(Comparison_own, Gene %in% pluripotency_markers)
Comparison_own <- Comparison_own[order(Comparison_own$Gene),]

# Check foldchanges of Lackner time course data
Comparison_timecourse_plot <- melt(Comparison_timecourse)
Comparison_timecourse_plot$variable <- factor(Comparison_timecourse_plot$variable, 
                                          levels = c("FC_0h", "FC_2h","FC_4h","FC_6h","FC_8h","FC_10h","FC_12h","FC_14h","FC_16h","FC_18h","FC_20h","FC_22h","FC_24h","FC_26h","FC_28h","FC_30h","FC_32h"),
                                          ordered = TRUE)
ggplot(Comparison_timecourse_plot, aes(x = variable, y= value, group = Gene, color = Gene)) + geom_line() # does not like smoothing is very necessary

# Calculate euclidean distance to own data
euclid <- data.frame(time = c(seq(2,32,2),0), wt = c(1:17), Med12KO = c(1:17))
for (i in 2:18){
  euclid$wt[i-1] <- sqrt(sum((Comparison_timecourse[,i] - Comparison_own$FC_wt)^2))
  euclid$Med12KO[i-1] <- sqrt(sum((Comparison_timecourse[,i] - Comparison_own$FC_Med12KO)^2))
}
# Normalize euclidean distance
euclid$wt <- (euclid$wt-min(euclid$wt))/(max(euclid$wt)-min(euclid$wt))
euclid$Med12KO <- (euclid$Med12KO-min(euclid$Med12KO))/(max(euclid$Med12KO)-min(euclid$Med12KO))

# Plot Heatmap (Figure 5C)
euclid <- euclid[order(euclid$time, decreasing = TRUE),]
rownames(euclid) <- (euclid$time-24)*-1
euclid$time <- NULL
my.breaks <- seq(0,1,0.01)
my.colors <- c(colorRampPalette(colors = c("red", "black"))(length(my.breaks)/20), colorRampPalette(colors = c("black", "white"))(length(my.breaks)*0.95))
pheatmap(t(euclid), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = my.colors, border_color = FALSE, 
         cellwidth = 40, cellheight = 10)





       