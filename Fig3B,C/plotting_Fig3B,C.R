# Create Figure panels Fig3B and C
# load libraries
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(ggh4x)
library(ggdendro)
library(eulerr)

# Load TPMs for PCA (exported with Seqmonk, normalized and quantitated as described in methods)
cts <- read.csv(file = "Annotated Probe Report for All Probes_Log2TPMs.txt",sep="\t")
# Make Probe names unique by including Positional information
rownames(cts) <- paste(cts$Probe, cts$Chromosome, cts$Start, cts$End, sep = "_") 
cts <- as.matrix(cts[,13:24])

# Calculate PCA and plot
PCAs <- prcomp(t(cts))
# print contributions of PCAs
PCAs$sdev^2 / sum(PCAs$sdev^2)

PCAs <- data.frame(PCAs$x)
PCAs$Replicate <- rep(c("1","2","3"), 4)
PCAs$Condition <- c(rep("wt in 2i", 3), rep("wt in N2B27", 3), 
                    rep("Med12 KO in 2i", 3), rep("Med12 KO in N2B27", 3))

# Plot PCA per sample
ggplot(PCAs, aes(x=PC1, y=PC2, color=Condition, shape=Replicate)) +
  geom_point(size=2) +
  xlab("PC1 (32.1%)") +
  ylab("PC2 (16.8%)") +
  scale_color_manual(values = c( "#E98624", "#F7D1A5","#2E6A67", "#91E2DA"))+
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_blank()) +
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(4, "cm"))
ggsave("./20231025_PCA_Plot_TPMs.pdf", width = 6, height = 5)

# For Differential gene expression analysis DESeq will be used, which takes unnormalized raw counts of reads as an input
# Load Raw counts
cts <- read.csv(file = "Annotated Probe Report for All Probes_Raw_Counts.txt",sep="\t")
# Make Probe names unique by including Positional information
rownames(cts) <- paste(cts$Probe, cts$Chromosome, cts$Start, cts$End, sep = "_") 
cts <- as.matrix(cts[,13:24])

# Define column with all metadata information (FGF state and Med12 genotype)
FGF_state = rep(c(rep("_2i", 3),rep("N2", 3)),2)
Med12_state = c(rep("wt", 6), rep("Med12_KO", 6))
coldata <- data.frame(State = paste(FGF_state, Med12_state, sep = "_"),
                      row.names = colnames(cts))

# Combine Counts and Metadata to DESeq object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ State)

# Filter lowly expressed genes (Raw count >= 100 in at least one sample)
keep <- rowSums(counts(dds) >= 100) >= 1
dds <- dds[keep,]

# Run Deseq2 and 
dds <- DESeq(dds)

# extract results as a list of all combinations+ applying log foldchange shrinkage
results <- list()
for (i in 1:6){
  cur_name <- paste(combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][1], 
                    combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][2], sep = "_vs_")
  results[cur_name] <- lfcShrink(dds, contrast = c("State", combn(unique(coldata$State)[c(2,4,1,3)],2)[,i]), type="ashr")
}

# Produce volcano plots
Comparison <- rbind(data.frame("Gene" = rownames(results$N2_wt_vs_N2_Med12_KO), 
                               "Comparison" = "wt vs. Med12 KO in N2B27",
                               "FC" = results$N2_wt_vs_N2_Med12_KO$log2FoldChange,
                               "p" = results$N2_wt_vs_N2_Med12_KO$padj),
                    data.frame("Gene" = rownames(results$N2_wt_vs_N2_Med12_KO),
                               "Comparison" = "wt vs. Med12 KO in 2i",
                               "FC" = results$`_2i_wt_vs__2i_Med12_KO`$log2FoldChange,
                               "p" = results$`_2i_wt_vs__2i_Med12_KO`$padj),
                    data.frame("Gene" = rownames(results$N2_wt_vs__2i_wt),
                               "Comparison" = "2i vs N2B27 in wt",
                               "FC" = results$N2_wt_vs__2i_wt$log2FoldChange,
                               "p" = results$N2_wt_vs__2i_wt$padj),
                    data.frame("Gene" = rownames(results$N2_Med12_KO_vs__2i_Med12_KO),
                               "Comparison" = "2i vs N2B27 in Med12 KO",
                               "FC" = results$N2_Med12_KO_vs__2i_Med12_KO$log2FoldChange,
                               "p" = results$N2_Med12_KO_vs__2i_Med12_KO$padj))
Comparison$Gene <- sub("_.*", "", Comparison$Gene) 
Comparison$Comparison <- factor(Comparison$Comparison, 
                                levels = c("wt vs. Med12 KO in 2i", "wt vs. Med12 KO in N2B27",
                                           "2i vs N2B27 in wt","2i vs N2B27 in Med12 KO"),
                                ordered = TRUE)

table(subset(Comparison, FC>1 & p<0.01)$Comparison)
table(subset(Comparison, FC<(-1) & p<0.01)$Comparison)

Comparison$signif <- ifelse((Comparison$FC<(-1) | Comparison$FC>1) & Comparison$p<0.01, 
                            "significant", "not significant")

DEG_2i_vs_N2B27_Med12KO <- subset(Comparison, (FC<(-1) | FC>1) & p<0.01 & Comparison == "2i vs N2B27 in Med12 KO")
DEG_2i_vs_N2B27_wt <- subset(Comparison, (FC<(-1) | FC>1) & p<0.01 & Comparison == "2i vs N2B27 in wt")
DEG_wt_vs_Med12KO_N2B27 <- subset(Comparison, (FC<(-1) | FC>1) & p<0.01 & Comparison == "wt vs. Med12 KO in N2B27")
DEG_wt_vs_Med12KO_2i <- subset(Comparison, (FC<(-1) | FC>1) & p<0.01 & Comparison == "wt vs. Med12 KO in 2i")

fill_colors <- c("#E98624", "#2E6A67","#A09A8D")

# Create plots
pdf("20231227_Euler_Plot_2i_vs_N2B27.pdf", width = 6, height = 6)
plot(euler(c("Med12 KO" = length(DEG_2i_vs_N2B27_Med12KO$Gene)-sum(DEG_2i_vs_N2B27_Med12KO$Gene %in% DEG_2i_vs_N2B27_wt$Gene),
      "wt" = length(DEG_2i_vs_N2B27_wt$Gene)-sum(DEG_2i_vs_N2B27_Med12KO$Gene %in% DEG_2i_vs_N2B27_wt$Gene),
      "Med12 KO&wt" = sum(DEG_2i_vs_N2B27_Med12KO$Gene %in% DEG_2i_vs_N2B27_wt$Gene))), 
     quantities = TRUE,
     fill = fill_colors)
dev.off()


fill_colors <- c("#D4D6D4","#A09A8D", "#BAB8B1")

pdf("20231227_Euler_Plot_wt_vs_Med12KO.pdf", width = 6, height = 6)
plot(euler(c("N2B27" = length(DEG_wt_vs_Med12KO_N2B27$Gene)-sum(DEG_wt_vs_Med12KO_N2B27$Gene %in% DEG_wt_vs_Med12KO_2i$Gene),
             "2i" = length(DEG_wt_vs_Med12KO_2i$Gene)-sum(DEG_wt_vs_Med12KO_N2B27$Gene %in% DEG_wt_vs_Med12KO_2i$Gene),
             "N2B27&2i" = sum(DEG_wt_vs_Med12KO_N2B27$Gene %in% DEG_wt_vs_Med12KO_2i$Gene))), 
     quantities = TRUE,
     fill = fill_colors)
dev.off()
