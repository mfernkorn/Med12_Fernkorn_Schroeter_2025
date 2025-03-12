# Plotting Figure 3H and I
# load libraries
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(ggh4x) # for force_panelsize option

# Load Raw counts
cts <- read.csv(file = "Annotated_Probe_Report_for_All_Probes.txt",sep="\t")
# Make Probe names unique by including Positional information
rownames(cts) <- paste(cts$Probe, cts$Chromosome, cts$Start, cts$End, sep = "_") 
cts <- as.matrix(cts[,13:24])

# Define column with all metadata information (FGF state and Med12 genotype)
FGF_state = rep(c(rep("unstimulated", 3),rep("stimulated", 3)),2)
Med12_state = c(rep("Med12_KO", 6), rep("Med12_wt", 6))
coldata <- data.frame(State = paste(FGF_state, Med12_state, sep = "_"),
                      row.names = colnames(cts))

# Combine Counts and Metadata to DESeq object
dds <- DESeqDataSetFromMatrix(countData = cts,
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

# Compare FGF target genes in Med12 wt and KO
# Extract all genes in both lists
Comparison <- rbind(data.frame("Gene" = rownames(results$stimulated_Med12_KO_vs_unstimulated_Med12_KO), 
                               "FC_Med12KOFGFvsMed12KO" = 2**results$stimulated_Med12_KO_vs_unstimulated_Med12_KO$log2FoldChange,
                               "FC_Med12wtFGFvsMed12wt" = 2**results$stimulated_Med12_wt_vs_unstimulated_Med12_wt$log2FoldChange,
                               "p_Med12KOFGFvsMed12KO" = results$stimulated_Med12_KO_vs_unstimulated_Med12_KO$padj,
                               "p_Med12wtFGFvsMed12wt" = results$stimulated_Med12_wt_vs_unstimulated_Med12_wt$padj))

Comparison$Gene <- sub("_.*", "", Comparison$Gene)
sum(Comparison$p_Med12KOFGFvsMed12KO < 0.01, na.rm = TRUE)
sum(Comparison$p_Med12wtFGFvsMed12wt < 0.01, na.rm = TRUE)

# Assign classes
Comparison$signific <- "non-sign"
Comparison$signific[Comparison$p_Med12KOFGFvsMed12KO<0.05] <- "sign_KO"
Comparison$signific[Comparison$p_Med12wtFGFvsMed12wt<0.05] <- "sign_wt"
Comparison$signific[Comparison$p_Med12wtFGFvsMed12wt<0.05 & Comparison$p_Med12KOFGFvsMed12KO<0.05] <- "sign_both"
Comparison$signific <- factor(Comparison$signific, levels = c("sign_both","sign_KO","sign_wt", "non-sign"), ordered = TRUE)

table(Comparison$signific)

# Plot for all target genes
plot <- scat_plot(dat = Comparison, geneselection = c("Spry4", Comparison$Gene[1:15]), include_density = FALSE, include_fit = FALSE, column_x = "FC_Med12wtFGFvsMed12wt", column_y = "FC_Med12KOFGFvsMed12KO",
          x_axis_lab = "Gene expression change\nupon FGF in wild-type", y_axis_lab = "Gene expression change\nupon FGF in Med12 mutant", log2_values = FALSE)
plot + geom_vline(xintercept = 1) + geom_hline(yintercept = 1)
ggsave(filename = "./20231228_FC_Comparison_FGFinMed12KOandwt_all_genes_labeled.pdf", width = 5, height = 5)
ggsave(filename = "./202300524_FC_Comparison_FGFinMed12KOandwt_all_genes.png", width = 5, height = 5)

# Subset and Plot for all upregulated genes, set threshhold dependent on preference
Comparison_pos <- subset(Comparison, Comparison$p_Med12wtFGFvsMed12wt<0.05 & Comparison$FC_Med12wtFGFvsMed12wt>1)

# Subset and Plot for all downregulated genes, set threshhold dependent on preference
Comparison_neg <- subset(Comparison, Comparison$p_Med12wtFGFvsMed12wt<0.05 & Comparison$FC_Med12wtFGFvsMed12wt<1)

# Of FGF target genes Plot Differences of Foldchanges in wt and Med12KO (with non log transformed foldchanges)
Comparison_pos$ratio <- Comparison_pos$FC_Med12KOFGFvsMed12KO/Comparison_pos$FC_Med12wtFGFvsMed12wt
Comparison_pos$direction <- "upregulated"
Comparison_neg$ratio <- Comparison_neg$FC_Med12KOFGFvsMed12KO/Comparison_neg$FC_Med12wtFGFvsMed12wt
Comparison_neg$direction <- "downregulated"

Comparison_target_genes <- rbind(Comparison_pos, Comparison_neg)

ggplot(Comparison_target_genes, aes(x=direction, y=ratio)) + 
  geom_violin(linewidth = 0.4) +
  geom_boxplot(width=0.2, outlier.size = 0.5) +
  ylab("Foldchange Ratio\nMed12 mutant vs wild type") +
  scale_y_continuous(trans = "log10") + 
  xlab("FGF target genes") +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(2.5, "cm"),
                   cols = unit(5, "cm")) +
  coord_flip()
ggsave("20231129_Comparison_FC_FGFtargetgenes_violinplot_p<0.01_non-log-transformed_coord_flip.pdf", width = 4, height = 5)


# Calculate Mode, Median and Mean of target gene foldchange ratio distributions
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Averages <- data.frame(Mode = c(Mode(round(Comparison_neg$ratio, digits = 3)), Mode(round(Comparison_pos$ratio, digits = 3))),
           Mean = c(mean(Comparison_neg$ratio),mean(Comparison_pos$ratio)),
           Median = c(median(Comparison_neg$ratio),median(Comparison_pos$ratio)))
rownames(Averages) <- c("Downregulated", "Upregulated")
print(Averages)
