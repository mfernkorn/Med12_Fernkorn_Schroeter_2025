# Generate Figure panels Fig1D and Fig2B

# Load libraries
library(ggplot2)
library(reshape2)
library(cowplot)
library(pheatmap)

# Load data
Low <- read.csv2("./20220113_low_all_data_spry4_enriched_AND2-gate_fdr<0.25_OR_fdr<0.05.csv", sep = ";")
High <- read.csv2("./20220113_high_all_data_spry4_enriched_AND2-gate_fdr<0.25_OR_fdr<0.05.csv", sep = ";")

# Fig 1D (low Spry4 expression)
Low_FDR <- data.frame(Low$L11_fdr, Low$L51_fdr, Low$L12_fdr, Low$L52_fdr, row.names = Low$id)
colnames(Low_FDR) <- c("L_1%_S1", "L_5%_S1", "L_1%_S2", "L_5%_S2")
Low_FDR <- as.data.frame(t(Low_FDR))
pheatmap(Low_FDR, main = "FDR", scale = "none",
         col=rev(viridis::cividis(100)), display_numbers =  TRUE, fontsize_number = 10,number_color = "black",
         cellwidth = 22, cellheight = 22, cluster_rows = FALSE,
         filename = "/Users/fernkorn/Desktop/20220331_Heatmap_FDR_Low.pdf", width = 12, height = 4)

# Fig 2B (high Spry4 expression)
High_FDR <- data.frame(High$H11_fdr, High$H51_fdr, High$H12_fdr, High$H52_fdr, row.names = High$id)
colnames(High_FDR) <- c("H_1%_S1", "H_5%_S1", "H_1%_S2", "H_5%_S2")
High_FDR <- as.data.frame(t(High_FDR))
pheatmap(High_FDR, main = "FDR", scale = "none",
         col=rev(viridis::cividis(100)), display_numbers =  TRUE, fontsize_number = 10,number_color = "black",
         cellwidth = 22, cellheight = 22, cluster_rows = FALSE,
         filename = "/Users/fernkorn/Desktop/20220331_Heatmap_FDR_High.pdf", width = 12, height = 4)

