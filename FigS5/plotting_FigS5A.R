# Plot Fig S5A
#load libraries
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(stringr)
library(ggpubr)
library(rstatix)
library(ggh4x)

# Load TPMs (exported from Seqmonk)
cts <- read.csv(file = "/Users/fernkorn/Documents/RNA-seq_data/Own_Data/bulk_RNA-seq/Med12KO_24h_Comparison/R_Downstream_Analysis/Annotated Probe Report for All Probes_Log2TPMs.txt",sep="\t", dec = ",")
# Make Probe names unique by including Positional information
rownames(cts) <- paste(cts$Probe, cts$Chromosome, cts$Start, cts$End, sep = "_") 
cts <- data.frame(as.matrix(cts[,c(1,13:24)]))

Med12like <- data.frame(t(subset(cts, Probe == "Med12l")[,2:13]))
colnames(Med12like) <- "Med12l"
Med12like$Med12l <- 2**as.numeric(Med12like$Med12l)
Med12like$Rep <- as.numeric(str_sub(rownames(Med12like), -1,-1))
Med12like$Condition <- str_sub(rownames(Med12like), -10,-3)
Med12like$Condition <- factor(Med12like$Condition, levels = c("wt_2i","M12KO_2i","wt_N2","M12KO_N2"), ordered = TRUE)

ggplot(Med12like, aes(x=Condition, y=Med12l, fill = Condition)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  geom_point(aes(shape = factor(Rep))) +
  scale_fill_manual(values=c('#2E6A67','#E98624','#91E2DA','#F7D1A5')) +
  # theme_bw(base_size = 12) +   theme(legend.position = "none", panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0,0,15)) + 
  scale_x_discrete(guide = guide_axis(angle = 90), labels=c("M12KO_2i" = "Med12 mutant in 2i", "M12KO_N2" = "Med12 mutant in N2B27",
                                                            "wt_2i" = "wt in 2i", "wt_N2" = "wt in N2B27")) + 
  ylab("Med12L expression (TPM)")+
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
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(4, "cm"))
ggsave("./20231129_Med12L_Expression_TPM_FGFwt_exp.pdf", height =4, width = 3)



