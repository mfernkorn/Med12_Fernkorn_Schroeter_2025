# Plot Figuer 5B left and right panel
# load libraries
library(ggplot2)
library("tidyr")
library(rstatix)
library(ggsignif)
library(ggh4x) # for force_panelsize option
library(reshape2)

# Load counts from all replicates
d_rep0_2iL <- read.csv("./20221013_2iL_counts.csv", sep = ";")
d_rep1_2iL <- read.csv("./20221122_2iL_counts.csv", sep = ";")
d_rep2_2iL <- read.csv("./20221128_2iL_counts.csv", sep = ";")
d_rep3_2iL <- read.csv("./20221130_2iL_counts.csv", sep = ";")
d_rep4_2iL <- read.csv("./20230213_2iL_counts_modified.csv", sep = ";")


d_rep0_2iL$Rep <- "0"
d_rep1_2iL$Rep <- "1"
d_rep2_2iL$Rep <- "2"
d_rep3_2iL$Rep <- "3"
d_rep4_2iL$Rep <- "4"

d_rep0_2iL$Colonies_relative <- d_rep0_2iL$Colonies/mean(subset(d_rep0_2iL, Condition == "Spry4Rep_2iL")$Colonies)
d_rep1_2iL$Colonies_relative <- d_rep1_2iL$Colonies/mean(subset(d_rep1_2iL, Condition == "Spry4Rep_2iL")$Colonies)
d_rep2_2iL$Colonies_relative <- d_rep2_2iL$Colonies/mean(subset(d_rep2_2iL, Condition == "Spry4Rep_2iL")$Colonies)
d_rep3_2iL$Colonies_relative <- d_rep3_2iL$Colonies/mean(subset(d_rep3_2iL, Condition == "Spry4Rep_2iL")$Colonies)
d_rep4_2iL$Colonies_relative <- d_rep4_2iL$Colonies/mean(subset(d_rep4_2iL, Condition == "Spry4Rep_2iL")$Colonies)

reps_2iL <- rbind(d_rep0_2iL, d_rep1_2iL, d_rep2_2iL, d_rep3_2iL, d_rep4_2iL)
write.csv2(reshape(reps_2iL, direction = "wide", timevar = "Rep", idvar = c("Condition", "Well")), file = "./20240105_Colony_numbers.csv")

# Average wells for statistical test
reps_2iL_Mean <- reshape2::dcast(reps_2iL, Rep + Condition ~ Well)
reps_2iL_Mean$Colonies_Mean <- rowMeans(reps_2iL_Mean[,c(3,4)])
reps_2iL_Mean$`1` <- NULL
reps_2iL_Mean$`2` <- NULL
# perform statistical tests
test <- kruskal_test(Colonies_Mean ~ Condition, data = reps_2iL_Mean)
test <- dunn_test(Colonies_Mean ~ Condition, data = reps_2iL_Mean, p.adjust.method = "bonferroni")

ggplot(reps_2iL, aes(x = Condition, y=Colonies_relative, color = factor(Well))) + 
  geom_point(aes(shape = Rep)) +
  scale_color_manual(values = c("black", "blue")) +
  scale_size_manual(values=c(2,2,2)) +
  scale_shape_manual(values=c(16,17,18,15,12)) +
  theme_bw() +theme(aspect.ratio=1, axis.text.x = element_text(angle = 90)) +
  stat_summary(aes(x = Condition), fun = mean, col = "darkred", shape = 3)
# ggsave("./20221122_2iL_Colony_numbers.pdf", width = 4, height = 4)

reps_2iL_sub_FGFKO <- subset(reps_2iL, Condition == "FGFKO" | Condition == "FGFKO+FGF" | Condition == "FGFKOMed12KO" | Condition == "FGFKOMed12KO+FGF")
reps_2iL_sub_FGFKO$Condition <- factor(reps_2iL_sub_FGFKO$Condition, levels =  c("FGFKO","FGFKOMed12KO", "FGFKO+FGF", "FGFKOMed12KO+FGF"), ordered = TRUE)

# Figure 5B right panel                                         
ggplot(reps_2iL_sub_FGFKO, aes(x = Condition, y=Colonies_relative)) + 
  geom_boxplot(aes(fill = Condition)) + geom_point(aes(shape = Rep, color = factor(Well)), size = 2) +
  scale_fill_manual(values = c("#2F6A67", "#E88624","#2F6A67", "#E88624")) +
  scale_color_manual(values = c("#164876", "#9CBCD3")) +
  scale_size_manual(values=c(2,2,2)) +
  scale_shape_manual(values=c(16,17,18,15,12)) +
  theme_bw(base_size = 12) +theme(aspect.ratio=1, ) +
  stat_signif(comparisons = list(c("FGFKO", "FGFKOMed12KO"), c("FGFKO+FGF", "FGFKOMed12KO+FGF"),c("FGFKO", "FGFKO+FGF"),c("FGFKOMed12KO", "FGFKOMed12KO+FGF")),
              map_signif_level = TRUE, y_position = c(1.4, 1.5, 1.6, 1.7)) +
  ylim(0,1.8) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90), 
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(1.8, "in"),
                   cols = unit(1.8, "in"))
ggsave("./20231128_2iL_FGFKO_Colony_numbers.pdf", width = 6, height = 6)

# Figure 5B left panel                                         
reps_2iL_sub_Spry4 <- subset(reps_2iL, Condition == "Spry4Rep" | Condition == "Spry4RepMed12KO")
ggplot(reps_2iL_sub_Spry4, aes(x = Condition, y=Colonies_relative)) + 
  geom_boxplot(aes(fill = Condition)) + geom_point(aes(shape = Rep, color = factor(Well)), size = 2) + 
  scale_color_manual(values = c("#164876", "#9CBCD3")) +
  scale_fill_manual(values = c("#2F6A67", "#E88624")) +
  scale_size_manual(values=c(2,2,2)) +
  scale_shape_manual(values=c(16,17,18,15,12)) +
  theme_bw(base_size = 12) + theme(aspect.ratio=1, axis.text.x = element_text(angle = 90)) +
  stat_signif(comparisons = list(c("Spry4Rep", "Spry4RepMed12KO")), map_signif_level = TRUE, y_position = c(1.4)) +
  ylim(0,1.8) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90), 
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(1.8, "in"),
                   cols = unit(1.08, "in"))
ggsave("./20231128_2iL_Spry4_Colony_numbers.pdf", width = 4, height = 6)


