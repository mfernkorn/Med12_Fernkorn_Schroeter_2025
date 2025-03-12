# Plot Figure 3H
# load library
library(ggplot2)
library(ggh4x) # for force_panelsize option
library(dplyr)
library(reshape2)
library(ggsignif)

# Load FACS data from different clones
filedir <- "./FACS_Data/"
file_names <- dir(filedir)
names_vec <- c("Med12 KO Cl2",  "Med12 KO Cl6","Med12 KO Cl7", "mock Cl2","mock Cl6","mock Cl7")
file_directories <- paste0(filedir, file_names)

cdf_data_rep1 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep1 <- rbind(cdf_data_rep1, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,17]), names_vec[i]))
}
colnames(cdf_data_rep1) <- c("Fluorescence","Condition")

# get summary statistics from FACS data per sample
summary_df <- cdf_data_rep1 %>%
  group_by(Condition) %>%
  summarize(median = median(Fluorescence), sd = sd(Fluorescence))
summary_df$Clone <- sub(".* ", "", summary_df$Condition)
summary_df$Gene_KO <- sub(" .*", "", summary_df$Condition)

summary_df$Gene_KO <- factor(summary_df$Gene_KO, 
                                  levels = c("mock","Med12"), 
                                  ordered = TRUE)

# Plot as point plot with path linking clones
ggplot(summary_df, aes(x= Gene_KO, y = median, group = Clone)) +
  geom_path(aes(color=Clone)) +
  # geom_errorbar(aes(ymin = median - sd, ymax = median + sd), width = 0.2, alpha = 0.7) +
  geom_point(aes(color=Clone)) +
  scale_color_manual(values = c("#C80000", "#329632", "#326496")) +
  xlab("Gene Knockout") + 
  ylab("Median Fluorescence") +
  ylim(3000, 17500)+
  scale_x_discrete(expand = c(0.1,0.11)) +
  geom_signif(comparisons = list(c("mock", "Med12")), test ="t.test",test.args = list(paired = TRUE), map_signif_level=TRUE, tip_length = 0, y_position = 15500) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.46),  # Adjust the size as needed
        text = element_text(color = 'black'),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_blank()) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(1, "cm"))
ggsave("./20231130_Induction_levels_median_ttest.pdf", width = 4, height = 5)