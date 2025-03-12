# Export Flowjo data by:
# selecting desired group of cells in upper panel of Flowjo, which should be exported
# File->Export/Concatinate->Export Group
# Export with default parameters to csv scale values
# Put all files for one cdf plot in one folder
setwd("/Users/fernkorn/Documents/FACS_Analysis/KO_short_list_high/")
library(ggplot2)
library(ggh4x) # for force_panelsize option

# Load Replicate 1
filedir <- "./20230216_Rep1/Single_Cell_Data/"
file_names <- dir(filedir)
names_vec <- c("Flcn",  "Lamtor4", "Lztr1","wt + mock","Smarcc1","Tsc1", "wt")
file_directories <- paste0(filedir, file_names)

cdf_data_rep1 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep1 <- rbind(cdf_data_rep1, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep1) <- c("Fluorescence","Condition")
cdf_data_rep1$Condition <- factor(cdf_data_rep1$Condition, 
                                  levels = c("wt", "wt + mock","Flcn",  "Lamtor4","Tsc1","Lztr1","Smarcc1" ), 
                                  ordered = TRUE)

# Load Replicate 2
filedir <- "./20230223_Rep2/Single_Cell_Data/"
file_names <- dir(filedir)
names_vec <- c("Flcn",  "Lamtor4", "Lztr1","wt + mock","Smarcc1","Tsc1", "wt")
file_directories <- paste0(filedir, file_names)

cdf_data_rep2 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep2 <- rbind(cdf_data_rep2, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep2) <- c("Fluorescence","Condition")
cdf_data_rep2$Condition <- factor(cdf_data_rep2$Condition, 
                             levels = c("wt", "wt + mock","Flcn",  "Lamtor4","Tsc1","Lztr1","Smarcc1" ), 
                             ordered = TRUE)

# Load Replicate 3
filedir <- "./20230303_Rep3/Single_Cell_Data/"
file_names <- dir(filedir)
names_vec <- c("Flcn",  "Lamtor4", "Lztr1","wt + mock","Smarcc1","Tsc1","Tsc2", "wt")
file_directories <- paste0(filedir, file_names)

cdf_data_rep3 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep3 <- rbind(cdf_data_rep3, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep3) <- c("Fluorescence","Condition")
cdf_data_rep3$Condition <- factor(cdf_data_rep3$Condition, 
                                  levels = c("wt", "wt + mock","Flcn",  "Lamtor4","Tsc1","Tsc2", "Lztr1","Smarcc1"), 
                                  ordered = TRUE)
# Load Replicate 4
filedir <- "./20230313_Rep4/"
file_names <- dir(filedir)
names_vec <- c("wt", "Flcn",  "Lamtor4", "Lztr1","wt + mock","Smarcc1","Tsc1","Tsc2")
file_directories <- paste0(filedir, file_names)

cdf_data_rep4 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep4 <- rbind(cdf_data_rep4, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep4) <- c("Fluorescence","Condition")
cdf_data_rep4$Condition <- factor(cdf_data_rep4$Condition, 
                                  levels = c("wt", "wt + mock","Flcn",  "Lamtor4","Tsc1","Tsc2", "Lztr1","Smarcc1"), 
                                  ordered = TRUE)

# Load Replicate 5
filedir <- "./20230314_Rep5/"
file_names <- dir(filedir)
names_vec <- c("wt", "Flcn",  "Lamtor4", "Lztr1","wt + mock","Smarcc1","Tsc1","Tsc2")
file_directories <- paste0(filedir, file_names)

cdf_data_rep5 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep5 <- rbind(cdf_data_rep5, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep5) <- c("Fluorescence","Condition")
cdf_data_rep5$Condition <- factor(cdf_data_rep5$Condition, 
                                  levels = c("wt", "wt + mock","Flcn",  "Lamtor4","Tsc1","Tsc2", "Lztr1","Smarcc1"), 
                                  ordered = TRUE)


# Plot Medians and Mean of Replicates of both Replicates
# Combine Dataframes
library(dplyr)
cdf_data_rep1$Rep <- "1"
cdf_data_rep2$Rep <- "2"
cdf_data_rep3$Rep <- "3"
cdf_data_rep4$Rep <- "4"
cdf_data_rep5$Rep <- "5"
cdf_data <- rbind(cdf_data_rep1, cdf_data_rep2, cdf_data_rep3,cdf_data_rep4,cdf_data_rep5)

# calculate median per group
cdf_data <- group_by(cdf_data, Condition, Rep)
cdf_data_median <- summarise(cdf_data, "Median Fluorescence"=median(Fluorescence))

# PLot relative median for both replicates
cdf_data_median$relfluorescence <- NA
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 1)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt + mock" & cdf_data_median$Rep == 1)])[which(cdf_data_median$Rep == 1)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 2)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt + mock" & cdf_data_median$Rep == 2)])[which(cdf_data_median$Rep == 2)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 3)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt + mock" & cdf_data_median$Rep == 3)])[which(cdf_data_median$Rep == 3)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 4)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt + mock" & cdf_data_median$Rep == 4)])[which(cdf_data_median$Rep == 4)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 5)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt + mock" & cdf_data_median$Rep == 5)])[which(cdf_data_median$Rep == 5)]


cdf_data_median_sub <- subset(cdf_data_median, Condition != "wt" & Condition != "wt + mock")
cdf_data_median_sub$Condition <- factor(cdf_data_median_sub$Condition, 
                                        levels = c("Flcn",  "Lamtor4","Tsc1","Tsc2", "Lztr1","Smarcc1"), 
                                        ordered = TRUE)
ggplot(cdf_data_median_sub, aes(x=Condition, y=relfluorescence,  fill = Condition)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", color = "#144976", width = 0.7, size = 0.48) +
  geom_jitter(aes(shape = factor(Rep)), position = position_jitter(0.1), color = "#4C4C4C") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "#144976") +
  # stat_summary(aes(x = Condition), fun = mean, col = "darkred", shape = 3) +
  # geom_point(aes(shape = factor(Rep), color = "all")) + 
  scale_fill_manual(values = c('#BF6AAA','#BF6AAA','#BF6AAA','#BF6AAA','#87A9B9','#CBDB2A')) +
                                          scale_shape_manual(values = c(15,16,17,18,19)) +
  scale_y_continuous(expand = c(0,0,0,0.1)) + 
  ylab("Median Expression\nrelative to mock transfected wt")+
  xlab("Gene Knockout")+
  scale_x_discrete(guide = guide_axis(angle = 90))+
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
  force_panelsizes(rows = unit(3.5, "cm"),
                   cols = unit(6, "cm"))
ggsave("./20230504_Fig2E_bar_Medians_five_Reps_relative_to_mock.pdf", height = 12, width = 12, dpi = 1000, units = "cm")

# Statistics
# Alternative approach: Multiple testing corrected t test (paired and one sided)
# single t tests
data_for_statistics_wide <- reshape(data.frame(cdf_data_median[,1:3]), direction = "wide", idvar = "Rep", timevar = "Condition")
data_for_statistics_wide <- data_for_statistics_wide[,2:9]

# Assuming df is your data frame
# Assuming controlCondition is the name of your control condition column

controlCondition <- "Median.Fluorescence.wt + mock"  # replace this with the actual name of your control condition column

# Create a list to store the results
results_list <- list()

# Perform paired t-test for each condition against the control
for (condition in colnames(data_for_statistics_wide)) {
  if (condition != controlCondition) {
    # Perform paired t-test
    t_test_result <- t.test(data_for_statistics_wide[, condition], data_for_statistics_wide[, controlCondition], paired = TRUE, alternative = "greater")
    
    # Store the results in the list
    results_list[[condition]] <- t_test_result
  }
}

# Extract p-values
p_values <- sapply(results_list, function(x) x$p.value)

# Perform multiple testing correction (e.g., Benjamini-Hochberg)
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Display the results
results_df <- data.frame(
  Condition = names(results_list),
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

print(results_df)














# Plot with only mTOR signaling
cdf_data_median_sub <- subset(cdf_data_median, Condition != "wt" & Condition != "wt + mock" & Condition != "Lztr1" & Condition != "Smarcc1")
cdf_data_median_sub$Condition <- factor(cdf_data_median_sub$Condition, 
                                        levels = c("Flcn",  "Lamtor4","Tsc1","Tsc2"), 
                                        ordered = TRUE)
ggplot(cdf_data_median_sub, aes(x=Condition, y=relfluorescence)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", fill = "#9CBCD3", color = "#144976", width = 0.8) +
  geom_jitter(aes(shape = factor(Rep)), position = position_jitter(0.1), color = "#4C4C4C") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "#144976") +
  # stat_summary(aes(x = Condition), fun = mean, col = "darkred", shape = 3) +
  # geom_point(aes(shape = factor(Rep), color = "all")) + 
  # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fb9a99')) +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), axis.text = element_text(color = "black")) +
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,2.95)) + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  ylab("Median Expression\nrelative to mock transfected wt")+
  xlab("Gene Knockout") +
  force_panelsizes(rows = unit(2, "in"),
                   cols = unit(1.6, "in"))
ggsave("./20230324_bar_Medians_five_Reps_relative_to_mock_onlymTOR.pdf", width = 5, height = 4)




ggplot(cdf_data_median_sub, aes(x=Condition, y=relfluorescence, color = "all")) +
  # geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  stat_summary(aes(x = Condition), fun = mean, col = "darkred", shape = 3) +
  geom_point(aes(shape = factor(Rep))) + theme_bw(base_size = 12) +
  # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fb9a99')) +
  scale_color_manual(values = "#144976") +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  force_panelsizes(rows = unit(2, "in"),
                   cols = unit(1.6, "in")) +
  theme_bw(base_size = 12) +   theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0.5,0,0.5)) + scale_x_discrete(guide = guide_axis(angle = 90)) + 
  ylab("Median Fluorescence\nrelative to mock transfected wt")
# ggsave("./20230322_Medians_five_Reps_relative_to_mock_onlymTOR.pdf", width = 5, height = 4)

ggplot(cdf_data_median_sub, aes(x=Condition, y=relfluorescence, color = "all")) +
  # geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  stat_summary(aes(x = Condition), fun = mean, col = "darkred", shape = 3) +
  geom_point(aes(shape = factor(Rep))) + theme_bw(base_size = 12) +
  # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fb9a99')) +
  scale_color_manual(values = "#144976") +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  force_panelsizes(rows = unit(2, "in"),
                   cols = unit(2, "in")) +
  theme_bw(base_size = 12) +   theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0,0.5,0,0.5)) + scale_x_discrete(guide = guide_axis(angle = 90)) + 
  ylab("Median Fluorescence\nrelative to mock transfected wt")

