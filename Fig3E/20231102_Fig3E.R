# Plot Figure 3H
# Load libraries
library(ggplot2)
library(ggh4x) # for force_panelsize option
library(dplyr)

# Load data from Med12 KO line with mTOR mutations
# Load Replicate 1
filedir <- "./20230306_Rep1/"
file_names <- dir(filedir)
names_vec <- c("Flcn",  "Lamtor4","Med12KO", "Med12KO + mock","Tsc1","Tsc2", "wt")
file_directories <- paste0(filedir, file_names)

cdf_data_rep1 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep1 <- rbind(cdf_data_rep1, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep1) <- c("Fluorescence","Condition")
cdf_data_rep1$Condition <- factor(cdf_data_rep1$Condition, 
                                  levels = c("wt", "Med12KO", "Med12KO + mock","Flcn","Lamtor4","Tsc1","Tsc2"), 
                                  ordered = TRUE)

# Load Replicate 2
filedir <- "./20230313_Rep2/"
file_names <- dir(filedir)
names_vec <- c("Med12KO","Flcn","Lamtor4","Med12KO + mock","Tsc1","Tsc2", "wt")
file_directories <- paste0(filedir, file_names)

cdf_data_rep2 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep2 <- rbind(cdf_data_rep2, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep2) <- c("Fluorescence","Condition")
cdf_data_rep2$Condition <- factor(cdf_data_rep2$Condition, 
                                  levels = c("wt", "Med12KO", "Med12KO + mock","Flcn","Lamtor4","Tsc1","Tsc2"), 
                                  ordered = TRUE)

# Load Replicate 3
filedir <- "./20230314_Rep3/"
file_names <- dir(filedir)
names_vec <- c("Med12KO","Flcn","Lamtor4","Lztr1", "Med12KO + mock","Smarcc1","Tsc1","Tsc2", "wt")
file_directories <- paste0(filedir, file_names)

cdf_data_rep3 <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data_rep3 <- rbind(cdf_data_rep3, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data_rep3) <- c("Fluorescence","Condition")
cdf_data_rep3$Condition <- factor(cdf_data_rep3$Condition, 
                                  levels = c("wt", "Med12KO", "Med12KO + mock","Flcn","Lamtor4","Tsc1","Tsc2","Lztr1","Smarcc1"), 
                                  ordered = TRUE)

# Group data for multiple replicates:
cdf_data_rep1$Rep <- "1"
cdf_data_rep2$Rep <- "2"
cdf_data_rep3$Rep <- "3"
cdf_data <- rbind(cdf_data_rep1, cdf_data_rep2, cdf_data_rep3)

# calculate median per group
cdf_data <- group_by(cdf_data, Condition, Rep)
cdf_data_median <- summarise(cdf_data, "Median Fluorescence"=median(Fluorescence))



# relative to wt
cdf_data_median$relfluorescence <- NA
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 1)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 1)])[which(cdf_data_median$Rep == 1)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 2)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 2)])[which(cdf_data_median$Rep == 2)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 3)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 3)])[which(cdf_data_median$Rep == 3)]

cdf_data_median_sub <- subset(cdf_data_median, Condition != "wt" & Condition != "Med12KO" & Condition != "Lztr1" & Condition != "Smarcc1")
cdf_data_median_sub$Condition <- factor(cdf_data_median_sub$Condition, 
                                        levels = c("Med12KO + mock", "Flcn",  "Lamtor4","Tsc1","Tsc2"), 
                                        labels = c("mock", "Flcn",  "Lamtor4","Tsc1","Tsc2"),
                                        ordered = TRUE)




# Load data from wt line same as for Fig2E 

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
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 1)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 1)])[which(cdf_data_median$Rep == 1)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 2)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 2)])[which(cdf_data_median$Rep == 2)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 3)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 3)])[which(cdf_data_median$Rep == 3)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 4)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 4)])[which(cdf_data_median$Rep == 4)]
cdf_data_median$relfluorescence[which(cdf_data_median$Rep == 5)] <- (cdf_data_median$`Median Fluorescence`/cdf_data_median$`Median Fluorescence`[which(cdf_data_median$Condition == "wt" & cdf_data_median$Rep == 5)])[which(cdf_data_median$Rep == 5)]

# Select genes of interest (mTOR related)
cdf_data_median_sub_2 <- subset(cdf_data_median, Condition != "wt" & Condition != "Lztr1" & Condition != "Smarcc1")

cdf_data_median_sub_2$Condition <- factor(cdf_data_median_sub_2$Condition, 
                                        levels = c("wt + mock", "Flcn",  "Lamtor4","Tsc1","Tsc2"), 
                                        labels = c("mock", "Flcn",  "Lamtor4","Tsc1","Tsc2"),
                                        ordered = TRUE)


# Barplot of both datasets (wt in the background)
ggplot(cdf_data_median_sub, aes(x=Condition, y=relfluorescence)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", fill = "#9CBCD3", color = "#144976", width = 0.8) +
  geom_bar(data = cdf_data_median_sub_2, aes(x=Condition, y=relfluorescence),
           position = "dodge", stat = "summary", fun = "mean", fill = "#9CBCD3",  width = 0.7, alpha = 0.5) +
  geom_jitter(data = cdf_data_median_sub_2, aes(x=Condition, y=relfluorescence, shape = factor(Rep)), 
              position = position_jitter(0.1), color = "#A09A8D", alpha = 0.5) +
  stat_summary(data = cdf_data_median_sub_2, aes(x=Condition, y=relfluorescence),
               geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "#144976", alpha=0.5) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", fill = "#9CBCD3", color = "#144976", width = 0.8) +
  geom_jitter(aes(shape = factor(Rep)), position = position_jitter(0.1), color = "#4C4C4C") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "#144976") +
  # stat_summary(aes(x = Condition), fun = mean, col = "darkred", shape = 3) +
  # geom_point(aes(shape = factor(Rep), color = "all")) + 
  # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c','#fb9a99')) +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  theme_bw(base_size = 12) + 
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,2.95)) + #ylim(0,2.8) + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  ylab("Median Fluorescence\nrelative to wt") +
  xlab("Gene Knockout") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black"),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.48),  # Adjust the size as needed
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_blank()) +
  force_panelsizes(rows = unit(1.5, "in"),
                   cols = unit(1.7, "in"))
ggsave("./20231102_Median_Fluorescent_relative_to_wt_in_Med12KO_and_wt_Fig3H.pdf", width = 5, height = 4)

  

