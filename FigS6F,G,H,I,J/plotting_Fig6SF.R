# Plot Fig S6G
# Load libraries
library(ggplot2)
library(ggh4x)
library(dplyr)

# Read in FACS data from one exploratory replicate
# Rep 9
path_list <- list.files(path = "./mCherry_Fluorescence_all_times/", recursive = FALSE)
mCherry_data_rep9 <- data.frame(X = c(), Area = c(), Mean= c(), StdDev = c(), IntDen = c(), RawIntDen = c(), Genotype = c(), Region = c())
for (i in 1:length(path_list)){
  temp <- read.csv(paste0("./mCherry_Fluorescence_all_times/", path_list[i]), sep = ",")
  # Assign Genotype
  if (grepl("Med12KO", path_list[i], fixed=TRUE)==TRUE){
    temp$Genotype <- "Med12 KO"} else {
      temp$Genotype <- "wt"
    }
  temp$Frame <- sub(".*tif_(.*?).csv", "\\1", path_list[i])
  temp$Region <- sub(".*MMStack_(.*?).ome.*", "\\1", path_list[i])
  
  mCherry_data_rep9 <- rbind(mCherry_data_rep9, temp)
}
mCherry_data_rep9$Rep <- 9 
mCherry_data_rep9$Frame <- as.numeric(mCherry_data_rep9$Frame)
mCherry_data_rep9 <- subset(mCherry_data_rep9, Area >30)

Means_rep9 <- mCherry_data_rep9 %>% group_by(Genotype, Frame) %>%
  summarise(Mean_over_Cells = mean(Mean), sd = sd(Mean,na.rm = TRUE))
Means_rep9$Rep <- 9


Means_rep9$Time <- Means_rep9$Frame/6
ggplot(Means_rep9, aes(x=Time, y=Mean_over_Cells, color = Genotype)) +
  geom_rect(aes(xmin = 0, xmax = 8, ymin = 4200, ymax = 4500),
            fill = "#E98624",
            alpha = 0.2, linewidth = 0
  )+
  geom_rect(aes(xmin = 4, xmax = 8, ymin = 4600, ymax = 4900),
            fill = "#2E6A67",
            alpha = 0.2, linewidth = 0
  )+
  geom_errorbar(aes(ymin=Mean_over_Cells-sd, ymax=Mean_over_Cells+sd),
                position=position_dodge(0.05), alpha = 0.4, linewidth = 0.3) +
  geom_point(size = 1) +
  # facet_wrap(vars(Rep), scales = "free_y") +
  xlab("Time [h]")+
  ylab("mCherry Expression")+
  scale_color_manual(values = c("#E98624", "#2E6A67"))+
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.24),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.24),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(3, "cm"))
ggsave("./Plots/20231129_Mean_Induction_Efficiency_rep9.pdf", width = 4, height = 4)
