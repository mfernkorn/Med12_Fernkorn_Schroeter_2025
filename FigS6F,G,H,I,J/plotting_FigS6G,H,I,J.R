# Plot FigS6G-J
#load libraries
library(reshape2)
library(ggplot2)
library(ggh4x)
library(zoo) # for rollmean function
library(dplyr)
library(pROC)


# Load Tracks analyzed with Trackmate
path_list <- list.files(path = "./Tracking_intensities/", recursive = FALSE)
for (i in 1:length(path_list)){
  # Read in CellRanger Output
  temp <- read.csv(paste0("./Tracking_intensities/", path_list[i]), sep = ",")[-(1:3),c("TRACK_ID", "FRAME", "MEAN_INTENSITY_CH1", "MEAN_INTENSITY_CH2", "MEAN_INTENSITY_CH3")]
  # Assign Genotype
  if (grepl("Med12KO", path_list[i], fixed=TRUE)==TRUE){
    temp$Genotype <- "Med12 mutant"} else {
      temp$Genotype <- "wild type"
    }
  temp$Region <- sub(".*MMStack_(.*?).csv.*", "\\1", path_list[i])
  # Return as individual objects
  assign(sub(".*MMStack_(.*?).csv.*", "\\1", path_list[i]), temp)
}

# Combine Data of all tracks and add Time
Tracks <- rbind(wt_1, wt_2, wt_3,wt_4,wt_5,wt_6, Med12KO_1, Med12KO_2, Med12KO_3,Med12KO_4, Med12KO_5, Med12KO_6, Med12KO_7)
rm(wt_1, wt_2, wt_3,wt_4,wt_5,wt_6, Med12KO_1, Med12KO_2, Med12KO_3,Med12KO_4, Med12KO_5, Med12KO_6, Med12KO_7)
Tracks <- data.frame(Genotype = Tracks$Genotype, Region = Tracks$Region, TRACK_ID = Tracks$TRACK_ID, FRAME = as.numeric(Tracks$FRAME), MEAN_INTENSITY_CH1 = as.numeric(Tracks$MEAN_INTENSITY_CH1), MEAN_INTENSITY_CH2 = as.numeric(Tracks$MEAN_INTENSITY_CH2))
Tracks$TIME <- Tracks$FRAME/6

# Load Intensities of staining
path_list <- list.files(path = "./Staining_Intensities/", recursive = FALSE)
for (i in 1:length(path_list)){
  # Read in CellRanger Output
  temp <- read.csv(paste0("./Staining_Intensities/", path_list[i]),  sep = ",", row.names = "X")
  # Assign Genotype
  if (grepl("Med12KO", path_list[i], fixed=TRUE)==TRUE){
    temp$Genotype <- "Med12 mutant"} else {
      temp$Genotype <- "wild type"
    }
  temp$Region <- sub(".*MMStack_(.*?)_Star.*", "\\1", path_list[i])
  # Return as individual objects
  assign(sub(".*MMStack_(.*?)_Star.*", "\\1", path_list[i]), temp)
}

# Combine Intensities from all regions and make TRACk_ID and Channel usable
Intensities <- rbind(wt_1, wt_2, wt_3,wt_4,wt_5,wt_6, Med12KO_1, Med12KO_2, Med12KO_3,Med12KO_4, Med12KO_5, Med12KO_6, Med12KO_7)
rm(wt_1, wt_2, wt_3,wt_4,wt_5,wt_6, Med12KO_1, Med12KO_2, Med12KO_3,Med12KO_4, Med12KO_5, Med12KO_6, Med12KO_7)
Intensities <- subset(Intensities, Area > 40) # Filter small cells picked up by stardist
Intensities$TRACK_ID <- sub(".*:", "", Intensities$Label)
Intensities$Channel <- sub("-.*", "", sub(".*C", "", Intensities$Label))
Intensities$Tracked <- "segmented"
Intensities$Tracked[as.numeric(Intensities$TRACK_ID) >= 0] <- "tracked"
# Plot intensities
# Rearrange Data for plotting
Intensities_wide <- reshape(Intensities, timevar = "Channel",
        idvar = c("Area","TRACK_ID","Region","Genotype","Tracked"),
        direction = "wide")

# Select Threshold to classify cells
SOX17threshold <- 1000
NANOGthreshold <- 1100


# Classify Cells
for (i in 1:length(Intensities_wide$Genotype)){
  if (Intensities_wide$Mean.2[i] > SOX17threshold && Intensities_wide$Mean.3[i] > NANOGthreshold){
    Intensities_wide$class[i] <- "double positive"
  } else if (Intensities_wide$Mean.2[i] < SOX17threshold && Intensities_wide$Mean.3[i] > NANOGthreshold){
    Intensities_wide$class[i] <- "Epi"
  } else if (Intensities_wide$Mean.2[i] > SOX17threshold && Intensities_wide$Mean.3[i] < NANOGthreshold){
    Intensities_wide$class[i] <- "PrE"
  } else {
    Intensities_wide$class[i] <- "double negative"
  }
}
# transfer Celltype Information to Tracks
Tracks_selected <- merge(Tracks, Intensities_wide[,c("class", "TRACK_ID", "Genotype","Region")], by = c("TRACK_ID","Genotype","Region"))
Tracks_selected$TRACK_ID_Genotype_Region <- paste(Tracks_selected$TRACK_ID, Tracks_selected$Region, sep = "_")

Tracks_selected <- Tracks_selected[order(Tracks_selected$FRAME),]
Tracks_selected <- Tracks_selected %>% group_by(class,TRACK_ID,Genotype,Region) %>% 
  mutate(rollavg=rollmean(MEAN_INTENSITY_CH2, k = 7, fill = c(NA,NA,NA)))

# Plot Induction efficiency over time by cell, split by Genotype (Fig S6G)
Tracks_selected <- subset(Tracks_selected, class == "Epi" | class == "PrE")
Tracks_selected$Genotype <- factor(Tracks_selected$Genotype, levels = c("wild type", "Med12 mutant"), ordered = TRUE)
ggplot(Tracks_selected, aes(x=TIME, y=rollavg, group = TRACK_ID_Genotype_Region, color = class)) +
  geom_line() +
  scale_color_manual(values = c("#00A76E", "#D86FAA"))+
  ylab("mCherry Fluorescence") + 
  xlab("Time /h") +
  geom_vline(xintercept = 10, linetype = "dotted", color = "darkgrey", linewidth = 0.5) +
  theme_bw(base_size = 10) +
  facet_wrap(vars(Genotype), nrow = 2) + 
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
  force_panelsizes(rows = unit(2.5, "cm"),
                   cols = unit(2.5, "cm"))
ggsave("./Plots/20231207_Tracks_2_Rep9.pdf", width =8, height = 6)

# Plot Induction levels at specific timepoint (Figure S6H)
Tracks_selected_timepoint <- subset(Tracks_selected, TIME == 10)
ggplot(Tracks_selected_timepoint, aes(x= class, y=MEAN_INTENSITY_CH2)) +
  geom_violin(aes(fill=class)) +
  scale_fill_manual(values = c("#00A76E", "#D86FAA"))+
  geom_jitter(size = 0.8, width = 0.3) +
  xlab("Celltype")+
  ylab("mCherry Fluorescence")+
  facet_wrap(~Genotype, nrow = 2) +
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
  force_panelsizes(rows = unit(2.5, "cm"),
                   cols = unit(2.5, "cm"))
ggsave("./Plots/20231207_Induction_Level_10h_Rep9.pdf", width =5, height = 4)



# Start ROC analysis by calculating TPR and FPR for different thresholds
a <- subset(Tracks_selected_timepoint, Genotype == "wild type")[,c("MEAN_INTENSITY_CH2", "class")]

myRoc <- roc(response = a$class, predictor = a$MEAN_INTENSITY_CH2, positive = 'PrE')
plot(myRoc)
auc(myRoc)
coords(myRoc, "best", best.method="youden")

b <- subset(Tracks_selected_timepoint, Genotype == "Med12 mutant")[,c("MEAN_INTENSITY_CH2", "class")]

myRoc <- roc(response = b$class, predictor = b$MEAN_INTENSITY_CH2, positive = 'PrE')
plot(myRoc)
auc(myRoc)
coords(myRoc, "best", best.method="youden")




# Calculate AUC over time
frames = as.numeric(levels(factor(Tracks_selected$FRAME)))
AUC <- data.frame(FRAME = c(), wt = c(), mutant= c(), wt_threshold = c(), mutant_threshold = c())
for (i in 1:length(frames)){
  Tracks_selected_timepoint <- subset(Tracks_selected, FRAME == frames[i])
  
  wt <- subset(Tracks_selected_timepoint, Genotype == "wild type")[,c("MEAN_INTENSITY_CH2", "class")]
  myRoc <- roc(response = wt$class, predictor = wt$MEAN_INTENSITY_CH2, positive = 'PrE')
  auc_wt <- auc(myRoc)
  treshold_wt <- coords(myRoc, "best", best.method="youden")$threshold[1]
  
  KO <- subset(Tracks_selected_timepoint, Genotype == "Med12 mutant")[,c("MEAN_INTENSITY_CH2", "class")]
  myRoc <- roc(response = KO$class, predictor = KO$MEAN_INTENSITY_CH2, positive = 'PrE')
  auc_KO <- auc(myRoc)
  treshold_KO <- coords(myRoc, "best", best.method="youden")$threshold[1]
  
  AUC <- rbind(AUC, data.frame(FRAME = i, wt = auc_wt, mutant= auc_KO, wt_threshold = treshold_wt, mutant_threshold = treshold_KO))
}
AUC_long <- melt(AUC[,1:3], id.vars = "FRAME")
AUC_long$TIME <- AUC_long$FRAME/6

# Figure S6I
ggplot(AUC_long, aes(x=TIME, y=value, color= variable)) +
  geom_point(size = 0.5) +
  xlab("Time [h]")+
  ylab("AUC")+
  xlim(0,19) +
  scale_color_manual(values = c("#2E6A67", "#E98624"))+
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
  force_panelsizes(rows = unit(2, "cm"),
                   cols = unit(3, "cm"))
ggsave("./Plots/20231207_AUC_over_Time_Rep9.pdf", width =8, height = 6)

# Plor Figure S6J
threshold_long <- melt(AUC[,c(1,4,5)], id.vars = "FRAME")
threshold_long$TIME <- threshold_long$FRAME/6

ggplot(threshold_long, aes(x=TIME, y=value, color= variable)) +
  geom_point(size = 0.5) +
  xlab("Time [h]")+
  ylab("Optimal Threshold")+
  scale_color_manual(values = c("#2E6A67", "#E98624"))+
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
  force_panelsizes(rows = unit(2, "cm"),
                   cols = unit(3, "cm"))
ggsave("./Plots/20231207_Optimal_Threshold_over_Time_Rep9.pdf", width =8, height = 6)
