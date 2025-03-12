# Plot FACS data for Fig 1G
# load libraries
library(ggplot2)
library(reshape2)
library(ggh4x) # for force_panelsize option
library(multcomp) # for statistics
library(dplyr)
library(rstatix)

# Load data of all replicates
# rep 1
wt <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_-pd03_002_Single Cells.csv"),"id"="wt")
wt_pd03 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_+pd03_001_Single Cells.csv"),"id"="wt_pd03")
Grb2 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Grb2_KO_003_Single Cells.csv"),"id"="Grb2_KO")
Ptpn11 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Ptpn11_KO_004_Single Cells.csv"),"id"="Ptpn11_KO")
Sox2 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Sox2_KO_005_Single Cells.csv"),"id"="Sox2_KO")
Tada1 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Tada1_KO_010_Single Cells.csv"),"id"="Tada1_KO")
Fam98b <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Fam98b_KO_006_Single Cells.csv"),"id"="Fam98b_KO")
Med12 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Med12_KO_007_Single Cells.csv"),"id"="Med12_KO")
Med24 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Med24_KO_008_Single Cells.csv"),"id"="Med24_KO")
Med25 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Med25_KO_009_Single Cells.csv"),"id"="Med25_KO")
Ikbkap <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Ikbkap_KO_011_Single Cells.csv"),"id"="Ikbkap_KO")
Elp3 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Elp3_KO_012_Single Cells.csv"),"id"="Elp3_KO")
Elp5 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Elp5_KO_013_Single Cells.csv"),"id"="Elp5_KO")
Kti12 <- data.frame(read.csv("./Data/export_ES Cells_Spry4Venus_Kti12_KO_014_Single Cells.csv"),"id"="Kti12_KO")
plot_data_rep1 <- rbind(wt[,c(10,7)], wt_pd03[,c(10,7)], Grb2[,c(10,7)], Ptpn11[,c(10,7)], Sox2[,c(10,7)], Tada1[,c(10,7)], Fam98b[,c(10,7)], 
                      Med12[,c(10,7)], Med24[,c(10,7)], Med25[,c(10,7)], Ikbkap[,c(10,7)], Elp3[,c(10,7)], Elp5[,c(10,7)], Kti12[,c(10,7)])
colnames(plot_data_rep1) <- c("Condition", "mVenus_Fluorescence")
plot_data_rep1$rep <- "1"

# rep2 
wt <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_-pd03_002_Single Cells.csv"),"id"="wt")
wt_pd03 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_+pd03_001_Single Cells.csv"),"id"="wt_pd03")
Grb2 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Grb2_KO_003_Single Cells.csv"),"id"="Grb2_KO")
Ptpn11 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Ptpn11_KO_004_Single Cells.csv"),"id"="Ptpn11_KO")
Sox2 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Sox2_KO_005_Single Cells.csv"),"id"="Sox2_KO")
Tada1 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Tada1_KO_010_Single Cells.csv"),"id"="Tada1_KO")
Fam98b <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Fam98b_KO_006_Single Cells.csv"),"id"="Fam98b_KO")
Med12 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Med12_KO_007_Single Cells.csv"),"id"="Med12_KO")
Med24 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Med24_KO_008_Single Cells.csv"),"id"="Med24_KO")
Med25 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Med25_KO_009_Single Cells.csv"),"id"="Med25_KO")
Ikbkap <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Ikbkap_KO_011_Single Cells.csv"),"id"="Ikbkap_KO")
Elp3 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Elp3_KO_012_Single Cells.csv"),"id"="Elp3_KO")
Elp5 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Elp5_KO_013_Single Cells.csv"),"id"="Elp5_KO")
Kti12 <- data.frame(read.csv("./Data/Rep2/export_ES Cells_Spry4Venus_Kti12_KO_014_Single Cells.csv"),"id"="Kti12_KO")
plot_data_rep2 <- rbind(wt[,c(10,7)], wt_pd03[,c(10,7)], Grb2[,c(10,7)], Ptpn11[,c(10,7)], Sox2[,c(10,7)], Tada1[,c(10,7)], Fam98b[,c(10,7)], 
                        Med12[,c(10,7)], Med24[,c(10,7)], Med25[,c(10,7)], Ikbkap[,c(10,7)], Elp3[,c(10,7)], Elp5[,c(10,7)], Kti12[,c(10,7)])
colnames(plot_data_rep2) <- c("Condition", "mVenus_Fluorescence")
plot_data_rep2$rep <- "2"

# rep3
wt <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_-pd03_002_Single Cells.csv"),"id"="wt")
wt_pd03 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_+pd03_001_Single Cells.csv"),"id"="wt_pd03")
Grb2 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Grb2_KO_003_Single Cells.csv"),"id"="Grb2_KO")
Ptpn11 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Ptpn11_KO_004_Single Cells.csv"),"id"="Ptpn11_KO")
Sox2 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Sox2_KO_005_Single Cells.csv"),"id"="Sox2_KO")
Tada1 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Tada1_KO_010_Single Cells.csv"),"id"="Tada1_KO")
Fam98b <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Fam98b_KO_006_Single Cells.csv"),"id"="Fam98b_KO")
Med12 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Med12_KO_007_Single Cells.csv"),"id"="Med12_KO")
Med24 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Med24_KO_008_Single Cells.csv"),"id"="Med24_KO")
Med25 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Med25_KO_009_Single Cells.csv"),"id"="Med25_KO")
Ikbkap <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Ikbkap_KO_011_Single Cells.csv"),"id"="Ikbkap_KO")
Elp3 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Elp3_KO_012_Single Cells.csv"),"id"="Elp3_KO")
Elp5 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Elp5_KO_013_Single Cells.csv"),"id"="Elp5_KO")
Kti12 <- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_Kti12_KO_014_Single Cells.csv"),"id"="Kti12_KO")
Mock<- data.frame(read.csv("./Data/Rep3/d6/export_ES Cells_Spry4Venus_wt_mock_transfected_015_Single Cells.csv"),"id"="wt_mock")
plot_data_rep3 <- rbind(wt[,c(10,7)], wt_pd03[,c(10,7)], Mock[,c(10,7)], Grb2[,c(10,7)], Ptpn11[,c(10,7)], Sox2[,c(10,7)], Tada1[,c(10,7)], Fam98b[,c(10,7)], 
                   Med12[,c(10,7)], Med24[,c(10,7)], Med25[,c(10,7)], Ikbkap[,c(10,7)], Elp3[,c(10,7)], Elp5[,c(10,7)], Kti12[,c(10,7)])
colnames(plot_data_rep3) <- c("Condition", "mVenus_Fluorescence")
plot_data_rep3$rep <- "3"

# Compute Medians
plot_data_median_rep3 <- aggregate(x= plot_data_rep3$mVenus_Fluorescence,
                               by = list(plot_data_rep3$Condition),
                               FUN = median)
plot_data_median_rep2 <- aggregate(x= plot_data_rep2$mVenus_Fluorescence,
                               by = list(plot_data_rep2$Condition),
                               FUN = median)
plot_data_median_rep1 <- aggregate(x= plot_data_rep1$mVenus_Fluorescence,
                                   by = list(plot_data_rep1$Condition),
                                   FUN = median)

# Statistically compare wt to all other conditions
# (t tests+ multiple testing correction)
plot_data_median_rep3$rep <- "3"
plot_data_median_rep2$rep <- "2"
plot_data_median_rep1$rep <- "1"
data_for_statistics <- rbind(plot_data_median_rep1, plot_data_median_rep2, plot_data_median_rep3)

colnames(data_for_statistics) <- c("Condition", "Mean", "Replicate")
data_for_statistics <- subset(data_for_statistics, Condition != "wt_mock")

#convert Condition variable to factor 
data_for_statistics$Condition <- factor(data_for_statistics$Condition, levels = c("wt", "wt_pd03", "Elp3_KO", "Elp5_KO","Fam98b_KO","Grb2_KO","Ikbkap_KO", "Kti12_KO","Med12_KO","Med24_KO","Med25_KO","Ptpn11_KO","Sox2_KO","Tada1_KO"), ordered = TRUE)

# single t tests
data_for_statistics_wide <- reshape(data_for_statistics, direction = "wide", idvar = "Replicate", timevar = "Condition")
data_for_statistics_wide <- data_for_statistics_wide[,2:15]
controlCondition <- "Mean.wt" 

# Create a list to store the results
results_list <- list()

# Perform paired t-test for each condition against the control
for (condition in colnames(data_for_statistics_wide)) {
  if (condition != controlCondition) {
    # Perform paired t-test
    t_test_result <- t.test(data_for_statistics_wide[, condition], data_for_statistics_wide[, controlCondition], paired = TRUE, alternative = "less")
    
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



# plot barcharts of the medians, normalized to wt expression
plot_data_median_rep3$x <- plot_data_median_rep3$x/plot_data_median_rep3[13,2]
plot_data_median_rep2$x <- plot_data_median_rep2$x/plot_data_median_rep2[13,2]
plot_data_median_rep1$x <- plot_data_median_rep1$x/plot_data_median_rep1[13,2]
plot_data <- rbind(plot_data_median_rep1, plot_data_median_rep2, plot_data_median_rep3)
colnames(plot_data) <- c("Condition", "Mean", "Replicate")
plot_data$Condition = factor(plot_data$Condition, levels = c('wt', 'wt_mock', 'wt_pd03','Grb2_KO', 'Ptpn11_KO','Sox2_KO','Fam98b_KO','Med12_KO','Med24_KO','Med25_KO','Tada1_KO','Ikbkap_KO','Elp3_KO','Elp5_KO','Kti12_KO'))

plot_data_sub <- subset(plot_data, Condition != "wt" & Condition != "wt_mock")

ggplot(plot_data_sub, aes(x=Condition, y=Mean,  fill = Condition)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", color = "#144976", width = 0.7, size = 0.48) +
  geom_jitter(aes(shape = factor(Replicate)), position = position_jitter(0.1), color = "#4C4C4C") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "#144976") +
  scale_fill_manual(values = c('#B9BBBD','#87A9B9','#87A9B9','#B9BBBD','#B9BBBD','#8C88A4',
                               '#8C88A4','#8C88A4','#8C88A4','#BCB998','#BCB998','#BCB998','#BCB998')) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  scale_y_continuous(expand = c(0,0,0,0.1)) + 
  ylab("Median Expression\nrelative to wt")+
  xlab("Gene Knockout")+
  scale_x_discrete(labels=c("wt_pd03" = "wt + PD03", "Grb2_KO" = "Grb2","Ptpn11_KO" = "Ptpn11", "Sox2_KO" = "Sox2",
                            "Fam98b_KO" = "Fam98b","Med12_KO" = "Med12","Med24_KO" = "Med24","Med25_KO" = "Med25",
                            "Tada1_KO" = "Tada1","Ikbkap_KO" = "Ikbkap","Elp3_KO" = "Elp3","Elp5_KO" = "Elp5","Kti12_KO" = "Kti12"),
                   guide = guide_axis(angle = 90))+
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
                   cols = unit(6, "cm"))
ggsave("./Plots/231129_Bar_Median_Spry4_Expression_allReplicates.pdf", width = 8, height = 4)


