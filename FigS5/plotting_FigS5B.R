# PLot FigS5B
# load libraries
library(ggplot2)
library(ggh4x) # for force_panelsize option

# Import FACS data, exported with flowjo
filedir <- "./Exported_Flowjo_Files/Med12KOline/"
file_names <- dir(filedir)
names_vec <- c("Med12 KO",  "Med12 KO + mock", "Med12 KO + Med12L KO","wt + mock", "wt")
file_directories <- paste0(filedir, file_names)


cdf_data <- data.frame()
for (i in 1:length(file_directories)){
  cdf_data <- rbind(cdf_data, data.frame(as.numeric(read.csv2(file_directories[i], sep = ",")[,7]), names_vec[i]))
}
colnames(cdf_data) <- c("Fluorescence","Condition")
cdf_data$Condition <- factor(cdf_data$Condition, 
                             levels = c("Med12 KO + Med12L KO", "Med12 KO",  "Med12 KO + mock", "wt", "wt + mock"), 
                             ordered = TRUE)

ggplot(cdf_data, aes(Fluorescence, colour = Condition)) + stat_ecdf() +
  scale_x_log10() +
  coord_cartesian(xlim = c(100, 50000)) +
  theme_bw(base_size = 12) + theme(aspect.ratio = 1) +
  xlab("Fluorescence [au]") + ylab("Relative abundance") +
  scale_color_manual(values=c('#6a3d9a','#a6cee3','#1f78b4','#b2df8a','#33a02c')) +
  force_panelsizes(rows = unit(2, "in"),
                   cols = unit(2, "in"))
ggsave("./20230228_CDF_Med12KOLine.pdf", width = 5, height = 2.5)
# Modified colors in Illustrator to fit final color scheme