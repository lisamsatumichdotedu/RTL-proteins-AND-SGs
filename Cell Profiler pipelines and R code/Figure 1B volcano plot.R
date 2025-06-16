library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(ggprism)
library(patchwork)
library(openxlsx)

# Reading and cleaning data based on Row 14 (log2(abundance ratio))
  MSdata <- read.xlsx(file.choose())
  MSdata <- MSdata[!is.na(MSdata[,14]),]
  
# Subsetting data to only keep gene description, log2(abundance ratio) and adjusted p-value
  MSdata <- MSdata[,c(3,14,17)]
  colnames(MSdata) <- c("Gene","Ratio","P_value")

# Mutating column 1 to containing only gene names
  genename <- sub(".*GN=", "", MSdata[,1])
  genename <- sub(" PE=.*", "", genename)
  MSdata <- MSdata %>% mutate(Gene = genename)
  
# Log10 transforming p-values
  MSdata <- MSdata %>% mutate(P_value = -log10(P_value))

# Choosing color scheme for data points (transparent black- no hits, red - upregulated, blue - downregulated)
  MSdata$label <- ""
  MSdata$color <- "grey"
  MSdata$alpha <- 0.2
  MSdata$boxlabel <- ""

  MSdata$label[MSdata$Ratio > log2(2.00) & MSdata$P_value > -log10(0.05)] <- MSdata$Gene[MSdata$Ratio > log2(2.00) & MSdata$P_value > -log10(0.05)]
  MSdata$color[MSdata$Ratio > log2(2.00) & MSdata$P_value > -log10(0.05)] <- "red"
  MSdata$alpha[MSdata$Ratio > log2(2.00) & MSdata$P_value > -log10(0.05)] <- 1.0

  MSdata$label[MSdata$Ratio < -log2(2.00) & MSdata$P_value > -log10(0.05)] <- MSdata$Gene[MSdata$Ratio < -log2(2.00) & MSdata$P_value > -log10(0.05)]
  MSdata$color[MSdata$Ratio < -log2(2.00) & MSdata$P_value > -log10(0.05)] <- "blue"
  MSdata$alpha[MSdata$Ratio < -log2(2.00) & MSdata$P_value > -log10(0.05)] <- 1.0

# Conditional formatting for PEG10
  MSdata$label[MSdata$Gene == "PEG10"] <- ""
  MSdata$boxlabel[MSdata$Gene == "PEG10"] <- "PEG10"
  
# Creating volcano plot without grid highlighting all significantly changed data points
  ggplot(MSdata, aes(x = Ratio, y = P_value)) + removeGrid() +
# Dot plot
  geom_point(color = MSdata$color, alpha = MSdata$alpha, size = 2.5) + 
# Drawing lines at 0 and other significant values
  geom_hline(yintercept = -log10(0.05) , linetype = "dotted", alpha = 0.6, size = 0.75) +
  geom_vline(xintercept = log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = -log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.65) +
# To prevent text overlap with points, set label alpha, face and size
  geom_text_repel(label = MSdata$label, max.iter = 20000, color = MSdata$color, size = 5.5) + 
# X-axis labels and ticks
  scale_x_continuous(expression(log[2]~(RTL8~KO/Control)),limits = c(-2,2), breaks = seq(-2,2,1), guide = guide_prism_minor()) +
# Y-axis labels and ticks
  scale_y_continuous(expression(-log[10]~(adjusted~p-value)), limits = c(0,15), breaks = seq(0,15,5), guide = guide_prism_minor()) +
# Specifying axes thickness, label positions, axis label size, major and minor tick length and sizes, and plot size
  theme(axis.line.x = element_line(color="black", size = 0.75), axis.title.x = element_text(vjust = 2, size = 20), axis.text = element_text(size = 17.5, face = "plain", color = "black"),
        axis.line.y = element_line(color="black", size = 0.75), axis.title.y = element_text(vjust = -0.5, size = 20),
        axis.ticks = element_line(size = 0.75), axis.ticks.length = unit(0.2,"cm"), prism.ticks.length = unit(0.15,"cm"),
        panel.background = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0.2, #top
                                                                                                0.15, #right
                                                                                               -0.3, #bottom
                                                                                               -0.25), #left
                                                                                               "cm")) +
# Drawing label around PEG10                                                                                                                                
  geom_label_repel(label = MSdata$boxlabel, max.iter = 20000, color = MSdata$color, label.size = 1, size = 5.5, fontface = "bold")    
# Saving plot as .jpg
  ggsave("TMT MS RTL8 KO vs Control.jpg", units="cm", width=15, height=15, dpi = 600)