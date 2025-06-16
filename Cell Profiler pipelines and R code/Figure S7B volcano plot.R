library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(ggprism)
library(patchwork)
library(openxlsx)
library(ggpubr) # For arranging graphs onto the same plot
library(ggvenn)

# Reading MS data and SG list
MSdata <- read.xlsx(file.choose())
SG_proteins <- read.csv(file.choose())
MSdata <- MSdata[MSdata$`#.PSMs`>=5,]

## Calculating enrichment values for SG proteins across data set ##

total_tier_1 <- count(SG_proteins[SG_proteins$Tier ==1,])
total_tier_2 <- count(SG_proteins[SG_proteins$Tier ==2,])
total_tier_3 <- count(SG_proteins[SG_proteins$Tier ==3,])

MSdata$SG_protein_tier <- case_when(
  MSdata$Gene.Symbol %in% SG_proteins$Protein[SG_proteins$Tier == 1] ~ 1,
  MSdata$Gene.Symbol %in% SG_proteins$Protein[SG_proteins$Tier == 2] ~ 2,
  MSdata$Gene.Symbol %in% SG_proteins$Protein[SG_proteins$Tier == 3] ~ 3,
  TRUE ~ 0
)

MS_tier_1 <- count(MSdata[MSdata$SG_protein_tier ==1,])
MS_tier_2 <- count(MSdata[MSdata$SG_protein_tier ==2,])
MS_tier_3 <- count(MSdata[MSdata$SG_protein_tier ==3,])

percent_tier_1 <- (MS_tier_1 / total_tier_1)*100
percent_tier_2 <- (MS_tier_2 / total_tier_2)*100
percent_tier_3 <- (MS_tier_3 / total_tier_3)*100

tier_1_in_data <- (MS_tier_1/nrow(MSdata))*100
tier_2_in_data <- (MS_tier_2/nrow(MSdata))*100
tier_3_in_data <- (MS_tier_3/nrow(MSdata))*100
total_SG_enrichment <- ((MS_tier_1+MS_tier_2+MS_tier_3)/nrow(MSdata))*100

print(paste("The percentage of Tier 1 SG proteins pulled down are ", round(percent_tier_1,2),"%. ", 
            "The percentage of Tier 2 SG proteins pulled down are ", round(percent_tier_2,2),"%. ",
            "The percentage of Tier 3 SG proteins pulled down are ", round(percent_tier_3,2),"%."))
print(paste("Tier 1 SG proteins represent ", round(tier_1_in_data,2),"% of all proteins in the dataset. ",
            "Tier 2 SG proteins represent ", round(tier_2_in_data,2),"% of all proteins in the dataset. ",
            "Tier 3 SG proteins represent ", round(tier_3_in_data,2),"% of all proteins in the dataset. ",
            "SG proteins represent ", round(total_SG_enrichment,2),"% of all proteins in the dataset."))

## The following block is for PEG10 KO/NTC analysis ##

# Cleaning data based on Column 7 (log2(abundance ratio PEG10 KO/NT)) and column 10 (p-value)
PEG10_NTC_data <- MSdata[complete.cases(MSdata[, c(7, 10)]), ]

# Subsetting data to only keep gene name, log2(abundance ratio) and adjusted p-value
PEG10_NTC_data <- PEG10_NTC_data[,c(1,7,10)]
colnames(PEG10_NTC_data) <- c("Gene","Ratio","P_value")

# Log10 transforming p-values
PEG10_NTC_data <- PEG10_NTC_data %>% mutate(P_value = -log10(P_value))

# Finding genes that are significantly up or down and then conditionally formatting (transparent black- no hits, red - upregulated, blue - downregulated)

PEG10_NTC_data$sig_change[PEG10_NTC_data$P_value >= -log10(0.05)]  <- if_else(
  PEG10_NTC_data$Ratio[PEG10_NTC_data$P_value >= -log10(0.05)] >= log2(2.00), "up",if_else(
  PEG10_NTC_data$Ratio[PEG10_NTC_data$P_value >= -log10(0.05)] <= -log2(2.00), "down", NA
  )
)

PEG10_NTC_data$color <- case_when(
  PEG10_NTC_data$sig_change == "up"  ~ "red",
  PEG10_NTC_data$sig_change == "down" ~ "blue",
  TRUE ~ "grey"
)

# Comparing against SG data set to assign tier and conditionally formatting labels (size & fontface) and points (size ,shape & alpha) based on tier

PEG10_NTC_data$SG_protein_tier[!is.na(PEG10_NTC_data$sig_change)] <- case_when(
  PEG10_NTC_data$Gene[!is.na(PEG10_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 1] ~ 1,
  PEG10_NTC_data$Gene[!is.na(PEG10_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 2] ~ 2,
  PEG10_NTC_data$Gene[!is.na(PEG10_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 3] ~ 3,
  PEG10_NTC_data$Gene[!is.na(PEG10_NTC_data$sig_change)] == "ATXN10" ~ 0,
  TRUE ~ NA
)

PEG10_NTC_data$label <- case_when(
  PEG10_NTC_data$SG_protein_tier > 0 ~ PEG10_NTC_data$Gene,
  PEG10_NTC_data$Gene == "ATXN10" ~ "ATXN10*",
  TRUE ~ NA
)

PEG10_NTC_data$font_face <- case_when(
  PEG10_NTC_data$SG_protein_tier == 1 ~ "bold",
  PEG10_NTC_data$SG_protein_tier == 2 ~ "italic",
  TRUE ~ "plain"
)

PEG10_NTC_data$font_size <- case_when(
  PEG10_NTC_data$SG_protein_tier == 1 ~ 7,
  PEG10_NTC_data$SG_protein_tier == 2 ~ 5.5,
  PEG10_NTC_data$SG_protein_tier == 3 ~ 4,
  PEG10_NTC_data$Gene == "ATXN10" ~ 4,
  TRUE ~ NA
)

PEG10_NTC_data$point_size <- case_when(
  PEG10_NTC_data$SG_protein_tier == 1 ~ 6.4,
  PEG10_NTC_data$SG_protein_tier == 2 ~ 3.2,
  PEG10_NTC_data$SG_protein_tier == 3 ~ 1.6,
  PEG10_NTC_data$Gene == "ATXN10" ~ 1.6,
  TRUE ~ 0.8
)

PEG10_NTC_data$point_shape <- case_when(
  PEG10_NTC_data$SG_protein_tier == 1 ~ 18,
  PEG10_NTC_data$SG_protein_tier == 2 ~ 17,
  PEG10_NTC_data$SG_protein_tier == 3 ~ 15,
  TRUE ~ 16
)

PEG10_NTC_data$alpha <- case_when(
  PEG10_NTC_data$SG_protein_tier == 1 ~ 1,
  PEG10_NTC_data$SG_protein_tier == 2 ~ 0.75,
  PEG10_NTC_data$SG_protein_tier == 3 ~ 0.5,
  PEG10_NTC_data$sig_change == TRUE ~ 0.25,
  PEG10_NTC_data$Gene == "ATXN10" ~ 0.5,
  TRUE ~ 0.1
)

# Creating volcano plot without grid highlighting all significantly changed data points
set.seed(42)
PEG10_NTC_plot <-ggplot(PEG10_NTC_data, aes(x = Ratio, y = P_value)) + removeGrid() +
  # Dot plot
  geom_point(color = PEG10_NTC_data$color, alpha = PEG10_NTC_data$alpha, size = PEG10_NTC_data$point_size, shape = PEG10_NTC_data$point_shape) + 
  # Drawing lines at 0 and other significant values
  geom_hline(yintercept = -log10(0.05) , linetype = "dotted", alpha = 0.6, size = 0.75) +
  geom_vline(xintercept = log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = -log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.65, alpha = 0.8) +
  # To prevent text overlap with points, set label properties and padding (to reduce collisions)
  geom_text_repel(label = PEG10_NTC_data$label, fontface = PEG10_NTC_data$font_face, size = PEG10_NTC_data$font_size, max.iter = 100000, 
                  color = PEG10_NTC_data$color, max.overlaps = Inf, box.padding = 1.2, point.padding = 0.2,) + 
  # X-axis labels and ticks
  scale_x_continuous(expression(log[2]~(PEG10~KO/NTC)),limits = c(-7,7), breaks = seq(-7,7,2), guide = guide_prism_minor()) +
  # Y-axis labels and ticks
  scale_y_continuous(expression(-log[10]~(adjusted~p-value)), limits = c(0,15.5), breaks = seq(0,15.5,5), guide = guide_prism_minor()) +
  # Specifying axes thickness, label positions, axis label size, major and minor tick length and sizes, and plot size
  theme(axis.line.x = element_line(color="black", size = 0.75), axis.title.x = element_text(vjust = 2, size = 20), axis.text = element_text(size = 17.5, face = "plain", color = "black"),
        axis.line.y = element_line(color="black", size = 0.75), axis.title.y = element_text(vjust = -0.5, size = 20),
        axis.ticks = element_line(size = 0.75), axis.ticks.length = unit(0.2,"cm"), prism.ticks.length = unit(0.15,"cm"),
        panel.background = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0.2, #top
                                                                                                 0.15, #right
                                                                                                 -0.3, #bottom
                                                                                                 -0.25), #left
                                                                                               "cm"))


## The following block is for RTL8 KO/NTC analysis ##

# Cleaning data based on Column 8 (log2(abundance ratio RTL8 KO/NT)) and 11 (p-value)
RTL8_NTC_data <- MSdata[complete.cases(MSdata[, c(8, 11)]), ]

# Subsetting data to only keep gene name, log2(abundance ratio) and adjusted p-value
RTL8_NTC_data <- RTL8_NTC_data[,c(1,8,11)]
colnames(RTL8_NTC_data) <- c("Gene","Ratio","P_value")

# Log10 transforming p-values
RTL8_NTC_data <- RTL8_NTC_data %>% mutate(P_value = -log10(P_value))

# Finding genes that are significantly up or down and then conditionally formatting (transparent black- no hits, red - upregulated, blue - downregulated)

RTL8_NTC_data$sig_change[RTL8_NTC_data$P_value >= -log10(0.05)]  <- if_else(
  RTL8_NTC_data$Ratio[RTL8_NTC_data$P_value >= -log10(0.05)] >= log2(2.00), "up",if_else(
    RTL8_NTC_data$Ratio[RTL8_NTC_data$P_value >= -log10(0.05)] <= -log2(2.00), "down", NA
  )
)

RTL8_NTC_data$color <- case_when(
  RTL8_NTC_data$sig_change == "up"  ~ "red",
  RTL8_NTC_data$sig_change == "down" ~ "blue",
  TRUE ~ "grey"
)

# Comparing against SG data set to assign tier and conditionally formatting labels (size & fontface) and points (size ,shape & alpha) based on tier

RTL8_NTC_data$SG_protein_tier[!is.na(RTL8_NTC_data$sig_change)] <- case_when(
  RTL8_NTC_data$Gene[!is.na(RTL8_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 1] ~ 1,
  RTL8_NTC_data$Gene[!is.na(RTL8_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 2] ~ 2,
  RTL8_NTC_data$Gene[!is.na(RTL8_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 3] ~ 3,
  RTL8_NTC_data$Gene[!is.na(RTL8_NTC_data$sig_change)] == "ATXN10" ~ 0,
  TRUE ~ NA
)

RTL8_NTC_data$label <- case_when(
  RTL8_NTC_data$SG_protein_tier > 0 ~ RTL8_NTC_data$Gene,
  RTL8_NTC_data$Gene == "ATXN10" ~ "ATXN10*",
  TRUE ~ NA
)

RTL8_NTC_data$font_face <- case_when(
  RTL8_NTC_data$SG_protein_tier == 1 ~ "bold",
  RTL8_NTC_data$SG_protein_tier == 2 ~ "italic",
  TRUE ~ "plain"
)

RTL8_NTC_data$font_size <- case_when(
  RTL8_NTC_data$SG_protein_tier == 1 ~ 7,
  RTL8_NTC_data$SG_protein_tier == 2 ~ 5.5,
  RTL8_NTC_data$SG_protein_tier == 3 ~ 4,
  RTL8_NTC_data$Gene == "ATXN10" ~ 4,
  TRUE ~ NA
)

RTL8_NTC_data$point_size <- case_when(
  RTL8_NTC_data$SG_protein_tier == 1 ~ 6.4,
  RTL8_NTC_data$SG_protein_tier == 2 ~ 3.2,
  RTL8_NTC_data$SG_protein_tier == 3 ~ 1.6,
  RTL8_NTC_data$Gene == "ATXN10" ~ 1.6,
  TRUE ~ 0.8
)

RTL8_NTC_data$point_shape <- case_when(
  RTL8_NTC_data$SG_protein_tier == 1 ~ 18,
  RTL8_NTC_data$SG_protein_tier == 2 ~ 17,
  RTL8_NTC_data$SG_protein_tier == 3 ~ 15,
  TRUE ~ 16
)

RTL8_NTC_data$alpha <- case_when(
  RTL8_NTC_data$SG_protein_tier == 1 ~ 1,
  RTL8_NTC_data$SG_protein_tier == 2 ~ 0.75,
  RTL8_NTC_data$SG_protein_tier == 3 ~ 0.5,
  RTL8_NTC_data$sig_change == TRUE ~ 0.25,
  RTL8_NTC_data$Gene == "ATXN10" ~ 0.5,
  TRUE ~ 0.1
)

# Creating volcano plot without grid highlighting all significantly changed data points
set.seed(35)
RTL8_NTC_plot <-ggplot(RTL8_NTC_data, aes(x = Ratio, y = P_value)) + removeGrid() +
  # Dot plot
  geom_point(color = RTL8_NTC_data$color, alpha = RTL8_NTC_data$alpha, size = RTL8_NTC_data$point_size, shape = RTL8_NTC_data$point_shape) + 
  # Drawing lines at 0 and other significant values
  geom_hline(yintercept = -log10(0.05) , linetype = "dotted", alpha = 0.6, size = 0.75) +
  geom_vline(xintercept = log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = -log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.65, alpha = 0.8) +
  # To prevent text overlap with points, set label properties and padding (to reduce collisions)
  geom_text_repel(label = RTL8_NTC_data$label, fontface = RTL8_NTC_data$font_face, size = RTL8_NTC_data$font_size, max.iter = 100000, 
                  color = RTL8_NTC_data$color, max.overlaps = Inf, box.padding = 1.2, point.padding = 0.2,) + 
  # X-axis labels and ticks
  scale_x_continuous(expression(log[2]~(RTL8~KO/NTC)),limits = c(-7,7), breaks = seq(-7,7,2), guide = guide_prism_minor()) +
  # Y-axis labels and ticks
  scale_y_continuous(expression(-log[10]~(adjusted~p-value)), limits = c(0,15.5), breaks = seq(0,15.5,5), guide = guide_prism_minor()) +
  # Specifying axes thickness, label positions, axis label size, major and minor tick length and sizes, and plot size
  theme(axis.line.x = element_line(color="black", size = 0.75), axis.title.x = element_text(vjust = 2, size = 20), axis.text = element_text(size = 17.5, face = "plain", color = "black"),
        axis.line.y = element_line(color="black", size = 0.75), axis.title.y = element_text(vjust = -0.5, size = 20),
        axis.ticks = element_line(size = 0.75), axis.ticks.length = unit(0.2,"cm"), prism.ticks.length = unit(0.15,"cm"),
        panel.background = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0.2, #top
                                                                                                 0.15, #right
                                                                                                 -0.3, #bottom
                                                                                                 -0.25), #left
                                                                                               "cm"))


## The following block is for UBQLN2 KO/NTC analysis ##

# Cleaning data based on Column 9 (log2(abundance ratio UBQLN2 KO/NT)) and 12 (p-value)
UBQLN2_NTC_data <- MSdata[complete.cases(MSdata[, c(9, 12)]), ]

# Subsetting data to only keep gene name, log2(abundance ratio) and adjusted p-value
UBQLN2_NTC_data <- UBQLN2_NTC_data[,c(1,9,12)]
colnames(UBQLN2_NTC_data) <- c("Gene","Ratio","P_value")

# Log10 transforming p-values
UBQLN2_NTC_data <- UBQLN2_NTC_data %>% mutate(P_value = -log10(P_value))

# Finding genes that are significantly up or down and then conditionally formatting (transparent black- no hits, red - upregulated, blue - downregulated)

UBQLN2_NTC_data$sig_change[UBQLN2_NTC_data$P_value >= -log10(0.05)]  <- if_else(
  UBQLN2_NTC_data$Ratio[UBQLN2_NTC_data$P_value >= -log10(0.05)] >= log2(2.00), "up",if_else(
    UBQLN2_NTC_data$Ratio[UBQLN2_NTC_data$P_value >= -log10(0.05)] <= -log2(2.00), "down", NA
  )
)

UBQLN2_NTC_data$color <- case_when(
  UBQLN2_NTC_data$sig_change == "up"  ~ "red",
  UBQLN2_NTC_data$sig_change == "down" ~ "blue",
  TRUE ~ "grey"
)

# Comparing against SG data set to assign tier and conditionally formatting labels (size & fontface) and points (size ,shape & alpha) based on tier

UBQLN2_NTC_data$SG_protein_tier[!is.na(UBQLN2_NTC_data$sig_change)] <- case_when(
  UBQLN2_NTC_data$Gene[!is.na(UBQLN2_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 1] ~ 1,
  UBQLN2_NTC_data$Gene[!is.na(UBQLN2_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 2] ~ 2,
  UBQLN2_NTC_data$Gene[!is.na(UBQLN2_NTC_data$sig_change)] %in% SG_proteins$Protein[SG_proteins$Tier == 3] ~ 3,
  UBQLN2_NTC_data$Gene[!is.na(UBQLN2_NTC_data$sig_change)] == "ATXN10" ~ 0,
  TRUE ~ NA
)

UBQLN2_NTC_data$label <- case_when(
  UBQLN2_NTC_data$SG_protein_tier > 0 ~ UBQLN2_NTC_data$Gene,
  UBQLN2_NTC_data$Gene == "ATXN10" ~ "ATXN10*",
  TRUE ~ NA
)

UBQLN2_NTC_data$font_face <- case_when(
  UBQLN2_NTC_data$SG_protein_tier == 1 ~ "bold",
  UBQLN2_NTC_data$SG_protein_tier == 2 ~ "italic",
  TRUE ~ "plain"
)

UBQLN2_NTC_data$font_size <- case_when(
  UBQLN2_NTC_data$SG_protein_tier == 1 ~ 7,
  UBQLN2_NTC_data$SG_protein_tier == 2 ~ 5.5,
  UBQLN2_NTC_data$SG_protein_tier == 3 ~ 4,
  UBQLN2_NTC_data$Gene == "ATXN10" ~ 4,
  TRUE ~ NA
)

UBQLN2_NTC_data$point_size <- case_when(
  UBQLN2_NTC_data$SG_protein_tier == 1 ~ 6.4,
  UBQLN2_NTC_data$SG_protein_tier == 2 ~ 3.2,
  UBQLN2_NTC_data$SG_protein_tier == 3 ~ 1.6,
  UBQLN2_NTC_data$Gene == "ATXN10" ~ 1.6,
  TRUE ~ 0.8
)

UBQLN2_NTC_data$point_shape <- case_when(
  UBQLN2_NTC_data$SG_protein_tier == 1 ~ 18,
  UBQLN2_NTC_data$SG_protein_tier == 2 ~ 17,
  UBQLN2_NTC_data$SG_protein_tier == 3 ~ 15,
  TRUE ~ 16
)

UBQLN2_NTC_data$alpha <- case_when(
  UBQLN2_NTC_data$SG_protein_tier == 1 ~ 1,
  UBQLN2_NTC_data$SG_protein_tier == 2 ~ 0.75,
  UBQLN2_NTC_data$SG_protein_tier == 3 ~ 0.5,
  UBQLN2_NTC_data$sig_change == TRUE ~ 0.25,
  UBQLN2_NTC_data$Gene == "ATXN10" ~ 0.5,
  TRUE ~ 0.1
)

# Creating volcano plot without grid highlighting all significantly changed data points
set.seed(25)
UBQLN2_NTC_plot <-ggplot(UBQLN2_NTC_data, aes(x = Ratio, y = P_value)) + removeGrid() +
  # Dot plot
  geom_point(color = UBQLN2_NTC_data$color, alpha = UBQLN2_NTC_data$alpha, size = UBQLN2_NTC_data$point_size, shape = UBQLN2_NTC_data$point_shape) + 
  # Drawing lines at 0 and other significant values
  geom_hline(yintercept = -log10(0.05) , linetype = "dotted", alpha = 0.6, size = 0.75) +
  geom_vline(xintercept = log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = -log2(2.00), linetype = "dotted", alpha = 0.6, size = 0.75) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.65, alpha = 0.8) +
  # To prevent text overlap with points, set label properties and padding (to reduce collisions)
  geom_text_repel(label = UBQLN2_NTC_data$label, fontface = UBQLN2_NTC_data$font_face, size = UBQLN2_NTC_data$font_size, max.iter = 100000, 
                  color = UBQLN2_NTC_data$color, max.overlaps = Inf, box.padding = 1.2, point.padding = 0.2,) + 
  # X-axis labels and ticks
  scale_x_continuous(expression(log[2]~(UBQLN2~KO/NTC)),limits = c(-7,7), breaks = seq(-7,7,2), guide = guide_prism_minor()) +
  # Y-axis labels and ticks
  scale_y_continuous(expression(-log[10]~(adjusted~p-value)), limits = c(0,15.5), breaks = seq(0,15.5,5), guide = guide_prism_minor()) +
  # Specifying axes thickness, label positions, axis label size, major and minor tick length and sizes, and plot size
  theme(axis.line.x = element_line(color="black", size = 0.75), axis.title.x = element_text(vjust = 2, size = 20), axis.text = element_text(size = 17.5, face = "plain", color = "black"),
        axis.line.y = element_line(color="black", size = 0.75), axis.title.y = element_text(vjust = -0.5, size = 20),
        axis.ticks = element_line(size = 0.75), axis.ticks.length = unit(0.2,"cm"), prism.ticks.length = unit(0.15,"cm"),
        panel.background = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0.2, #top
                                                                                                 0.15, #right
                                                                                                 -0.3, #bottom
                                                                                                 -0.25), #left
                                                                                               "cm"))

# Making a blank legend plot

UBQLN2_NTC_data <- UBQLN2_NTC_data %>% add_row(Gene = "blank", point_shape = 0)
legend_plot <- ggplot(UBQLN2_NTC_data, aes(Ratio, P_value, shape = as.factor(point_shape))) + geom_point() +
  scale_shape_manual(values = c(16,15,17,18,0), 
                     labels = c("Non-SG protein (unlabeled, ATXN10 marked with *)", "Tier 3 SG protein ",
                                expression(paste(italic("Tier 2 "), "SG protein")),expression(paste(bold("Tier 1 "),"SG protein")), "hits of interest")) +
  guides(shape = guide_legend(override.aes = list(size = c(0.8,1.6,3.2,6.4,10)))) +
  theme_void() + # Remove axes and background
  theme(legend.position = "bottom", legend.text = element_text(size = 20, margin = margin(r=0.3, unit ="in")), 
        legend.spacing.x = unit(0.1, "cm")) + # Place legend at the bottom
  labs(shape = "")
  
legend <- get_legend(legend_plot)

# Arranging and saving plots as .jpg

plots <- ggarrange(PEG10_NTC_plot, RTL8_NTC_plot, UBQLN2_NTC_plot, nrow = 1, widths = c(1, 1,1))
final_plot <- ggarrange(plots, legend, ncol = 1, heights = c(10,1))
ggsave("Volcano plots.jpg", units ="in", width = 20, height = 11, dpi = 600) # these are arbitrary to ensure minimal overlap of labels


## Writing significantly changed SG proteins to a new tibble ##

# Downregulated genes

PEG10_down <- na.omit(PEG10_NTC_data$label[PEG10_NTC_data$sig_change == "down" & !is.na(PEG10_NTC_data$SG_protein_tier)])
RTL8_down <- na.omit(RTL8_NTC_data$label[RTL8_NTC_data$sig_change == "down" & !is.na(RTL8_NTC_data$SG_protein_tier)])
UBQLN2_down <- na.omit(UBQLN2_NTC_data$label[UBQLN2_NTC_data$sig_change == "down" & !is.na(UBQLN2_NTC_data$SG_protein_tier)])

down_genes <- list("RTL8 KO" = RTL8_down, "PEG10 KO" = PEG10_down, "UBQLN2 KO" = UBQLN2_down)

venn_down <- ggvenn(down_genes, fill_color = c("#D4F6FF","#FFF1D4","#FFD5F9"), show_elements = TRUE, label_sep = "\n")

# Upregulated genes

PEG10_up <- na.omit(PEG10_NTC_data$label[PEG10_NTC_data$sig_change == "up" & !is.na(PEG10_NTC_data$SG_protein_tier)])
RTL8_up <- na.omit(RTL8_NTC_data$label[RTL8_NTC_data$sig_change == "up" & !is.na(RTL8_NTC_data$SG_protein_tier)])
UBQLN2_up <- na.omit(UBQLN2_NTC_data$label[UBQLN2_NTC_data$sig_change == "up" & !is.na(UBQLN2_NTC_data$SG_protein_tier)])

up_genes <- list("RTL8 KO" = RTL8_up, "PEG10 KO" = PEG10_up, "UBQLN2 KO" = UBQLN2_up)

venn_up <- ggvenn(up_genes, fill_color = c("#D4F6FF","#FFF1D4","#FFD5F9"), show_elements = TRUE, label_sep = "\n")

final_venn <- ggarrange(venn_up, venn_down, nrow = 1)
ggsave("Venn diagram.jpg", units ="in", width = 15, height = 15, dpi = 600)
