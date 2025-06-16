library(openxlsx)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)

# Reading MS data with >=5 PSMs and list of PEG10 interactors and EV proteins
MSdata <- read.xlsx(file.choose())
MSdata <- MSdata[MSdata$`#.PSMs`>=5,] %>% distinct(Gene.Symbol, .keep_all = TRUE)
PEG10_interactors <- read.xlsx(file.choose())
SG_proteins <- read.csv(file.choose())

# Creating a blank table with the required column names to capture PEG10 interactors in MS data set as well as corresponding log2FC values
columns <- c("Gene.Symbol", "RTL8 KO", "PEG10 KO", "UBQLN2 KO","Font")
heat_map_list_PEG10_interactors <- data.frame(matrix(ncol = length(columns), nrow = 126))
colnames(heat_map_list_PEG10_interactors) <- columns

# Copy gene names of PEG10 interactors that are >=5 PSMs in MS data set
heat_map_list_PEG10_interactors$Gene.Symbol <- if_else(PEG10_interactors$`Gene.Name` %in% MSdata$Gene.Symbol, PEG10_interactors$`Gene.Name`, NA)
heat_map_list_PEG10_interactors <- heat_map_list_PEG10_interactors[!is.na(heat_map_list_PEG10_interactors$Gene.Symbol),]

# Join heat map and MS datasets, copy values from columns in MS data set to respective columns in heat map data set and only keep selected columns
heat_map_list_PEG10_interactors <- heat_map_list_PEG10_interactors %>% left_join(MSdata, by = "Gene.Symbol") %>%
  mutate(`RTL8 KO` = coalesce(`RTL8 KO`, `Abundance.Ratio.(log2):.(RTL8.KO.Arsenite)./.(NT.Arsenite)`)) %>%
  mutate(`PEG10 KO` = coalesce(`PEG10 KO`, `Abundance.Ratio.(log2):.(PEG10.KO.Arsenite)./.(NT.Arsenite)`)) %>%
  mutate(`UBQLN2 KO` = coalesce(`UBQLN2 KO`, `Abundance.Ratio.(log2):.(UBQLN2.KO.arsenite)./.(NT.Arsenite)`)) %>%
  select(Gene.Symbol, `RTL8 KO`, `PEG10 KO`, `UBQLN2 KO`)

# Mark SG proteins with bold font
heat_map_list_PEG10_interactors$Font <- if_else(heat_map_list_PEG10_interactors$Gene.Symbol %in% SG_proteins$Protein, "bold","plain")

# Creating a matrix for heatmap() function keeping values from cell lines and converting gene names to rows
heat_map_matrix_PEG10_interactors <- as.matrix(heat_map_list_PEG10_interactors[,c(2:4)])
rownames(heat_map_matrix_PEG10_interactors) <- as.character(heat_map_list_PEG10_interactors$Gene.Symbol)

# Exporting as jpeg with 10x10 in size and 600 dpi resolution
jpeg("PEG10 interactors heatmap.jpg", width = 6000, height = 6000, res = 600)

# Creating a heatmap with custom parameters
heat_map_PEG10 <- Heatmap(heat_map_matrix_PEG10_interactors,
        col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")), # Color gradient from red to white to blue
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        border_gp = gpar(col = "black", lty = 2),
        row_dend_width = unit(0.75, "cm"),  # Adjust the width of row and column dendrogram
        column_dend_height = unit(0.375, "cm"), 
        row_title = expression(paste(bolditalic("PEG10 interactors"), " (", paste(bold("Bold")), "- SG proteins)")),
        row_title_side = "left",
        row_title_gp = gpar(fontsize = 16),
        row_names_gp = gpar(fontsize = 12, fontface = heat_map_list_PEG10_interactors$Font),
        row_gap = unit(5,"mm"),
        column_title = "Cell line",
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_rot = 45,
        heatmap_legend_param = list(
          title = expression(log[2]~(KO/NTC)~"in SG IP"),
          title_position = "topcenter", 
          title_gp = gpar(fontsize = 12, fontface = "plain"),
          legend_position = "bottom",
          legend_direction = "horizontal",
          labels_gp = gpar(fontsize = 11, fontface = "plain")
        )
)

draw(heat_map_PEG10, heatmap_legend_side = "top")

dev.off() # Turn off export

# Reading list of EV proteins
EV_proteins <- read.xlsx(file.choose())

# Creating a blank table with the required column names to capture EV proteins in MS data set as well as corresponding log2FC values
columns <- c("Gene.Symbol", "RTL8 KO", "PEG10 KO", "UBQLN2 KO")
 heat_map_list_EV_proteins <- data.frame(matrix(ncol = length(columns), nrow = 30))
colnames(heat_map_list_EV_proteins) <- columns

# Copy gene names of EV proteins that are >=5 PSMs in MS data set
 heat_map_list_EV_proteins$Gene.Symbol <- if_else(EV_proteins$`Gene.Name` %in% MSdata$Gene.Symbol, EV_proteins$`Gene.Name`, NA)
 heat_map_list_EV_proteins <-  heat_map_list_EV_proteins[!is.na(heat_map_list_EV_proteins$Gene.Symbol),]
 heat_map_list_EV_proteins <- heat_map_list_EV_proteins %>% add_row(Gene.Symbol = "ATXN10") %>% add_row(Gene.Symbol = "TCAF1")

# Join heat map and MS datasets, copy values from columns in MS data set to respective columns in heat map data set and only keep selected columns
 heat_map_list_EV_proteins <- heat_map_list_EV_proteins %>% left_join(MSdata, by = "Gene.Symbol") %>%
  mutate(`RTL8 KO` = coalesce(`RTL8 KO`, `Abundance.Ratio.(log2):.(RTL8.KO.Arsenite)./.(NT.Arsenite)`)) %>%
  mutate(`PEG10 KO` = coalesce(`PEG10 KO`, `Abundance.Ratio.(log2):.(PEG10.KO.Arsenite)./.(NT.Arsenite)`)) %>%
  mutate(`UBQLN2 KO` = coalesce(`UBQLN2 KO`, `Abundance.Ratio.(log2):.(UBQLN2.KO.arsenite)./.(NT.Arsenite)`)) %>%
  select(Gene.Symbol, `RTL8 KO`, `PEG10 KO`, `UBQLN2 KO`)
 
# Mark SG proteins with bold font
 heat_map_list_EV_proteins$Font <- if_else(heat_map_list_EV_proteins$Gene.Symbol %in% SG_proteins$Protein, "bold","plain")

# Creating a matrix for heatmap() function keeping values from cell lines and converting gene names to rows
heat_map_matrix_EV_proteins <- as.matrix(heat_map_list_EV_proteins[,c(2:4)])
rownames(heat_map_matrix_EV_proteins) <-  heat_map_list_EV_proteins$Gene.Symbol

# Exporting as jpeg with 6x6 in size and 600 dpi resolution
jpeg("EV proteins heatmap.jpg", width = 3600, height = 3600, res = 600)

# Creating a heatmap with custom parameters
heat_map_EVs <- Heatmap(heat_map_matrix_EV_proteins,
                        col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")), # Color gradient from red to white to blue
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        border_gp = gpar(col = "black", lty = 2),
                        row_dend_width = unit(0.75, "cm"),  # Adjust the width of row and column dendrogram
                        column_dend_height = unit(0.375, "cm"), 
                        row_title = expression(paste(bolditalic("EV proteins"), " (", paste(bold("Bold")), "- also SG proteins)")),
                        row_title_side = "left",
                        row_title_gp = gpar(fontsize = 16),
                        row_names_gp = gpar(fontsize = 12, fontface = heat_map_list_EV_proteins$Font),
                        row_gap = unit(1,"mm"),
                        column_title = "Cell line",
                        column_title_side = "bottom",
                        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                        column_names_rot = 45,
                        heatmap_legend_param = list(
                          title = expression(log[2]~(KO/NTC)~"in SG IP"),
                          title_position = "topcenter", 
                          title_gp = gpar(fontsize = 12, fontface = "plain"),
                          legend_position = "bottom",
                          legend_direction = "horizontal",
                          labels_gp = gpar(fontsize = 11, fontface = "plain")
                        )
                      )

draw(heat_map_EVs, heatmap_legend_side = "top")

dev.off() # Turn off export
