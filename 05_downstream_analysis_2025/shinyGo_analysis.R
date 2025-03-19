# GO enrichment and KEGG pathway analysis is a pain with the (slightly) older Roslin genome assembly of Pacific oyster, solely because of poor annotation
# luckily, ShinyGo webtool has GO annotation for this version of the genome, only downside is can't do batch-processing, copy paste one set at a time!
# another option is to use g:Profiler, but that seems a bit too stringent
# I will go ahead with shinyGO


############ Step1: find all marker genes #############################################

# Find markers genes and annoatate using ORSON supl data from the French group
# All we need is gene ID but joining with published ORSON/cg_science datasets to retrieve gene names/symbols

all_markers_test_seu <- FindAllMarkers(object = test_seu, min.pct = 0.20, log2fc.threshold = 0.25)

# Load ORSON annotation
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv")

# Load cg_science paper annotation
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  as.data.frame()

colnames(cg_science)[c(1, 2)] <- c("Gene.ID", "Gene.symbol.science")

# Annotate using both annotations above
all_markers_test_seu_ann <- left_join(all_markers_test_seu, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

################ step 2: Get background gene list ###########################################################
# background genes for shinygo analysis, this is the list of ALL genes detected in our single-nucleus data
all_genes_detected <- as.data.frame(rownames(test_seu))
write_tsv(all_genes_detected, "all_genes_detected.txt")

################ step 3: Get upregulated genes for each cluster      #########################################
all_up <- all_markers_test_seu_ann %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.10)

# Save genes for each cluster as tsv files for shinygo
all_up_clean <- all_up %>% 
  dplyr::select(cluster, gene)

# Create a new directory if it doesn't exist
dir.create("signif_genes_shinygo", showWarnings = FALSE)

# Get unique cluster values
unique_clusters <- unique(all_up_clean$cluster)

# Iterate through each cluster and save as a .tsv file inside "signif_genes_shinygo"
for (cluster_id in unique_clusters) {
  subset_data <- filter(all_up_clean, cluster == cluster_id)  # Subset data
  filename <- paste0("signif_genes_shinygo/up_cluster_", cluster_id, ".tsv")  # Create path & filename
  write.table(subset_data, file = filename, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

################ Step 4: Get downregulated genes for each cluster      #########################################
all_down <- all_markers_test_seu_ann %>% 
  filter(avg_log2FC < -0.50 & p_val_adj < 0.10)

all_down_clean <- all_down %>% 
  dplyr::select(cluster, gene)

# Create a new directory if it doesn't exist
dir.create("signif_genes_shinygo", showWarnings = FALSE)

# Get unique cluster values
unique_clusters <- unique(all_down_clean$cluster)

# Iterate through each cluster and save as a .tsv file inside "signif_genes_shinygo"
for (cluster_id in unique_clusters) {
  subset_data <- filter(all_down_clean, cluster == cluster_id)  # Subset data
  filename <- paste0("signif_genes_shinygo/down_cluster_", cluster_id, ".tsv")  # Create path & filename
  write.table(subset_data, file = filename, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

################# Step 5: Run shinyGO web tool ####################################################
# https://bioinformatics.sdstate.edu/go/
# ShinyGO 0.82
# Select species: Crassostrea gigas Pacific oyster - GCA_902806645.1
# In the input gene list, e.g. paste all up genes for cluster 0, you can query only one gene list at a time
# In the background gene list, paste background genes from step 2 above
# Submit job and select 'GO Biological Process' under Pathway database menu, select 30 pathways to show, check 'Remove redundancy' and 'Show pathway IDs'; uncheck 'Abbreviate pathways'
# Sort results by 'Select by FDR, sort by Fold Enrichment'
# Download 'Top Pathways shown above' in shinygo directory, we are interested in GO BP only.
# Fold enrichment is defined as the percentage of genes in your list that are in a pathway divided by the corresponding percentage in the background genes
# Fold Enrichment: Measures the magnitude of enrichment. Higher values indicate stronger enrichment and are an important metric of effect size
# FDR q-values adjust P-values for multiple testing to control the proportion of type I errors

###################### Step 6: Plot the dot plot using ggplot2 ###############################

#.................... 6a: manual one plot at a time
# Load up and down genes shinyGO results for cluster 7
go_up <- read_csv("shinygo/up_cl_7.csv") %>% 
  mutate(Condition = "up_genes")

go_down <-  read_csv("shinygo/down_cl_7.csv") %>% 
  mutate(Condition = "down_genes")

# Combine into one
go_combined <- bind_rows(go_up, go_down) %>% 
  mutate(Condition = factor(Condition, levels = c("up_genes", "down_genes")))

# Plot the dot plot
ggplot(data = go_combined, aes(x = Condition, y = Pathway, 
                        color = `Enrichment FDR`, size = `Fold Enrichment`)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("Cluster 7") + 
  ggtitle("GO Biological Process enrichment")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#.................... 6b: batch processing
# Requires shinyGO csv files for up and down in shinygo directory
# You need to amend number of clusters in the script below

library(tidyverse)

# Ensure the output directory exists
output_dir <- "shinygo/shiny_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Extract unique cluster numbers from filenames
clusters <- 0:20  # Update this if you have more clusters

# Loop over each cluster and generate plots
for (cl in clusters) {
  # Read data for the current cluster
  go_up <- read_csv(paste0("shinygo/up_cluster_", cl, ".csv")) %>%
    mutate(Condition = "up_genes")
  
  go_down <- read_csv(paste0("shinygo/down_cluster_", cl, ".csv")) %>%
    mutate(Condition = "down_genes")
  
  # Combine data
  go_combined <- bind_rows(go_up, go_down) %>%
    mutate(Condition = factor(Condition, levels = c("up_genes", "down_genes")))
  
  # Plot
  plot <- ggplot(data = go_combined, aes(x = Condition, y = Pathway, 
                                         color = `Enrichment FDR`, size = `Fold Enrichment`)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + 
    ylab("") + 
    xlab(paste("Cluster", cl)) + 
    ggtitle(paste("GO Biological Process Enrichment - Cluster", cl)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.y=element_text(size=8), axis.text.x=element_text(size=8))
  
  # Save plot in shinygo/shiny_plots
  ggsave(filename = paste0(output_dir, "/GO_Enrichment_Cluster_", cl, ".png"), 
         plot = plot, width = 10, height = 12)
}


