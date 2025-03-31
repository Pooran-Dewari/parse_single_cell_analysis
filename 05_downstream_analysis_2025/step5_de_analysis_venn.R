library(ggvenn)
library(Seurat)
library(dplyr)

# Create a parent directory to store all results
parent_dir <- "venn_plots"
dir.create(parent_dir, showWarnings = FALSE)

# Initialize the global list structure for upregulated and downregulated genes
all_up_genes_by_condition <<- list(
  early = list(),
  mid = list(),
  late = list()
)

all_down_genes_by_condition <<- list(
  early = list(),
  mid = list(),
  late = list()
)

# Define function to perform DE analysis for each condition
perform_comparison <- function(condition) {
  # Create a subdirectory for the specific condition inside the parent directory
  #condition_dir <- file.path(parent_dir, condition)
  #dir.create(condition_dir, showWarnings = FALSE)
  
  control_clusters <- unique(seu_obj_clean$sample_subtypes[grepl("control", seu_obj_clean$sample_subtypes)])
  target_clusters <- unique(seu_obj_clean$sample_subtypes[grepl(condition, seu_obj_clean$sample_subtypes)])
  
  # Ensure matching cluster order
  control_clusters <- control_clusters[order(gsub("_control", "", control_clusters))]
  target_clusters <- target_clusters[order(gsub(paste0("_", condition), "", target_clusters))]
  
  # Print control and target clusters to check their values
  print(paste("Control clusters for", condition, ":"))
  print(control_clusters)
  print(paste("Target clusters for", condition, ":"))
  print(target_clusters)
  
  # Loop through each pair of control and target clusters
  for (i in 1:length(control_clusters)) {
    control_label <- control_clusters[i]
    target_label <- target_clusters[i]
    
    # Extract the readable tissue/type name from control_label
    tissue_name <- gsub("_control", "", control_label)  # e.g., "adductor_muscle"
    
    #comparison_dir <- file.path(condition_dir, paste0(tissue_name, "_control_vs_", condition))
    #dir.create(comparison_dir, showWarnings = FALSE)
    
    # Perform differential expression analysis
    cluster_control_target <- FindMarkers(seu_obj_clean, logfc.threshold = 1,
                                          min.pct = 0.25,
                                          ident.1 = control_label,
                                          ident.2 = target_label) %>%
      filter(p_val_adj < 0.01) %>%
      rownames_to_column(var = "gene")
    
    # Print the number of DE genes found
    print(paste("Found", nrow(cluster_control_target), "differentially expressed genes between", control_label, "and", target_label))
    
    # Get upregulated and downregulated genes
    valid_up_genes <- cluster_control_target %>% filter(avg_log2FC < -1) %>% pull(gene)
    valid_down_genes <- cluster_control_target %>% filter(avg_log2FC > 1) %>% pull(gene)
    
    # Remove NA and empty values
    valid_up_genes <- valid_up_genes[!is.na(valid_up_genes) & valid_up_genes != ""]
    valid_down_genes <- valid_down_genes[!is.na(valid_down_genes) & valid_down_genes != ""]
    
    # Add upregulated genes to the global list for the current condition and cluster
    if (length(valid_up_genes) > 0) {
      print(paste("Adding", length(valid_up_genes), "upregulated genes for", tissue_name, condition, "to the list"))
      all_up_genes_by_condition[[condition]][[tissue_name]] <<- valid_up_genes
    }
    
    # Add downregulated genes to the global list for the current condition and cluster
    if (length(valid_down_genes) > 0) {
      print(paste("Adding", length(valid_down_genes), "downregulated genes for", tissue_name, condition, "to the list"))
      all_down_genes_by_condition[[condition]][[tissue_name]] <<- valid_down_genes
    }
  }
}

# Run function for different comparisons (early, mid, late)
comparison_conditions <- c("early", "mid", "late")
lapply(comparison_conditions, perform_comparison)

# Create necessary directories for storing the Venn diagrams
venn_up_dir <- file.path(parent_dir, "venn_diagram", "venn_up")
venn_down_dir <- file.path(parent_dir, "venn_diagram", "venn_down")
dir.create(venn_up_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(venn_down_dir, showWarnings = FALSE, recursive = TRUE)

# Generate Venn diagrams for upregulated genes
for (condition in c("early", "mid", "late")) {
  for (cluster_name in names(all_up_genes_by_condition[[condition]])) {
    up_genes <- all_up_genes_by_condition[[condition]][[cluster_name]]
    
    if (length(up_genes) > 0) {
      # Filter out empty clusters and create Venn diagrams for each cluster's upregulated genes
      list_up <- list(
        early_up = all_up_genes_by_condition$early[[cluster_name]],
        mid_up = all_up_genes_by_condition$mid[[cluster_name]],
        late_up = all_up_genes_by_condition$late[[cluster_name]]
      )
      list_up <- lapply(list_up, function(x) x[!is.na(x) & x != ""])  # Filter out NA and empty values
      
      # Only create Venn diagram if the list is not empty
      if (length(list_up$early_up) > 0 && length(list_up$mid_up) > 0 && length(list_up$late_up) > 0) {
        # Create the Venn diagram plot
        venn_up_plot <- ggvenn(list_up, fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
                               stroke_size = 0.5, set_name_size = 4) +
          ggtitle(paste("Upregulated Genes -", cluster_name))
        
        # Save the Venn diagram plot
        venn_up_file <- file.path(venn_up_dir, paste0(cluster_name, "_venn_up.pdf"))
        dir.create(dirname(venn_up_file), showWarnings = FALSE)
        
        # Save the plot as a PDF
        ggsave(venn_up_file, plot = venn_up_plot)
      }
    }
  }
}

# Generate Venn diagrams for downregulated genes
for (condition in c("early", "mid", "late")) {
  for (cluster_name in names(all_down_genes_by_condition[[condition]])) {
    down_genes <- all_down_genes_by_condition[[condition]][[cluster_name]]
    
    if (length(down_genes) > 0) {
      # Filter out empty clusters and create Venn diagrams for each cluster's downregulated genes
      list_down <- list(
        early_down = all_down_genes_by_condition$early[[cluster_name]],
        mid_down = all_down_genes_by_condition$mid[[cluster_name]],
        late_down = all_down_genes_by_condition$late[[cluster_name]]
      )
      list_down <- lapply(list_down, function(x) x[!is.na(x) & x != ""])  # Filter out NA and empty values
      
      # Only create Venn diagram if the list is not empty
      if (length(list_down$early_down) > 0 && length(list_down$mid_down) > 0 && length(list_down$late_down) > 0) {
        # Create the Venn diagram plot
        venn_down_plot <- ggvenn(list_down, fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
                                 stroke_size = 0.5, set_name_size = 4) +
          ggtitle(paste("Downregulated Genes -", cluster_name))
        
        # Save the Venn diagram plot
        venn_down_file <- file.path(venn_down_dir, paste0(cluster_name, "_venn_down.pdf"))
        dir.create(dirname(venn_down_file), showWarnings = FALSE)
        
        # Save the plot as a PDF
        ggsave(venn_down_file, plot = venn_down_plot)
      }
    }
  }
}


# manually test one if need be
cluster_name <- "cluster1"  # Replace with your specific cluster if needed

# Get upregulated genes for 'cluster_1' across all conditions
list_up <- list(
  early_up = all_up_genes_by_condition$early[[cluster_name]],
  mid_up = all_up_genes_by_condition$mid[[cluster_name]],
  late_up = all_up_genes_by_condition$late[[cluster_name]]
)

list_up

dev.off()
dev.off()

ggvenn(list_up, fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4) +
  ggtitle(paste("Upregulated Genes -", cluster_name))
