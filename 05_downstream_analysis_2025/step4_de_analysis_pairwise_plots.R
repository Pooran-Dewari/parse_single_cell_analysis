# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dittoSeq)

# Create a parent directory to store all results
parent_dir <- "de_plots"
dir.create(parent_dir, showWarnings = FALSE)

# Define function to perform DE analysis for each condition
perform_comparison <- function(condition) {
  # Create a subdirectory for the specific condition inside the parent directory
  condition_dir <- file.path(parent_dir, condition)
  dir.create(condition_dir, showWarnings = FALSE)
  
  control_clusters <- unique(seu_obj_clean$sample_subtypes[grepl("control", seu_obj_clean$sample_subtypes)])
  target_clusters <- unique(seu_obj_clean$sample_subtypes[grepl(condition, seu_obj_clean$sample_subtypes)])
  
  # Ensure matching cluster order
  control_clusters <- control_clusters[order(gsub("_control", "", control_clusters))]
  target_clusters <- target_clusters[order(gsub(paste0("_", condition), "", target_clusters))]
  
  for (i in 1:length(control_clusters)) {
    control_label <- control_clusters[i]
    target_label <- target_clusters[i]
    
    # Extract the readable tissue/type name from control_label
    tissue_name <- gsub("_control", "", control_label)  # e.g., "adductor_muscle"
    
    comparison_dir <- file.path(condition_dir, paste0(tissue_name, "_control_vs_", condition))
    dir.create(comparison_dir, showWarnings = FALSE)
    
    # Perform differential expression analysis
    cluster_control_target <- FindMarkers(seu_obj_clean, logfc.threshold = 1,
                                          min.pct = 0.25,
                                          ident.1 = control_label,
                                          ident.2 = target_label) %>%
      filter(p_val_adj < 0.01) %>%
      rownames_to_column(var = "gene")
    
    # Get upregulated and downregulated genes
    valid_up_genes <- cluster_control_target %>% filter(avg_log2FC < -1) %>% pull(gene)
    valid_down_genes <- cluster_control_target %>% filter(avg_log2FC > 1) %>% pull(gene)
    
    # Remove NA and empty values
    valid_up_genes <- valid_up_genes[!is.na(valid_up_genes) & valid_up_genes != ""]
    valid_down_genes <- valid_down_genes[!is.na(valid_down_genes) & valid_down_genes != ""]
    
    print(paste("Valid upregulated genes for", tissue_name, "control vs", condition, ":"))
    print(valid_up_genes)
    
    print(paste("Valid downregulated genes for", tissue_name, "control vs", condition, ":"))
    print(valid_down_genes)
    
    # Function to generate plots with readable filenames
    plot_genes <- function(valid_genes, plot_type, suffix) {
      if (length(valid_genes) > 2) {
        valid_genes_in_seurat <- valid_genes[valid_genes %in% rownames(seu_obj_clean)]
        
        if (length(valid_genes_in_seurat) > 0) {
          print(paste("Plotting", plot_type, "genes for", tissue_name, "control vs", condition, ":"))
          print(valid_genes_in_seurat)
          
          file_name <- paste0(tissue_name, "_control_vs_", condition, "_de_genes_", suffix, ".pdf")
          pdf(file.path(comparison_dir, file_name))
          
          if (plot_type == "VlnPlot") {
            print(VlnPlot(seu_obj_clean, features = valid_genes_in_seurat[1:min(10, length(valid_genes_in_seurat))], 
                          stack = TRUE, flip = TRUE, group.by = "condition_new"))
          } else {
            print(DotPlot(seu_obj_clean, features = rev(valid_genes_in_seurat[1:min(10, length(valid_genes_in_seurat))]), 
                          group.by = "condition_new") +
                    coord_flip() +
                    scale_x_discrete(position = "top"))
          }
          
          dev.off()
        } else {
          print(paste("No valid", plot_type, "genes found in Seurat object for", tissue_name, "control vs", condition))
        }
      } else {
        print(paste("Not enough", plot_type, "genes found for", tissue_name, "control vs", condition, "(more than 2 required)"))
      }
    }
    
    # Generate plots
    plot_genes(valid_up_genes, "VlnPlot", "vlnplot_up")
    plot_genes(valid_down_genes, "VlnPlot", "vlnplot_down")
    plot_genes(valid_up_genes, "DotPlot", "dotplot_up")
    plot_genes(valid_down_genes, "DotPlot", "dotplot_down")
    
    # Generate heatmap
    cluster_aggregate_exp <- AggregateExpression(seu_obj_clean, return.seurat = TRUE, group.by = "condition_new")
    
    heatmap_plot <- dittoHeatmap(cluster_aggregate_exp, genes = cluster_control_target$gene,
                                 annot.by = "condition_new",
                                 heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
                                 main = paste(tissue_name, "control vs", condition, "differentially expressed genes"),
                                 cluster_cols = TRUE, cluster_rows = TRUE, scale = "row", show_colnames = TRUE)
    
    heatmap_file <- paste0(tissue_name, "_control_vs_", condition, "_de_genes_heatmap.pdf")
    
    pdf(file.path(comparison_dir, heatmap_file))
    
    print(heatmap_plot)
    
    dev.off()
    
    # Save gene lists with readable filenames
    base_filename <- paste0(tissue_name, "_control_vs_", condition, "_de_genes_")
    
    write.table(cluster_control_target$gene, file.path(comparison_dir, paste0(base_filename, "all.txt")), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(valid_up_genes, file.path(comparison_dir, paste0(base_filename, "up.txt")), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(valid_down_genes, file.path(comparison_dir, paste0(base_filename, "down.txt")), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Run function for different comparisons
comparison_conditions <- c("late", "mid", "early")
lapply(comparison_conditions, perform_comparison)
