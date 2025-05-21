# load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)
library(dittoSeq)
library(Seurat.utils)


##### step 1: load Seurat object & prepare columns for group comparison #######

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

# load Seurat object with pca, harmony, umap done beforehand!
seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")
seu_obj 

new.cluster.ids <- c("cluster0", "cluster1", "gill_ciliary",
                     "hepatopancreas", "gill_nec", "gill_type1",
                     "hyalinocytes", "haemocytes_type1",
                     "cluster8", "cluster9", "mantle_vesicular",
                     "immature_haemocytes", "macrophage_like",
                     "adductor_muscle", "mantle", "digestive_gland",
                     "gill_type2", "small_granules_cells")

names(new.cluster.ids) <- levels(seu_obj)

seu_obj <- RenameIdents(seu_obj, new.cluster.ids) # add cluster info to Idents

#add the annotations to metadata column
#seu_obj$seurat_annotations <- Idents(seu_obj)


#add new condition column to metadata for grouping
#24-hpiA shows bad qPCR for viral dna, so putting "unsure" there
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "early",
  seu_obj$sample %in% "24-hpiA" ~ "mid?",
  seu_obj$sample %in% "24-hpiJ" ~ "mid",
  seu_obj$sample %in% c("72-hpiJ", "96-hpiE") ~ "late",
  TRUE ~ "unknown"  # Optional: Handles unexpected values
)

unique(seu_obj$condition_new)

# reorder levels in the condition_new
seu_obj$condition_new <- factor(seu_obj$condition_new, 
                                levels = c("control", "early", "mid?", "mid", "late"))


# subset the Seurat object to exclude 'mid?' , i.e. 24hA
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")

rm(seu_obj)

# drop unused levels
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

# create a column in meta data to hold cluster info & condition
# we will use info from this column to compare groups
seu_obj_clean$sample_subtypes <- paste(Idents(seu_obj_clean), seu_obj_clean$condition_new, sep = "_")
Idents(seu_obj_clean) <- "sample_subtypes" # seurat object now has new Idents = old idents + condition_new
unique(seu_obj_clean$sample_subtypes)


### optional: we could replace gene IDs with gene symbols for readability in heatmap plots..#####

# load ORSON annotation
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv") %>% 
  select(Gene.ID, Description = Sequence.Description)

#load cg science annotation
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  rename(Gene.ID = "X1", Gene.symbol.science = "X2" )

# Perform a full join to keep all rows
merged_df <- full_join(ORSON, cg_science, by = "Gene.ID")


# Perform a full join and concatenate Gene.symbol.science first
merged_df <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""), 
      coalesce(Description, ""), 
      sep = " : "
    ) %>% trimws()  # Remove extra spaces
  ) %>% 
  mutate(
    gene_symbol = gsub("[^[:alnum:] :]", " ", gene_symbol)  # Keep "-" but replace other special characters
  )

gene_map2 <- merged_df %>%
  select(gene_id = Gene.ID, gene_symbol)



# duplicates are a PAIN and will generate error below
#Error in .which_data(assay, slot, object)[genes, cells.use] : 
#subscript out of bounds
#Called from: as.matrix(.which_data(assay, slot, object)[genes, cells.use])

#Custom suffix for duplicates
gene_map2 <- gene_map2 %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_symbol = if(n() > 1) {
      paste0(gene_symbol, LETTERS[seq_along(gene_symbol)])
    } else {
      gene_symbol
    }
  ) %>%
  ungroup()


# add gene ID to gene symbol so that we can plot later using gene id
gene_map2 <- gene_map2 %>%
  mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " "))

# Filter the mapping to include only genes present in Seurat object
gene_map_filtered <- gene_map2 %>% 
  filter(gene_id %in% rownames(seu_obj_clean))

# Create a named vector for renaming
gene_symbols <- gene_map_filtered$gene_symbol[match(rownames(seu_obj_clean), gene_map_filtered$gene_id)]
rename_vector <- setNames(gene_symbols, rownames(seu_obj_clean))

# Handle missing gene symbols by keeping their original IDs
rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]

# Ensure uniqueness of gene symbols
rename_vector <- make.unique(rename_vector)

# Rename genes in Seurat object
seu_obj_clean <- RenameGenesSeurat(
  obj = seu_obj_clean,
  newnames = rename_vector
)

rownames(seu_obj_clean)
################# step 2: find pairwise DE genes across conditions & plot ######

# Create a parent directory to store all results
parent_dir <- "de_plots_full_ann_v8"
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
    
    # Calculate progress percentage
    progress <- round((i / length(control_clusters)) * 100, 1)
    
    # Print progress update
    cat("\rProgress:", progress, "% - Processing", tissue_name, "(control vs", condition, ") ...\n", sep=" ")
    
    
    # Get upregulated and downregulated genes
    valid_up_genes <- cluster_control_target %>% filter(avg_log2FC < -1) %>% pull(gene)
    valid_down_genes <- cluster_control_target %>% filter(avg_log2FC > 1) %>% pull(gene)
    
    # Remove NA and empty values
    valid_up_genes <- valid_up_genes[!is.na(valid_up_genes) & valid_up_genes != ""]
    valid_down_genes <- valid_down_genes[!is.na(valid_down_genes) & valid_down_genes != ""]
    
    #print(paste("Valid upregulated genes for", tissue_name, "control vs", condition, ":"))
    #print(valid_up_genes)
    
    #print(paste("Valid downregulated genes for", tissue_name, "control vs", condition, ":"))
    #print(valid_down_genes)
    
    # Function to generate plots with readable filenames
    plot_genes <- function(valid_genes, plot_type, suffix) {
      if (length(valid_genes) > 2) {
        valid_genes_in_seurat <- valid_genes[valid_genes %in% rownames(seu_obj_clean)]
        
        if (length(valid_genes_in_seurat) > 0) {
          #print(paste("Plotting", plot_type, "genes for", tissue_name, "control vs", condition, ":"))
          #print(valid_genes_in_seurat)
          cat("Plotting", length(valid_genes_in_seurat), plot_type, "genes for", tissue_name, "(control vs", condition, ")...\n")
          
          file_name <- paste0(tissue_name, "_control_vs_", condition, "_de_genes_", suffix, ".pdf")
          #pdf(file.path(comparison_dir, file_name))
          pdf(file.path(comparison_dir, file_name), width = 12, height = 10)  # **Increased figure size**
          
          plot_title <- paste0("Top 10 markers: ", tissue_name, " (control vs ", condition, ")")
          
          if (plot_type == "VlnPlot") {
            print(VlnPlot(seu_obj_clean, features = valid_genes_in_seurat[1:min(10, length(valid_genes_in_seurat))], 
                          stack = TRUE, flip = TRUE, group.by = "condition_new") + ggtitle(plot_title) + NoLegend()) 
          } else {
            print(DotPlot(seu_obj_clean, features = rev(valid_genes_in_seurat[1:min(10, length(valid_genes_in_seurat))]), 
                          group.by = "condition_new") +
                    coord_flip() +
                    scale_x_discrete(position = "top") + ggtitle(plot_title))
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
    
    #..............save the matrix for use later.................
    
    # Save the aggregated expression matrix used for the heatmap
    expr_mat <- GetAssayData(cluster_aggregate_exp, slot = "data")
    
    # Subset to only genes used in heatmap (those that passed DE filters)
    expr_mat_subset <- expr_mat[cluster_control_target$gene, ]
    
    # Convert to data frame
    expr_df <- as.data.frame(as.matrix(expr_mat_subset))
    expr_df <- tibble::rownames_to_column(expr_df, var = "gene")
    
    # Prefix condition columns with comparison info (e.g., "gill_nec_vs_early.control")
    comparison_label <- paste0(tissue_name, "_vs_", condition)
    colnames(expr_df)[-1] <- paste0(comparison_label, ".", colnames(expr_df)[-1])
    
    # Save to TSV
    expr_table_file <- paste0(comparison_label, "_heatmap_expr_table.tsv")
    write_tsv(expr_df, file.path(comparison_dir, expr_table_file))
    
    
    #..............save the matrix for use later.................
    
    heatmap_plot <- dittoHeatmap(cluster_aggregate_exp, genes = cluster_control_target$gene,
                                 annot.by = "condition_new",
                                 heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
                                 main = paste0("Total ", length(cluster_control_target$gene), " DE genes: ", 
                                               tissue_name, " (control vs ", condition, ")"),
                                 cluster_cols = TRUE, cluster_rows = TRUE, scale = "row", show_colnames = TRUE)
    
    heatmap_file <- paste0(tissue_name, "_control_vs_", condition, "_de_genes_heatmap.pdf")
    
    pdf(file.path(comparison_dir, heatmap_file), width = 12, height = 10)  # **Increased figure size**
    
    #pdf(file.path(comparison_dir, heatmap_file))
    
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

####################### ends #################################################


# code below will create a unified list of DE genes for each cluster by combining...
# DE genes from early, mid, and late; and create a single heatmap

# Ensure necessary directories exist
unified_dir <- file.path(parent_dir, "unified_gene_lists")
dir.create(unified_dir, showWarnings = FALSE)

heatmap_dir <- file.path(unified_dir, "heatmap_plots")
dir.create(heatmap_dir, showWarnings = FALSE)

# Generate aggregated expression data once
cluster_aggregate_exp <- AggregateExpression(seu_obj_clean, return.seurat = TRUE, group.by = "condition_new")
expr_mat <- GetAssayData(cluster_aggregate_exp, slot = "data")

# Update condition_new with alphabetical prefixes
seu_obj_clean$condition_new <- factor(
  paste0(
    c("A-", "B-", "C-", "D-")[match(seu_obj_clean$condition_new, c("control", "early", "mid", "late"))],
    seu_obj_clean$condition_new
  ),
  levels = c("A-control", "B-early", "C-mid", "D-late")
)


# Loop through each cluster
processed_clusters <- unique(sub("_vs_.*", "", names(cluster_gene_list)))

for (base_name in processed_clusters) {
  # Remove '_control' suffix if present
  base_name_clean <- sub("_control$", "", base_name)
  
  combined_genes <- character()
  
  for (condition in c("early", "mid", "late")) {
    full_cluster_name <- paste(base_name, condition, sep = "_vs_")
    if (full_cluster_name %in% names(cluster_gene_list)) {
      combined_genes <- unique(c(combined_genes, cluster_gene_list[[full_cluster_name]]))
    }
  }
  
  # Save unified gene list
  write.table(combined_genes,
              file = file.path(unified_dir, paste0(base_name_clean, "_unified_genes.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Filter for genes that exist in expression data
  valid_genes <- combined_genes[combined_genes %in% rownames(expr_mat)]
  if (length(valid_genes) == 0) {
    warning(paste("No valid genes found in expression matrix for:", base_name_clean))
    next
  }
  
  
  
  # Generate heatmap
  heatmap_plot <- dittoHeatmap(
    object = cluster_aggregate_exp,
    genes = valid_genes,
    annot.by = "condition_new",
    order.by = "condition_new",
    heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste0("Unified DE genes heatmap: ", base_name_clean, " (", length(valid_genes), " genes)"),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    scale = "row",
    show_colnames = TRUE
  )
  
  
  # Save heatmap to PDF
  heatmap_file <- paste0(base_name_clean, "_unified_genes_heatmap.pdf")
  pdf(file.path(heatmap_dir, heatmap_file), width = 12, height = 10)
  print(heatmap_plot)
  dev.off()
}
