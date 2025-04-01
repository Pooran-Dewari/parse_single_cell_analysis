
# load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)
library(dittoSeq)
library(ggvenn)


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



##### step 2: find DE genes, append to global lists, and plot venn diagram  ####

# Create a parent directory to store results
parent_dir <- "venn_plots"
dir.create(parent_dir, showWarnings = FALSE)

# Initialize global lists for upregulated and downregulated genes
all_up_genes_by_condition <<- list(early = list(), mid = list(), late = list())
all_down_genes_by_condition <<- list(early = list(), mid = list(), late = list())

# Function to perform DE analysis for each condition
perform_comparison <- function(condition) {
  control_clusters <- unique(seu_obj_clean$sample_subtypes[grepl("control", seu_obj_clean$sample_subtypes)])
  target_clusters <- unique(seu_obj_clean$sample_subtypes[grepl(condition, seu_obj_clean$sample_subtypes)])
  
  # Ensure matching cluster order
  control_clusters <- control_clusters[order(gsub("_control", "", control_clusters))]
  target_clusters <- target_clusters[order(gsub(paste0("_", condition), "", target_clusters))]
  
  for (i in seq_along(control_clusters)) {
    control_label <- control_clusters[i]
    target_label <- target_clusters[i]
    tissue_name <- gsub("_control", "", control_label)
    
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
    valid_up_genes <- valid_up_genes[valid_up_genes != "" & !is.na(valid_up_genes)]
    valid_down_genes <- valid_down_genes[valid_down_genes != "" & !is.na(valid_down_genes)]
    
    if (length(valid_up_genes) > 0) 
      all_up_genes_by_condition[[condition]][[tissue_name]] <<- valid_up_genes
    
    if (length(valid_down_genes) > 0) 
      all_down_genes_by_condition[[condition]][[tissue_name]] <<- valid_down_genes
  }
}

# Run function for different conditions
lapply(c("early", "mid", "late"), perform_comparison)

# Create directories only if genes exist
venn_up_dir <- file.path(parent_dir, "venn_up")
venn_down_dir <- file.path(parent_dir, "venn_down")

if (any(sapply(all_up_genes_by_condition, function(x) length(unlist(x)) > 0))) 
  dir.create(venn_up_dir, showWarnings = FALSE, recursive = TRUE)

if (any(sapply(all_down_genes_by_condition, function(x) length(unlist(x)) > 0))) 
  dir.create(venn_down_dir, showWarnings = FALSE, recursive = TRUE)

# Function to generate Venn diagrams
generate_venn <- function(gene_list, condition, type, output_dir) {
  for (cluster_name in names(gene_list[[condition]])) {
    genes <- gene_list[[condition]][[cluster_name]]
    
    if (length(genes) > 0) {
      list_genes <- list(
        early = all_up_genes_by_condition$early[[cluster_name]],
        mid = all_up_genes_by_condition$mid[[cluster_name]],
        late = all_up_genes_by_condition$late[[cluster_name]]
      )
      list_genes <- lapply(list_genes, function(x) x[x != "" & !is.na(x)])
      
      if (all(sapply(list_genes, length) > 0)) {
        venn_plot <- ggvenn(list_genes, 
                            fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
                            stroke_size = 0.6, 
                            set_name_size = 6) +
          ggtitle(paste(type, "genes in", cluster_name))+
          theme(plot.title = element_text(size = 20,  face = "bold", hjust = 0.6))  
        
        venn_file <- file.path(output_dir, paste0(cluster_name, "_venn_", type, ".pdf"))
        ggsave(venn_file, plot = venn_plot,  width = 8, height = 8, units = "in")
      }
    }
  }
}

# Generate Venn diagrams
lapply(c("early", "mid", "late"), function(cond) generate_venn(all_up_genes_by_condition, cond, "up", venn_up_dir))
lapply(c("early", "mid", "late"), function(cond) generate_venn(all_down_genes_by_condition, cond, "down", venn_down_dir))
