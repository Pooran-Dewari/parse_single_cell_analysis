
# script to make main figures: 19/05/2025

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
unique(Idents(seu_obj))
# [1] 16 1  6  2  3  8  12 0  13 5  4  10 9  11 14 15 7  17
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17

new.cluster.ids <- c("cluster0", "cluster1", "gill_ciliary",
                                          "hepatopancreas", "gill_nec", "gill_type1",
                                          "hyalinocytes", "haemocytes_type1",
                                          "cluster8", "cluster9", "mantle_vesicular",
                                          "immature_haemocytes", "macrophage_like",
                                          "adductor_muscle", "mantle", "digestive_gland",
                                          "gill_type2", "small_granules_cells")

names(new.cluster.ids) <- levels(seu_obj)

seu_obj <- RenameIdents(seu_obj, new.cluster.ids) # add cluster info to Idents
unique(Idents(seu_obj))
#[1] gill_type2           cluster1             hyalinocytes   ............

DimPlot(seu_obj, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("Starsolo harmony umap : 18 dim x 6 res x 3k variable features")+
  theme(plot.title = element_text(hjust = 0.5))


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


# drop unused levels
seu_obj_clean$condition_new <- droplevels(seu_obj_clean$condition_new)

rm(seu_obj)

###### 13/05/2025
####### clean seurat object without the problematic 24h-A sample,
####### we will use clean seurat object now onwards for all the downstream analysis



############# Figure 1 starts #######################***************************

#Figure 1A
VlnPlot(seu_obj_clean, features = "percent.mt", group.by = "sample", pt.size = 0.125)+
  NoLegend()+
  ggtitle("Mitochondrial reads contribution")+
  xlab("") +
  ylab("% Mitochondrial reads")+
  theme(
    axis.text.x = element_text(size = 18),      # X-axis tick labels
    axis.text.y = element_text(size = 18),      # Y-axis tick labels
    axis.title.y = element_text(size = 18),     # Y-axis title
    axis.title.x = element_text(size = 18),     # X-axis title
    plot.title = element_text(size = 20, face = "bold")  # Plot title
  )+
  theme(axis.text.x = element_blank())

VlnPlot(seu_obj_clean, features = "nFeature_RNA", group.by = "sample", pt.size = 0.125)+
  NoLegend()+
  ggtitle("Number of genes(features) detected")+
  xlab("") +
  ylab("Number of genes detected")+
  theme(
    axis.text.x = element_text(size = 18),      # X-axis tick labels
    axis.text.y = element_text(size = 18),      # Y-axis tick labels
    axis.title.y = element_text(size = 18),     # Y-axis title
    axis.title.x = element_text(size = 18),     # X-axis title
    plot.title = element_text(size = 20, face = "bold")  # Plot title
  )+
  theme(axis.text.x = element_blank())

VlnPlot(seu_obj_clean, features = "nCount_RNA", group.by = "sample", pt.size = 0.125, y.max = 10000)+
  NoLegend()+
  ggtitle("Number of transcripts (UMIs) detected")+
  xlab("") +
  theme(
    axis.text.x = element_text(size = 18),      # X-axis tick labels
    axis.text.y = element_text(size = 18),      # Y-axis tick labels
    axis.title.y = element_text(size = 18),     # Y-axis title
    axis.title.x = element_text(size = 18),     # X-axis title
    plot.title = element_text(size = 20, face = "bold")  # Plot title
  )+
  ylab("Number of transcripts detected")

#median UMI per sample
seu_obj_clean@meta.data %>%
  group_by(sample) %>%
  summarize(median_UMI = median(nCount_RNA))

#median UMI all dataset
seu_obj_clean@meta.data %>%
  summarize(median_UMI = median(nCount_RNA))


#median genes per sample
seu_obj_clean@meta.data %>%
  group_by(sample) %>%
  summarize(median_transcripts = median(nFeature_RNA))

#median genes all dataset
seu_obj_clean@meta.data %>%
  summarize(median_transcripts = median(nFeature_RNA))

###

#figure 1B
DimPlot(seu_obj_clean, reduction = "umap", label = TRUE, 
        pt.size = 0.4, repel = T, label.box = F, 
        sizes.highlight = T, label.size = 4) + 
  #NoLegend()+
  #ggtitle("UMAP showing 18 distinct transcriptomic clusters")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Title
    axis.text = element_text(size = 12),      # Axis tick labels
    axis.title = element_text(size = 12)      # Axis titles (if shown)
  )


###


#cellProportions
seu_obj_clean$Cluster <- Idents(seu_obj_clean)

meta_df <- seu_obj_clean@meta.data

prop_df <- meta_df %>%
  group_by(sample, Cluster) %>%   # Use "Cluster" instead of "seurat_clusters"
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

palette_18 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#a6cee3", "#fb9a99",
  "#984ea3", "#ffff33", "#4daf4a", "#f781bf", "#999999", "#ffcc00"
)

ggplot(prop_df, aes(x = sample, y = proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Proportion of cells") +
  xlab("Sample") +
  ggtitle("Cluster composition across samples") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  scale_fill_manual(values = palette_18)

############# Figure 1 ends   #######################***************************
#








#






############# Figure 3 starts #######################***************************

#figure 3: marker genes

#list of two key marker genes for each cluster
gene_list <- c(
  "G21268", "G1208", "G11503", "G7179", "G21910", "G12025", "G3157", "G13652",
  "G6819", "G31525", "G2457", "G32588", "G384", "G22071", "G10190", "G25152",
  "G10388", "G4381", "G18439", "G25991", "G22661", "G26978", "G6983", "G5032",
  "G25637", "G3128", "G4180", "G7827", "G15965", "G15964", "G22673", "G27503",
  "G12733", "G5864"
)

# Figure 3A: DotPlot
DotPlot(seu_obj_clean, features = gene_list) +
  scale_color_gradient(low = "salmon", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  ggtitle("Expression of Marker Genes Across Clusters")


# figure3B: featureplot key marker genes, can't plot for all clusters, just pick key cell types

selected_genes <- c(
  "G21910", "G32588", "G384", 
  "G6983","G25637", "G7827", 
  "G15965", "G22673", "G12733"
)


FeaturePlot(seu_obj_clean, features = selected_genes, ncol = 3)

############# Figure 3 ends ##########################***************************
#
#


############# Figure 4 starts #######################***************************

# Figure 4: viral transcripts expression across clusters and infection stage


library(Seurat)
library(stringr)
library(forcats)

# get all viral gene names
viral <- rownames(seu_obj_clean) %>% 
  as.data.frame() %>% 
  filter(str_detect(., "^ORF" )) %>% 
  filter(!str_detect(., "\\." )) %>%
  rename(gene = ".")


# need clusters and infection stages in order (control, early, mid, late)..
# create a column in meta data to hold cluster info & condition
# we will use info from this column to compare groups
seu_obj_clean$sample_subtypes <- paste(Idents(seu_obj_clean), seu_obj_clean$condition_new, sep = "_")
Idents(seu_obj_clean) <- "sample_subtypes" # seurat object now has new Idents = old idents + condition_new
unique(seu_obj_clean$sample_subtypes)

# Define base order
base_order <- c("cluster0", "cluster1", "gill_ciliary",
                "hepatopancreas", "gill_nec", "gill_type1",
                "hyalinocytes", "haemocytes_type1",
                "cluster8", "cluster9", "mantle_vesicular",
                "immature_haemocytes", "macrophage_like",
                "adductor_muscle", "mantle", "digestive_gland",
                "gill_type2", "small_granules_cells")

# Define time order
time_order <- c("control", "early", "mid", "late")

# Combine to full desired order
desired_order <- as.vector(outer(base_order, time_order, paste, sep = "_"))

# Relevel the factor levels of the identities
Idents(seu_obj_clean) <- factor(Idents(seu_obj_clean), levels = desired_order)


# quick dotplot
DotPlot(seu_obj_clean, features = viral$gene, cols = c("#4575b4", "#d73027"), dot.scale = 8) +
  scale_color_gradientn(colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027")) +  # Smooth viridis-like gradient
  theme_minimal(base_size = 14) +  # Clean minimal theme with larger font
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    color = "Avg. Expression",
    size = "Pct. Expressed"
  )


#proper figure 4

gene_names <- viral$gene 

alt_labels <- ifelse(seq_along(gene_names) %% 2 == 0, gene_names, "") 
# Extract Idents (clusters) and clean labels 
cluster_labels <- levels(Idents(seu_obj_clean)) 

clean_cluster_labels <- str_remove(cluster_labels, "_(control|early|mid|late)$") 

DotPlot(seu_obj_clean, features = gene_names, cols = c("#4575b4", "#d73027"), dot.scale = 8) + 
  scale_color_gradientn(colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027")) + 
  scale_x_discrete(labels = alt_labels) + scale_y_discrete(labels = clean_cluster_labels) +  #  Removes suffix from y-axis ticks 
  theme_minimal(base_size = 14) + 
  theme( axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"), 
         axis.text.y = element_text(face = "bold"), axis.title = element_text(face = "bold"), 
         legend.position = "right", panel.grid.major = element_line(color = "gray90"), 
         panel.grid.minor = element_blank(), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) ) + 
  coord_flip(clip = "off") + 
  labs( x = NULL, y = NULL, color = "Avg. Expression", size = "Pct. Expressed")

############# Figure 4 ends #########################***************************

#
#



############# Figure 5 starts#########################***************************
# DGE analysis




############# Figure 5 ends #########################***************************


# for Figure 2
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
#......... gene symbols incorporated in gene IDs...
#
############# Figure 2 starts   #######################***************************

# Haemocytes sub-types in our data, based on Divonne et al

#haemocyte marker genes from Divonne et al

haem_markers <- read_csv("haemocytes_markers_cgigas.csv") %>% 
  select(gene, cluster) # cluster 1-7 in their data are haemocytes sub-types

# split by cluster, so that we can check each of their cluster markers in our data
haem_mark_list <- split(haem_markers$gene, haem_markers$cluster)


# plot top 30 genes for each sub-type Divonne et al, specify using haem_mark_list[[]]; 1-7
DotPlot(seu_obj_clean, features = haem_mark_list[[5]][1:30])+
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Marker genes (Divonne et al 2024)')+
  ggtitle("Expression of ")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12)) 

# summary Divonne et al markers in our clusters
#cluster1: not great; can't see any specific markers
#cluster2 : some good markers seen in our cluster6 (2182 nuclei); c("G32588", "G2457"); Hyalinocytes
#cluster3: c("G15881", "G12733", "G5864") in our cluster6; small granule cells
#cluster4: not great; mostly in our cluster 11 (named immature haemoctyes) but not very specific?
#cluster5: good markers in our cluster 7; c("G21444", "G22986", "G384"), need to find what this one is labelled in Divonne et al
#cluster6: not great
#cluster7: okayish, can see c("G1202", "G28618", "G13140", "G10862") in our cluster 10,12; vesicular cells



# plot top two marker genes from Divonne et al

# divonne et al cluster 1; macrophage-like ->  professional phagocytes
#" Cluster 1 demonstrated enrichment in viral processing and endocytosis"
DotPlot(seu_obj_clean, features = c("G8994", "G8841"))+
  coord_flip()+ RotatedAxis()

# divonne et al cluster 2; hyalinocytes
# Our findings suggest that hyalinocytes and/or the blast-like cells may be
# a cellular target of the OsHV-1 virus, the causal agent of POMS, which dampens the
# expression of certain AMPs
DotPlot(seu_obj_clean, features = c("G24846", "G17277"))+
  coord_flip()+ RotatedAxis()

# divonne et al cluster 3; small granule cells ->  professional phagocytes
DotPlot(seu_obj_clean, features = c("G22387", "G21091"))+
  coord_flip()+ RotatedAxis()

# divonne et al cluster 4
DotPlot(seu_obj_clean, features = c("G35227", "G31074"))+
  coord_flip()+ RotatedAxis()

# divonne et al cluster 5
DotPlot(seu_obj_clean, features = c("G14441", "G12472"))+
  coord_flip()+ RotatedAxis()

# divonne et al cluster 6
DotPlot(seu_obj_clean, features = c("G29512", "G3767"))+
  coord_flip()+ RotatedAxis()


# divonne et al cluster 7; vesicular cells
DotPlot(seu_obj_clean, features = c("G13140", "G32859"))+
  coord_flip()+ RotatedAxis()


#..........
# list markers here
divonne_markers <- c("G12639", "G12733", "G22387", "G5864", "G32289",
                     "G2459", "G29330", "G32588", "G2457",
                     "G4972", "G7023", "G1172", "G16773",
                     "G21444", "G22986", "G384",
                     "G28618", "G1202")

# get gene symbols
divonne_symbols <- gene_map_filtered %>%
  filter(gene_id %in% divonne_markers) %>%
  mutate(symbol_trimmed = sub("^\\S+\\s+", "", gene_symbol))  # removes gene ID for plot

# custom ggplot theme
theme_pub <- theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.4),
    plot.margin = margin(10, 10, 10, 10)
  )

# Create the dot plot
divonne <- DotPlot(seu_obj_clean, features = divonne_symbols$gene_symbol, cols = c("gray80", "firebrick3"), dot.scale = 6) +
  theme_pub +
  scale_x_discrete(labels = divonne_symbols$symbol_trimmed) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = NULL, color = "Avg. Expression", size = "Pct. Expressed")

ggsave("Figure 2A dotplot.pdf", plot = divonne, width = 10, height = 12, units = "in", dpi = 600)

#figure 2B
# Prepare gene list
features <- divonne_symbols$gene_symbol
names(features) <- divonne_symbols$symbol_trimmed

# Create stacked violin plot
stacked_vln <- VlnPlot(
  seu_obj_clean,
  features = features,
  stack = TRUE,
  flip = TRUE,
  pt.size = 0
) +
  NoLegend()

ggsave("Figure2B_stacked_violin.pdf", plot = stacked_vln, width = 10, height = 12, dpi = 600)

#...........
############# Figure 2 ends   #######################***************************
#
#
#
#
#
#
#
#
#



















# doHeatmap_plot <- DoHeatmap(seu_obj, features = gene_list, angle = 90, size = 3.2)+
#   scale_fill_gradientn(colors = c("black", "yellow")) + 
#   #NoLegend()+
#   theme(plot.margin = margin(t=60, 10, 10, 10)) +
#   theme(axis.text.y = element_blank()) #remove gene names on the left

# Clustered_DotPlot(seu_obj_clean, features = gene_list,
#                   row_label_size = 0,  cluster_ident = F, cluster_feature = F, x_lab_rotate = 90)
# 

#Stacked_VlnPlot(seu_obj_clean, features = gene_list,  x_lab_rotate = 90)


#######################################################
##
##
##
##













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

DimPlot(seu_obj_clean, reduction = "umap", label = F, 
        pt.size = 0.5, repel = T, label.box = F, 
        sizes.highlight = T, label.size = 6, group.by = "seurat_clusters") + 
  NoLegend()+
  ggtitle("Starsolo harmony umap : 18 dim x 6 res x 3k variable features")+
  theme(plot.title = element_text(hjust = 0.5))


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
      coalesce(Gene.ID, ""),
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
parent_dir <- "de_plots_full_ann_v2"
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

