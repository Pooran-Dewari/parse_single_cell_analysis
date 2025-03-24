

########## load harmony umap rds  ################################

#start here if coming here after few days

library(Seurat)
library(tidyverse)
library(Matrix)
#library (ggpubr)
#library(harmony)
library(scCustomize)

#####################################...........................

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

# load seurat object that we had created in previous step with pca & umap done!
seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")
seu_obj
head(seu_obj)


DimPlot(seu_obj, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("Starsolo harmony umap : 18 dim x 6 res x 3k variable features")+
  theme(plot.title = element_text(hjust = 0.5))


#################### find markers and annotate genes   ##########################

### find markers genes and annoatate using ORSON supl data from the French group

all_markers <- FindAllMarkers(object = seu_obj, min.pct = 0.20, log2fc.threshold = 0.25)
seu_obj

# convert all markers to df
df_all_markers <- as.data.frame(all_markers) %>% 
  filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "Gene.ID")

# annotate all_markers with ORSON
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv")

df_all_markers_ann <- left_join(df_all_markers, ORSON, by = "Gene.ID")

# now annotate with cg_science paper Piovani et al 2023
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  rename(Gene.ID = "X1", Gene.symbol.science = "X2" )

df_all_markers_ann_v2 <- left_join(df_all_markers_ann, cg_science, by = "Gene.ID") %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

#write
write_tsv(df_all_markers_ann_v2, "df_all_markers_ann_v2.tsv" )

#################### annotation done!!! #######################################



#.............. plot TF expression across clusters
df_all_markers_ann_v2_tf <- df_all_markers_ann_v2 %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("transcription factor|sox", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.25)

DotPlot(seu_obj, features = df_all_markers_ann_v2_tf$Gene.ID) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('TF')+
  ggtitle("Expression of TFs across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12)) 

VlnPlot(seu_obj, features = "G13515") 
AverageExpression(seu_obj, features = "G13515")

#........................................................



#..........................................................

# get lit of top 30 markers for each cluster
top30_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 30, named_vector = FALSE,
                                     data_frame = TRUE) 

top30_markers_ann <- left_join(top30_markers, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

write_tsv(top30_markers_ann, "star-strict-top30_markers_ann.tsv")

top5_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE,
                                     data_frame = TRUE) 
Clustered_DotPlot(seurat_object = seu_obj, features = top5_markers$gene, plot_km_elbow=F)


########## START: plot top markers per cluster ###############################

#method 1: manual

#get top 30 markers as a list
top30_mark_list <- split(top30_markers_ann$gene, top30_markers_ann$cluster)

#some checks
sum(is.na(top30_markers_ann$Sequence.Description)) #how many NA in ORSON
sum(is.na(top30_markers_ann$Gene.symbol.science)) #how many NA in science paper annotation in top 30

#dotplot; note some mito genes come up for cluster 0 in this case
#this is because our mito cut off was 5%, this cluster has high mito contribution
DotPlot(seu_obj, features = top30_mark_list[17]) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Top 30 marker genes')+
  ggtitle("Expression of top 30 marker genes, cluster 0")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))+
  theme(axis.line.y = element_line(arrow = arrow(type='closed', length = unit(10,'pt'))))





#method 2:batch fashion dotplots
# Create the directory if it doesn't exist
if (!dir.exists("clusters_dotplot_3k_harmony")) {
  dir.create("clusters_dotplot_3k_harmony")
}

# Loop through each cluster (from 1 to 20, with each index in `top30_mark_list`)
for (i in 1:18) {
  # Generate the dot plot for each cluster
  plot <- DotPlot(seu_obj, features = top30_mark_list[[i]]) +
    coord_flip() +
    RotatedAxis() +
    ylab('Cluster') +  
    xlab('Top 30 marker genes') +
    ggtitle(paste("Expression of top 30 marker genes, cluster", i - 1)) + # Adjust title for 0-based cluster indexing
    #theme_minimal() + # Sets a white background as base
    theme(
      plot.title = element_text(hjust = 0.5, size = 16), # Adjust title font size
      axis.text.y = element_text(size = 9), # Adjust y-axis text size
      axis.text.x = element_text(size = 10), # Adjust x-axis text size
      axis.title.x = element_text(size = 12), # Adjust x-axis label size
      axis.title.y = element_text(size = 12), # Adjust y-axis label size
      panel.grid = element_blank(), # Removes grid lines if needed
      panel.border = element_blank(), # Removes any borders around plot
      panel.background = element_rect(fill = "white", color = NA), # Ensures panel background is white
      plot.background = element_rect(fill = "white", color = NA), # Ensures entire plot background is white
      axis.line.y = element_line(arrow = arrow(type = 'closed', length = unit(10, 'pt')))
    )
  
  # Save each plot to the "clusters_dotplot" directory
  ggsave(filename = paste0("clusters_dotplot_3k_harmony/Cluster_", i - 1, "_DotPlot.png"), plot = plot, width = 8, height = 6)
  
  # Display the plot (if running interactively)
  print(plot)
}

########## END: plot top markers per cluster ###############################




# visualise cells by cluster to find infection stage-specific clusters, bubble plot
x <- prop.table(table(seu_obj@meta.data$sample, seu_obj@meta.data$seurat_clusters), margin = 2)
x <- as.data.frame(x) %>% 
  rename(sample = Var1, clusterID = Var2) %>% 
  mutate(percent = Freq * 100)

x0 <- x %>% 
  filter(clusterID == 1)

sum(x0$percent)

# Create a bubble plot without the fill legend but keeping the size legend
ggplot(x, aes(x = sample, y = factor(clusterID), size = percent, fill = factor(clusterID))) +
  geom_point(alpha = 0.7, shape = 21, color = "black") +
  scale_size(range = c(3, 15), name = "percent") +  # Keep the size legend
  guides(fill = "none") +  # Hide the fill legend
  theme_minimal(base_size = 18) +
  labs(title = "cluster size across samples", x = "Sample", y = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#.................................
##########################################################################################
##########################################################################################















### below, probably not very useful unless gene contributes to >30% reads!!

## find out highly expressed genes and decide if they need to be removed!!
sce <- as.SingleCellExperiment(seu_obj)
library(scater)
#identify top 20 expressed genes.. if >5-10%, they might be driving lot of variation
#remove them before clustering!
plotHighestExprs(sce, exprs_values = "counts", n=20)
#G32887
#G32877
#G14031
#G14032
#G6643

#G32881
#G32174
#G3716



# calculate and plot % of counts coming from a gene 
seu_obj[["rna_G32887"]] <- PercentageFeatureSet(seu_obj, "G32887")
head(seu_obj)
VlnPlot(seu_obj , "rna_G32887")


# batch process: calculate and plot % of counts coming from a gene
# Define the list of genes
gene_list <- c("G32887", "G32877", "G32174", "G14031", "G14032", "G6643", "G32881")

# Iterate through the gene list and apply the operations
for (gene in gene_list) {
  gene_rna <- paste0("rna_", gene)
  seu_obj[[gene_rna]] <- PercentageFeatureSet(seu_obj, features = gene)
}

VlnPlot(seu_obj, features = "rna_G32174")


#.....................................................................................
#......................................................................................
