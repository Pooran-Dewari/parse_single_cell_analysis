# working on this script: 8 Jan 2025


######################### STEP1: loading data and filtering ###################
#load seurat objects, merge, filter out rRNA genes, discard high mito cells


#load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(scCustomize)

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")


# load filtered seu obj
seu1 <- read_rds("seu_1.rds")
seu2 <- read_rds("seu_2.rds")
seu3 <- read_rds("seu_3.rds")
seu4 <- read_rds("seu_4.rds")
seu5 <- read_rds("seu_5.rds")
seu6 <- read_rds("seu_6.rds")
seu7 <- read_rds("seu_7.rds")
seu8 <- read_rds("seu_8.rds")


# merge 8 objects
merged_8 <- merge(x=seu1, 
                  y=list(seu2, seu3, seu4, seu5, seu6,seu7, seu8),
                  add.cell.ids = c("seu1", "seu2", "seu3", "seu4", "seu5", "seu6", "seu7", "seu8"))



# join the 8 layers into one for easy normalisation, this will bring all 8 objects into one layer of 'counts'
seu_obj <- JoinLayers(merged_8)


#dimensions in the object, number of genes and number of cells
dim(seu_obj)
#[1]   33534 1179648

# filtering by number of genes
seu_obj<-  subset(seu_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

# check what left after filtering
dim(seu_obj)
#[1] 33534 26459

#qc plot
QC_Plot_UMIvsGene(seurat_object = seu_obj, low_cutoff_gene = 200, high_cutoff_gene = 3000, low_cutoff_UMI = 300,
                  high_cutoff_UMI = 10000)



### add sample info
# Define a mapping from orig.ident to sample categories
orig_to_sample <- c(
  "C9" = "Uninfected", "C10" = "Uninfected", "C11" = "Homogenate", "C12" = "Homogenate",
  "D1" = "6-hpiA", "D2" = "6-hpiA", "D3" = "6-hpiD", "D4" = "6-hpiD",
  "D5" = "24-hpiA", "D6" = "24-hpiA", "D7" = "24-hpiJ", "D8" = "24-hpiJ",
  "D9" = "72-hpiJ", "D10" = "72-hpiJ", "D11" = "96-hpiE", "D12" = "96-hpiE"
)

# Assign sample labels to a new 'sample' column based on the 'orig.ident' column
seu_obj@meta.data$sample <- orig_to_sample[seu_obj@meta.data$orig.ident]


head(seu_obj)
tail(seu_obj)



table(seu_obj@meta.data$sample)

#get sample names
unique(seu_obj@meta.data$sample)

#  order the sample names
seu_obj$sample <- factor(seu_obj$sample, levels=c("Uninfected", "Homogenate",
                                                  "6-hpiA", "6-hpiD",
                                                  "24-hpiA", "24-hpiJ",
                                                  "72-hpiJ", "96-hpiE"))

table(seu_obj@meta.data$sample)




# test script to remove ribo RNA genes from the matrix
# first get ribosomal genes from gff using the bash below
# grep 'ribosomal RNA;' Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 | awk '{print $9}'| awk -F";|:" '{print $2}' | grep -v "MZ*" > ribo_rRNA_genes
# add 28S rRNA, can't find 18S rRNA
# grep '=28S' Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 | awk '{print $9}'| awk -F";|:" '{print $2}' | grep -v "MZ*" >> ribo_rRNA_genes 


ribo_genes <- read_table("ribo_rRNA_genes", col_names = F) 
ribo_genes_to_remove <- as.character(ribo_genes$X1)

#optional: add genes that are highly expressed to this list, adding G32889 which I discovered later in analysis that it ..
# ... accounts to >40% counts in each cell

ribo_genes_to_remove <- c(ribo_genes_to_remove, "G32889")


##### test on one of the seu objects/sub-lib first
# dim(seu1)
# seu1_filt <- seu1[!(rownames(seu1) %in% ribo_genes_to_remove), ]
# dim(seu1_filt)
# # check if removed correctly, should return TRUE if removed!
# dim(seu1)[1]-dim(seu1_filt)[1]==length(ribo_genes_to_remove)

dim(seu_obj)
#[1] 33534 26459
seu_obj_filt <- seu_obj[!(rownames(seu_obj) %in% ribo_genes_to_remove), ]

dim(seu_obj_filt)
# check if removed correctly, should return TRUE if removed!
dim(seu_obj)[1]-dim(seu_obj_filt)[1]==length(ribo_genes_to_remove)

#remove seu_obj to avoid mistakenly using this instead of the new filt object
rm(seu_obj)
###from now on, use seu_obj_filt object, this one has rRNA genes removed!

VlnPlot(seu_obj_filt, features = c("nFeature_RNA", "nCount_RNA"))

#get list of ALL mito genes

x <- as.data.frame(rownames(seu_obj_filt))

x <- x %>% 
  filter(str_detect(`rownames(seu_obj_filt)`, "MZ" ))

MZ_genes <- as.character(x$`rownames(seu_obj_filt)`)

mt.genes = c("ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", 
             "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB", MZ_genes) # got the MZ ones from a cluster that looks very weird!

mt.genes.present <- mt.genes[mt.genes %in% rownames(seu_obj_filt)]

#print(paste0(length(mt.genes.present),"/",length(mt.genes), " mito genes present in the count matrix"))
seu_obj_filt <- PercentageFeatureSet(seu_obj_filt,features = mt.genes.present, col.name = "percent.mt")

head(seu_obj_filt)
tail(seu_obj_filt)

# plot percent.mt per sample, can see some cells have >20% mito contribrution, jeez!!
VlnPlot(seu_obj_filt, features = c("percent.mt"),pt.size = 0, group.by = "sample") + 
  geom_hline(yintercept = 20) + NoLegend()

# FIlter out high mito cells
dim(seu_obj_filt)
#[1] 33466 26459
seu_obj_filt <- subset(seu_obj_filt, subset = percent.mt < 5)
dim(seu_obj_filt)
#[1] 33467 23740 ########### approx 3k nuclei removed..!!
VlnPlot(seu_obj_filt, features = c("percent.mt"),pt.size = 0, group.by = "sample") + 
  geom_hline(yintercept = 20) + NoLegend()


# code below is just a sanity check that ribo rRNA genes were removed successfully earlier

ribo.genes.present <- ribo_genes_to_remove[ribo_genes_to_remove %in% rownames(seu_obj_filt)]
seu_obj_filt <- PercentageFeatureSet(seu_obj_filt,features = ribo.genes.present, col.name = "percent.ribo")
head(seu_obj_filt) # should display percent.ribo as value zero across all rows
VlnPlot(seu_obj_filt, features = c("percent.ribo"),pt.size = 0, group.by = "sample") + geom_hline(yintercept = 20) + NoLegend()
#shouldn't see anything in the plot, all values are 0


#save seu object, this one has no ribo rRNA genes, plus ALL mito >5% removed!!
saveRDS(seu_obj_filt, file = "seu_obj_filt.rds")

rm(list=ls())
#...............................................................................
















































########### STEP2: Data normalisation, scaling, PCA #####################
#...............................................................................
# got this script from thomas
library(Seurat)
library(tidyverse)
library(Matrix)
#library (ggpubr)
#library(harmony)
library(scCustomize)

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

#standard seurat workflow until elbow plot

# load the filtered seu obj
seu_obj <- read_rds("seu_obj_filt.rds")

# log-transform (normalisation)
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)

# linear transformation (‘scaling’) 
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)
head(seu_obj)
# perform linear dimensional reduction (PCA)
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
DimPlot(seu_obj, reduction = "pca") + NoLegend()

seu_obj
saveRDS(seu_obj, file = "seu_obj_pca_3k.rds")

# can end session here and load the rds later to save on memory
rm(list=ls())


#################################################################################
library(Seurat)
library(tidyverse)
library(Matrix)
#library (ggpubr)
#library(harmony)
library(scCustomize)

#####################################...........................

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

# load seurat object that we had created in previous step with pca done!
seu_obj <- read_rds("seu_obj_pca_3k.rds")
seu_obj
# An object of class Seurat 
# 33466 features across 23740 samples within 1 assay 
# Active assay: RNA (33467 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 dimensional reduction calculated: pca

DimPlot(seu_obj, reduction = "pca") + NoLegend()



# decide on dims
#DimHeatmap(seu_obj, dims = 1:25, cells = 500, balanced = TRUE)

#use this to assess number of PC's to use in downstream analysis
ElbowPlot(seu_obj, ndims=50)

#looking at the elbo plot and heatmap, dims 15,20, 25,30 look reasonable to test!

#test dims=15

# going ahead with dims = 20 and res 0.6

seu_obj<- FindNeighbors(seu_obj, dims = 1:20)
seu_obj<- FindClusters(seu_obj, resolution = 0.6)


seu_obj <- RunUMAP(seu_obj, dims = 1:20)

#plot the UMAP
DimPlot(seu_obj, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("Starsolo umap: 20 dim x 6 res x 3k variable features")+
  theme(plot.title = element_text(hjust = 0.5))

seu_obj #should now show umap and pca in reduction layer



########### highlight clusters

# Define the directory name
dir_name <- "20dimx6resx3kFeatres"

# Source the script and pass `dir_name` as an argument
source("highlight.R")

# Call the `run_script()` function with the `dir_name`
run_script(dir_name)

saveRDS(seu_obj, file = "seu_obj_umap_20d_6r_3kRes.rds")

rm(list=ls())
######
#####
####
##
#















########## load umap rds  ################################

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
seu_obj <- read_rds("seu_obj_umap_20d_6r_3kRes.rds")
seu_obj

DimPlot(seu_obj, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("Starsolo umap original: 20 dim x 6 res x 3k variable features")+
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

VlnPlot(seu_obj, features = "G5010") 
AverageExpression(seu_obj, features = "G31799")

#........................................................



#..........................................................

# get lit of top 30 markers for each cluster
top30_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 30, named_vector = FALSE,
                                     data_frame = TRUE) 

top30_markers_ann <- left_join(top30_markers, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

write_tsv(top30_markers_ann, "star-strict-top30_markers_ann.tsv")

#compare with original seurat object
top_markers_original_seu <- Extract_Top_Markers(marker_dataframe = all_markers,
                                            num_genes = 5, named_vector = FALSE,
                                            make_unique = TRUE)
Clustered_DotPlot(seurat_object = seu_obj, features = top_markers_original_seu, plot_km_elbow=F)


########## START: plot top markers per cluster ###############################

#method 1: manual

#get top 30 markers as a list
top30_mark_list <- split(top30_markers_ann$gene, top30_markers_ann$cluster)

#some checks
sum(is.na(top30_markers_ann$Sequence.Description)) #how many NA in ORSON
sum(is.na(top30_markers_ann$Gene.symbol.science)) #how many NA in science paper annotation in top 30

#dotplot; note some mito genes come up for cluster 0 in this case
#this is because our mito cut off was 5%, this cluster has high mito contribution
DotPlot(seu_obj, features = top30_mark_list[2]) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Top 30 marker genes')+
  ggtitle("Expression of top 30 marker genes, cluster 0")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))+
  theme(axis.line.y = element_line(arrow = arrow(type='closed', length = unit(10,'pt'))))





#method 2:batch fashion dotplots
# Create the directory if it doesn't exist
if (!dir.exists("clusters_dotplot_3k")) {
  dir.create("clusters_dotplot_3k")
}

# Loop through each cluster (from 1 to 20, with each index in `top30_mark_list`)
for (i in 1:22) {
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
  ggsave(filename = paste0("clusters_dotplot_3k/Cluster_", i - 1, "_DotPlot.png"), plot = plot, width = 8, height = 6)
  
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

## find out highly expressed genes and decide if they need to be removed!!
sce <- as.SingleCellExperiment(seu_obj)
library(scater)
#identify top 20 expressed genes.. if >5-10%, they might be driving lot of variation
#remove them before clustering!
plotHighestExprs(sce, exprs_values = "counts", n=20)
#G32887
#G32877
#G32174
#G14031
#G14032
#G6643
#G32881


# calculate and plot % of counts coming from a gene 
seu_obj[["rna_G32887"]] <- PercentageFeatureSet(seu_obj, "GG32887")
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

VlnPlot(seu_obj, features = "rna_G32881")


#.....................................................................................
#......................................................................................












##################################################################
###################################################################
####### test remove odd cluster, i.e. 1 in our case and create test_seu
test_seu <- subset(x = seu_obj,invert=T, idents ="1")

# log-transform (normalisation)
test_seu <- NormalizeData(test_seu, normalization.method = "LogNormalize", scale.factor = 10000)

test_seu <- FindVariableFeatures(test_seu, selection.method = "vst", nfeatures = 3000)

# linear transformation (‘scaling’) 
all.genes <- rownames(test_seu)
test_seu <- ScaleData(test_seu, features = all.genes)
head(test_seu)
# perform linear dimensional reduction (PCA)
test_seu <- RunPCA(test_seu, features = VariableFeatures(object = test_seu))
DimPlot(test_seu, reduction = "pca") + NoLegend()

test_seu
saveRDS(test_seu, file = "test_seu_pca_3k.rds")

# can end session here and load the rds later to save on memory
rm(list=ls())


#################################################################################
library(Seurat)
library(tidyverse)
library(Matrix)
#library (ggpubr)
#library(harmony)
library(scCustomize)

#####################################...........................

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

# load seurat object that we had created in previous step with pca done!
test_seu <- read_rds("test_seu_pca_3k.rds")
test_seu
# An object of class Seurat 
# 33466 features across 21189 samples within 1 assay 
# Active assay: RNA (33466 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

DimPlot(test_seu, reduction = "pca") + NoLegend()



# decide on dims
#DimHeatmap(test_seu, dims = 1:25, cells = 500, balanced = TRUE)

#use this to assess number of PC's to use in downstream analysis
ElbowPlot(test_seu, ndims=50)

#looking at the elbo plot and heatmap, dims 15,20, 25,30 look reasonable to test!

#test dims=15

# going ahead with dims = 20 and res 0.6

test_seu <- FindNeighbors(test_seu, dims = 1:20)
test_seu <- FindClusters(test_seu, resolution = 0.6)


test_seu <- RunUMAP(test_seu, dims = 1:20)

DimPlot(test_seu, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("cluster 1 removed and reanalysed")+
  theme(plot.title = element_text(hjust = 0.5))

test_seu
saveRDS(test_seu, file = "test_seu_umap_20dim_6res_3kFeat.rds")

###########














#start here next day

library(Seurat)
library(tidyverse)
library(Matrix)
#library (ggpubr)
#library(harmony)
library(scCustomize)

#####################################...........................

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

# load seurat object that we had created in previous step with pca done!
test_seu <- read_rds("test_seu_umap_20dim_6res_3kFeat.rds")
test_seu
table(test_seu@meta.data$sample)


############# qc plots ###################################

VlnPlot(test_seu, features = "percent.mt", group.by = "sample", pt.size = 0, )+
  NoLegend()+
  ggtitle("Mitochondrial reads contribution across samples")+
  xlab("Sample") +
  ylab("% Mitochondrial reads")

VlnPlot(test_seu, features = "nFeature_RNA", group.by = "sample", pt.size = 0)+
  NoLegend()+
  ggtitle("Number of genes(features) detected across samples")+
  xlab("Sample") +
  ylab("Number of genes detected")

VlnPlot(test_seu, features = "nCount_RNA", group.by = "sample", pt.size = 0, y.max = 10000)+
  NoLegend()+
  ggtitle("Number of transcripts (UMIs) detected across samples")+
  xlab("Sample") +
  ylab("Number of transcripts detected")

#median UMI per sample
test_seu@meta.data %>%
  group_by(sample) %>%
  summarize(median_UMI = median(nCount_RNA))

#median UMI all dataset
test_seu@meta.data %>%
  summarize(median_UMI = median(nCount_RNA))


#median genes per sample
test_seu@meta.data %>%
  group_by(sample) %>%
  summarize(median_transcripts = median(nFeature_RNA))

#median genes all dataset
test_seu@meta.data %>%
  summarize(median_transcripts = median(nFeature_RNA))

########################################################
  
DimPlot(test_seu, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("cluster 1 removed and reanalysed: a lot cleaner umap")+
  theme(plot.title = element_text(hjust = 0.5))

### find markers genes and annoatate using ORSON supl data from the French group

all_markers_test_seu <- FindAllMarkers(object = test_seu, min.pct = 0.20, log2fc.threshold = 0.25)



##
# annotate all_markers with ORSON
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv")

# now annotate with cg_science paper
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  as.data.frame()

colnames(cg_science)[c(1, 2)] <- c("Gene.ID", "Gene.symbol.science")

all_markers_test_seu_ann <- left_join(all_markers_test_seu, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

# all_markers_test_seu_ann_13 <- all_markers_test_seu_ann %>% 
#   filter(cluster == 13)

all_markers_test_seu_ann_up_significant <- all_markers_test_seu_ann %>% 
  filter(avg_log2FC > 1.5 & pct.2 < 0.05 & pct.1 > 0.2)
    

  write_tsv(all_markers_test_seu_ann_up_significant, "all_markers_test_seu_ann_up_significant.tsv")
  

##


dim(test_seu)
    
top_markers_test_seu <- Extract_Top_Markers(marker_dataframe = all_markers_test_seu, 
                                            num_genes = 5, named_vector = FALSE,
                                  make_unique = TRUE)

library(scCustomize)
Clustered_DotPlot(seurat_object = test_seu, features = top_markers_test_seu, 
                  plot_km_elbow=F, k=1)

#compare with original seurat object
# top_markers_original_seu <- Extract_Top_Markers(marker_dataframe = all_markers, 
#                                             num_genes = 5, named_vector = FALSE,
#                                             make_unique = TRUE)
# Clustered_DotPlot(seurat_object = seu_obj, features = top_markers_original_seu, plot_km_elbow=F)


# get lit of top 30 markers for each cluster
top30_markers_test_seu <- Extract_Top_Markers(marker_dataframe = all_markers_test_seu, num_genes = 30, named_vector = FALSE,
                                     data_frame = TRUE) 

# annotate all_markers with ORSON
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv")

# now annotate with cg_science paper
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  as.data.frame()

colnames(cg_science)[c(1, 2)] <- c("Gene.ID", "Gene.symbol.science")





top30_markers_test_seu_ann <- left_join(top30_markers_test_seu, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

write_tsv(top30_markers_test_seu_ann, "star-strict-top30_markers_ann_test_seu.tsv")
#


# Default plot spacing (plot_spacing = 0.15 and spacing_unit = 'cm')
Stacked_VlnPlot(test_seu , features = c("G19458","G17783",
                                        "G26785", "G2482","G32688","G32686",
                                        "G29679"), x_lab_rotate = TRUE)


#
#
# get lit of top 50 markers for each cluster
top50_markers_test_seu <- Extract_Top_Markers(marker_dataframe = all_markers_test_seu, num_genes = 50, named_vector = FALSE,
                                              data_frame = TRUE) 

# annotate all_markers with ORSON
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv")

# now annotate with cg_science paper
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  as.data.frame()

colnames(cg_science)[c(1, 2)] <- c("Gene.ID", "Gene.symbol.science")



library(xlsx)
top50_markers_test_seu_ann <- left_join(top50_markers_test_seu, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description)

write.xlsx(top50_markers_test_seu_ann, "top50_markers_test_seu_ann.xlsx", sheetName="top 50 markers")

#################################




#17 Feb 2025

test_seu$celltype <- Idents(test_seu)

#combined/blended featureplot
FeaturePlot(test_seu, features = c( 'G5996', 'G5983' ),  blend = TRUE, label = F)

#if need to extract them
plots <- FeaturePlot(test_seu, features = c('G27279', 'G20558' ), combine = FALSE, blend = TRUE, label = TRUE)
plots[[3]] + NoLegend()  # Get just the co-expression plot, built-in legend is meaningless for this plot
plots[[4]] # Get just the key
CombinePlots(plots[3:4], legend = 'none') # Stitch the co-expression and key plots together

## we can calculate and plot coexpression of two genes
# Define genes of interest


gene1 <- "G5996"
gene2 <- "G5983"

# Fetch expression data, including cluster names
expr_data <- FetchData(test_seu, vars = c(gene1, gene2))  # No need to fetch "seurat_clusters"
expr_data$celltype <- test_seu$celltype  # Assign cell type labels

# Identify co-expressing cells
expr_data$coexpressed <- (expr_data[[gene1]] > 0) & (expr_data[[gene2]] > 0)

# Compute co-expression percentage per cluster (now using names)
coexp_pct <- expr_data %>%
  group_by(celltype) %>%
  summarize(pct_coexpressed = sum(coexpressed) / n() * 100)



colors <- c(
  "red", "blue", "green", "purple", "orange", "pink", "cyan", "brown", "magenta", "yellow",
  "gray", "darkgreen", "navy", "gold", "coral", "lightblue", "darkred", "violet", "seagreen2", "black", "salmon"
)

ggplot(coexp_pct, aes(x = celltype, y = pct_coexpressed, fill = factor(celltype))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = paste("Percentage of cells co-expressing", gene1, "and", gene2, "in each cluster"),
       x = "Cell Type",
       y = "Co-expression Percentage") +
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate labels for readability
  scale_fill_manual(values = colors) +  # Assign custom colors
  NoLegend()

#################







VlnPlot(test_seu, "G15660")

cluster_7_top50 <- top50_markers_test_seu_ann %>% 
  filter(cluster == 7)
library(xlsx)
write.xlsx(cluster_7_top50, "cluster_7_top50.xlsx", sheetName="top 50 markers cluster7")

#plot some cluster 7 genes
VlnPlot(test_seu, features = c("G1208", "G1206", "G4113", "G31799", "G8546", "G17926", "G8077"))
FeaturePlot(test_seu, features = "G8546")
#

#tf
#.............. plot all TFs expression across clusters
all_markers_test_seu_ann_tf <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("transcription factor|sox", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.25)

DotPlot(test_seu, features = unique(all_markers_test_seu_ann_tf$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('TF')+
  ggtitle("Expression of TFs across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))


#plot only those TFs that are in top 30 marker genes in clusters
top30_markers_test_seu_ann_tf <- top30_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("transcription factor|sox", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.25)

DotPlot(seu_obj, features = unique(top30_markers_test_seu_ann_tf$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('TF')+
  ggtitle("Expression of TFs across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#
#.............. plot HSP expression across clusters
#hsp
df_all_markers_test_seu_ann_hsp <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("heat shock", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(df_all_markers_test_seu_ann_hsp$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('HSP')+
  ggtitle("Expression of HSPs across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12)) 

#

#.............. plot chitin synthase expression across clusters
#chitin
df_all_markers_test_seu_ann_chitin <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("chitin", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(df_all_markers_test_seu_ann_chitin$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('chitin')+
  ggtitle("Expression of chitin* across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12)) 

#
#.............. plot vwf domain-containing expression across clusters
#vwf domain
df_all_markers_test_seu_ann_vwf <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("vwf", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(df_all_markers_test_seu_ann_vwf$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('vwf domain')+
  ggtitle("Expression of vwf domain across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#
#.............. plot egf-like synthase expression across clusters
#egf
#df_all_markers_test_seu_ann_egf <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("egf", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(df_all_markers_test_seu_ann_egf$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('egf')+
  ggtitle("Expression of egf* across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12)) 

#
#.............. plot beta lactamse domain-containing expression across clusters
#lactamase
df_all_markers_test_seu_ann_lact <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("lactamase", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(df_all_markers_test_seu_ann_lact$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('lactamase domain')+
  ggtitle("Expression of beta lactamase domain across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))
#
#fibrocystin

test_seu_fibrocystin <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("fibrocystin", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_fibrocystin$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('fibrocystin-L')+
  ggtitle("Expression of fibrocystin-L across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))
#
#c1q

test_seu_Lactase_phlorizin <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("Lactase-phlorizin hydrolase", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_Lactase_phlorizin$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Lactase-phlorizin hydrolase')+
  ggtitle("Expression of Lactase-phlorizin hydrolase across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))
#
#hexosaminidase
test_seu_Lactase_hexosaminidase <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("hexosaminidase", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_Lactase_hexosaminidase$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('beta hexosaminidase')+
  ggtitle("Expression of beta hexosaminidase across clusters")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#gills
#cathepsin L
test_seu_cathepsin <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("cathepsin L", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_cathepsin$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('cathepsin')+
  ggtitle("Expression of cathepsin L across clusters: Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#c-type lectin
test_seu_c_type_lectin <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("c-type lectin", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_lectin$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('c-type lectin: Gills')+
  ggtitle("Expression of C-type lectin clusters: Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#solute carrier
test_seu_c_type_SLC <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("solute carrier", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_SLC$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('solute carrier genes')+
  ggtitle("Expression of solute carrier across all clusters: Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))


#
#slc6
test_seu_c_type_slc6 <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("slc6", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_slc6$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('slc6')+
  ggtitle("Expression of slc6 across all clusters: Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#mucin-5
test_seu_c_type_mucin5 <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("mucin-5", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_mucin5$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('mucin-5')+
  ggtitle("Expression of mucin-5 across all clusters: Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))


##mucin
test_seu_c_type_mucin <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("mucin", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_mucin$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('mucin')+
  ggtitle("Expression of mucins across all clusters: Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#slc6
test_seu_c_type_slc6 <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("slc6", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_slc6$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('slc6')+
  ggtitle("Expression of serotonin transporter, SERT (slc6a4a) across all clusters: Neuroepithelial cells Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))


#rgs, regulator of G-protein signaling
test_seu_c_type_rgs <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("rgs", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_rgs$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('rgs')+
  ggtitle("Expression of rgs across all clusters: Neuroepithelial cells Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#kcn: potassium channel
test_seu_c_type_kcn <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("kcn", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_kcn$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('kcn')+
  ggtitle("Expression of kcn (potassium channel) across all clusters: Neuroepithelial cells Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#synaptotagmin
test_seu_c_type_synapto <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("synaptotagmin", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_synapto$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('synaptotagmin')+
  ggtitle("Expression of synaptotagmin across all clusters: Neuroepithelial cells Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

#meis2
test_seu_c_type_meis2 <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("meis2", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_meis2$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('meis2')+
  ggtitle("Expression of meis2 TF across all clusters: Neuroepithelial cells Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))


#trpv6
test_seu_c_type_trpv6 <- all_markers_test_seu_ann %>% 
  filter(apply(., 1, function(row) any(str_detect(row, regex("trpv6", ignore_case = TRUE))))) %>% 
  filter(avg_log2FC > 0.5)

DotPlot(test_seu, features = unique(test_seu_c_type_trpv6$gene)) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('trpv6')+
  ggtitle("Expression of trpv6across all clusters: Neuroepithelial cells Gills")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12))

########## START: plot top markers per cluster ###############################


#get top 30 markers as a list
top30_mark_list_test_seu <- split(top30_markers_test_seu_ann$gene, top30_markers_test_seu_ann$cluster)

#some checks
sum(is.na(top30_markers_test_seu_ann$Sequence.Description)) #how many NA in ORSON
sum(is.na(top30_markers_test_seu_ann$Gene.symbol.science)) #how many NA in science paper annotation in top 30

#################################################################

#batch fashion dotplots
# Create the directory if it doesn't exist
if (!dir.exists("clusters_dotplot_3k_subset1")) {
  dir.create("clusters_dotplot_3k_subset1")
}

# Loop through each cluster (from 1 to 20, with each index in `top30_mark_list`)
for (i in 1:21) {
  # Generate the dot plot for each cluster
  plot <- DotPlot(test_seu, features = top30_mark_list_test_seu[[i]]) +
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
  ggsave(filename = paste0("clusters_dotplot_3k_subset1/Cluster_", i - 1, "_DotPlot.png"), plot = plot, width = 8, height = 6)
  
  # Display the plot (if running interactively)
  print(plot)
}

########## END: plot top markers per cluster ###############################


########### highlight clusters

# Define the directory name
dir_name <- "20dimx6resx3kFeatres_test_seu"

# Source the script and pass `dir_name` as an argument
source("highlight_test_seu.R")

# Call the `run_script()` function with the `dir_name`
run_script(dir_name)
#################################################################

VlnPlot(test_seu, features = "G3629")

############# done up to this mark on 10 Jan 2025
#########
#######
#####
###
#





## for shiny go
# get background genes
#all_genes <- as.data.frame(rownames(seu_obj))
#write_tsv(all_genes, "background_genes.txt")


## for rank-based go analysis
ORSON <- read_csv("../ORSON_french_group_2024_biorxiv_suppl2.csv")

custom_ann <- ORSON %>% 
  
  select(Gene.ID, Annotation.GO.ID) %>% 
  
  filter(Annotation.GO.ID != "NA")


write.table(custom_ann, "custom_annot.tab", row.names = FALSE, quote = FALSE, sep = "\t")
#write_tsv(custom_ann, "table_go_ann.tab", col_names = F)

gene_df <- read_tsv("../star-top20_markers_ann.tsv") 

gene_list <- gene_df %>% 
  #ilter(cluster == 5) %>% 
  select(gene, avg_log2FC) 

write.csv(gene_list, "genes_of_int.csv", quote = FALSE, row.names = FALSE)

########### highlight clusters

# Define the directory name
dir_name <- "20dimx6resx3000feat"

# Source the script and pass `dir_name` as an argument
source("highlight.R")

# Call the `run_script()` function with the `dir_name`
run_script(dir_name)

###################...................................................



seu_obj
# An object of class Seurat 
# 33467 features across 40191 samples within 1 assay 
# Active assay: RNA (33467 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap





#haemocyte marker genes from Divonne et al

haem_markers <- read_csv("haemocytes_markers_cgigas.csv") %>% 
  select(gene, cluster)

haem_mark_list <- split(haem_markers$gene, haem_markers$cluster)


# plot top 30 genes
DotPlot(test_seu, features = haem_mark_list[[7]][1:30])+
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Marker genes (Divonne et al 2024)')+
  ggtitle("Expression of ")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=12)) 
############### summary Divonne et al markers in our clusters
#cluster1: not great; can't see any specific markers; macrophage-like cells
#cluster2 : some good markers seen in our cluster6 (2182 nuclei); c("G32588", "G2457"); Hyalinocytes
#cluster3: c("G15881", "G12733", "G5864") in our cluster6; small granule cells
#cluster4: not great; mostly in cluster9 but not very specific?
#cluster5: good markers in our cluster9; c("G21444", "G22986", "G384")
#cluster6: not great
#cluster7: okayish, can see c("G1202", "G28618", "G13140") in our cluster 8,15; vesicular cells

gene_of_int =  c("G5864", "G4920", "G6647", "G10977", "G29519", "G16036", "G13885", "G16034", "G15342", "G31986", "G24730", "G8152", "G25378")
gene_of_int =  "G1202"
VlnPlot(seu_obj, features = gene_of_int)+
  xlab('Cluster') +  ylab('Expression level')

DotPlot(seu_obj, features = gene_of_int)
VlnPlot(test_seu, features = "G4972")


as.data.frame(table(seu_obj@meta.data$seurat_clusters))
# Var1 Freq
# 1     0 5451
# 2     1 3549
# 3     2 3414
# 4     3 3137
# 5     4 3106
# 6     5 2947
# 7     6 2182
# 8     7 2106
# 9     8 2074
# 10    9 1582
# 11   10 1580
# 12   11 1484
# 13   12 1402
# 14   13 1305
# 15   14 1027
# 16   15 1017
# 17   16  940
# 18   17  784
# 19   18  623
# 20   19  481

#VlnPlot(seu_obj, features = "G21444", split.by = "Samples")

#
#


############## 14 Jan 2025
#assign and label cell type based on markers
DimPlot(test_seu, reduction = "umap", label=T, label.box = F)

# new.cluster.ids <- c(new.cluster.ids <- c("Gill_ciliary (0)", "Hepatopancreas (1)", 
#                                           "Cluster 2 (2)",
#                                           "Gill_NEC (3)", "Gill_type 1 (4)", "Hyalinocyte (5)", 
#                                           "Haemocyte 1 (6)","Infection_response (7)", 
#                                           "Cluster 8 (8)",
#                                           "Immature_haemocyte (9)", 
#                                           "Cluster 10 (10)", 
#                                           "Adductor (11)", 
#                                           "Cluster 12 (12)", 
#                                           "Macrophage_like (13)", "Mantle (14)",
#                                           "Gill_type 2 (15)", "Gill_type 3 (16)", 
#                                           "Digestive_gland (17)", "Small_granule_cell (18)",
#                                           "Cell_cycle_related (19)", "Cluster 20 (20)"))

new.cluster.ids <- c(new.cluster.ids <- c("Gill_ciliary (0)", "Hepatopancreas (1)", 
                                          "Cluster 2",
                                          "Gill_NEC (3)", "Gill_type 1 (4)", "Hyalinocyte (5)", 
                                          "Haemocyte 1 (6)","Infection_response (7)", 
                                          "Cluster 8",
                                          "Immature_haemocyte (9)", 
                                          "Cluster 10", 
                                          "Adductor (11)", 
                                          "Cluster 12", 
                                          "Macrophage_like (13)", "Mantle (14)",
                                          "Gill_type 2 (15)", "Gill_type 3 (16)", 
                                          "Digestive_gland (17)", "Small_granule_cell (18)",
                                          "Cell_cycle (19)", "Cluster 20"))

names(new.cluster.ids) <- levels(test_seu)

test_seu <- RenameIdents(test_seu, new.cluster.ids)
table(test_seu$seurat_clusters)
DimPlot(test_seu, reduction = "umap", label = T, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = F, label.size = 6) + 
  #NoLegend()+
  ggtitle("Single-nucleus atlas of Paciifc oyster")+
  theme(plot.title = element_text(hjust = 0.5))

my_colors <- c("Gill_ciliary_cells(0)"= "deeppink1", "Hepatopancreas(1)"="aquamarine4", 
               "Mantle?(2)"= "royalblue1", "Gill_NEC(3)"= "deeppink1", 
               "Gill_type1(4)"= "deeppink1", "Hyalinocytes(5)"="darkorange", 
               "Macrophage_like_BBL(6)"="darkorange","Infection_response_7"="darkorange", 
               "Mantle_or_Vesicular_cells?(8)"="darkorange","Immature_haemocytes(9)"="darkorange", 
               "Mantle_sensory?(10)"= "royalblue1", "Adductor(1)"="pink3", 
               "Gill_NEC?(12)"= "deeppink1", "Macrophage_like(13)"="darkorange", 
               "Mantle(14)"= "royalblue1", "Vesicular_cells(15)"="darkorange", 
               "Gill_type2?(16)"= "deeppink1", "Haemocytes?(17)"="darkorange",
               "Small_granule_cells(18)"="darkorange", "Cell_cycle_related(19)"="cyan1",
               "Hepatopancreas_gills?(20)"="brown")


DimPlot(test_seu, reduction = "umap", label = T, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = F, cols = my_colors) + 
  #NoLegend()+
  ggtitle("Single-nucleus atlas of Paciifc oyster")+
  theme(plot.title = element_text(hjust = 0.5))



######## 20 Feb 2025
#key markers

# clear cut cell type assigned for these
VlnPlot(test_seu, features = c("G2457", "G32588")) #cluster5 Hyalinocytes
VlnPlot(test_seu, features = c("G12733", "G5864")) #cluster18 SGC
VlnPlot(test_seu, features = c("G384", "G22986")) #cluster 6, haemocyte 1
VlnPlot(test_seu, features = c("G13140", "G25309")) #cluster 15 gills 2
VlnPlot(test_seu, features = c("G32748", "G20782")) #cluster 9 immature haemocytes
VlnPlot(test_seu, features = c("G6983", "G2310")) #cluster 13 macrophage-like
VlnPlot(test_seu, features = c("G8077", "G8546")) #cluster 7 immune reponse
VlnPlot(test_seu, features = c("G7179", "G11503" )) #cluster 0 gill ciliary
VlnPlot(test_seu, features = c("G3157", "G13652" )) #cluster 3 gill NEC
VlnPlot(test_seu, features = c("G27156", "G6819" )) #cluster 4 gill 1
VlnPlot(test_seu, features = c("G22673", "G27503" )) #cluster 16 gill 3
VlnPlot(test_seu, features = c("G7827", "G4180" )) #cluster 14 mantle
VlnPlot(test_seu, features = c("G25637", "G3128" )) #cluster 11 adductor muscle
VlnPlot(test_seu, features = c("G21910", "G12025" )) #cluster 1 hepatopancreas
VlnPlot(test_seu, features = c("G15965", "G15964" )) #cluster 17 digestive gland

#these clusters couldn't be assigned to any specific cell type/tissue
VlnPlot(test_seu, features = c("G1208", "G21268" ))  #cluster 2; not very specific markers though
VlnPlot(test_seu, features = c("G18439", "G25991" ))  #cluster 8; look highly specific
VlnPlot(test_seu, features = c("G10190", "G25152"))  #cluster 10
VlnPlot(test_seu, features = c("G10388", "G4381" ))  #cluster 12
# cluster 19 is cyclin genes, so no point getting markers for this cluster
VlnPlot(test_seu, features = c("G27937","G27023"))  #cluster 20

#1st gene IDs of marker genes from above
first_marker <- c("G2457", "G12733", "G384", "G13140", "G32748", "G6983", "G8077", "G7179", 
  "G3157", "G27156", "G22673", "G7827", "G25637", "G21910", "G15965", 
  "G1208", "G18439", "G10190", "G10388", "G27937")

#2n gene IDs of marker genes from above
second_marker <- c("G32588", "G5864", "G22986", "G25309", "G20782", "G2310", "G8546", "G11503", 
  "G13652", "G6819", "G27503", "G4180", "G3128", "G12025", "G15964", 
  "G21268", "G25991", "G25152", "G4381", "G27023")

#both markers
both_markers <- c(first_marker, second_marker)
#not make the stacked violin plot
#Stacked_VlnPlot(test_seu , features = first_marker)

# Generate the stacked violin plot and store it in a variable
stacked_plot <- Stacked_VlnPlot(test_seu, features = first_marker,  x_lab_rotate = TRUE)+
  coord_flip()

# Save the plot as a high-quality PNG
ggsave(filename = "stacked_vlnplot.png", plot = stacked_plot, width = 10, height = 8, dpi = 300)

Clustered_DotPlot(test_seu, features = both_markers)
# Or save as PDF for vector quality
#ggsave(filename = "stacked_vlnplot.pdf", plot = stacked_plot, width = 10, height = 8)





#################################################################################################################
#doheatmap figure1 for PUBLICATION
doHeatmap_plot <- DoHeatmap(test_seu, features = top_markers_test_seu, angle = 90, size = 3.2)+
  scale_fill_gradientn(colors = c("black", "yellow")) + 
  #NoLegend()+
  theme(plot.margin = margin(t=60, 10, 10, 10)) +
  theme(axis.text.y = element_blank()) #remove gene names on the left

#save plot
ggsave("figures/doHeatmap_plot.pdf", plot = doHeatmap_plot, width = 12, height = 8, dpi = 300)
ggsave("figures/doHeatmap_plot.png", plot = doHeatmap_plot, width = 12, height = 8, dpi = 300)

#......................

feature_plot <- FeaturePlot(test_seu, features = first_marker)
ggsave("figures/feature_plot.pdf", plot = feature_plot, width = 12, height = 8, dpi = 300)
ggsave("figures/feature_plot.png", plot = feature_plot, width = 12, height = 8, dpi = 300)

#...................


#...................
#............................
#clustered dotplot for the top 5 genes per cluster
clustered_dotplot_list <- Clustered_DotPlot(test_seu, features = top_markers_test_seu, show_ident_colors = F,
                  row_label_size = 0,  cluster_ident = F, cluster_feature = F, x_lab_rotate = 90)

clustered_dotplot <- clustered_dotplot_list[[2]] #second item in the list is the plot we want

#note: ggsave can't save this class of object!!
#png("figures/clustered_dotplot.png", width = 20, height = 30, units="cm", res = 300) #if png needed
pdf("figures/clustered_dotplot.pdf")
library(ComplexHeatmap)

draw(clustered_dotplot)
dev.off()

#................................................

#stacked violin plots for two markers

stacked_plot <- Stacked_VlnPlot(test_seu, features = first_marker,  x_lab_rotate = 90)

# Save the plot as a high-quality PNG
ggsave(filename = "figures/stacked_vlnplot.png", plot = stacked_plot, width = 10, height = 8, dpi = 300)
ggsave(filename = "figures/stacked_vlnplot.pdf", plot = stacked_plot, width = 10, height = 8, dpi = 300)

#try two genes each for all haemocyte populations
stacked_plot_haem <- Stacked_VlnPlot(test_seu, features = c("G2457", "G32588",
                                                        "G12733", "G5864",
                                                        "G384", "G22986", "G32748", "G20782",
                                                        "G6983", "G2310","G8077", "G8546"), x_lab_rotate = 90)

ggsave("figures/stacked_plot_haem.pdf", plot = stacked_plot_haem, width = 16, height = 8, dpi = 300)

#gill
stacked_plot_gill <- Stacked_VlnPlot(test_seu, features = c("G13140", "G25309","G7179", 
                                                            "G11503" ,"G3157", "G13652",
                                                            "G27156", "G6819" ,"G22673", 
                                                            "G27503"), x_lab_rotate = 90)

ggsave("figures/stacked_plot_gill.pdf", plot = stacked_plot_gill, width = 16, height = 8, dpi = 300)

#other tissues
stacked_plot_tissues <- Stacked_VlnPlot(test_seu, features = c("G7827", "G4180",
                                                               "G25637", "G3128",
                                                               "G21910", "G12025",
                                                               "G15965", "G15964"), x_lab_rotate = 90)

ggsave("figures/stacked_plot_tissues.pdf", plot = stacked_plot_tissues, width = 16, height = 8, dpi = 300)


#unassigned clusters
stacked_plot_unassign <- Stacked_VlnPlot(test_seu, features = c("G1208", "G21268", "G18439",
                                                                "G25991", "G10190", "G25152",
                                                                "G10388", "G4381", "G27937",
                                                                "G27023"), x_lab_rotate = 90)

ggsave("figures/stacked_plot_unassign.pdf", plot = stacked_plot_unassign, width = 16, height = 8, dpi = 300)


#..............................

#umap with label box
umap_plot <- DimPlot(test_seu, reduction = "umap", label = T, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = F, label.size = 6) + 
  #NoLegend()+
  ggtitle("Single-nucleus atlas of Paciifc oyster")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/umap.pdf", plot = umap_plot, width = 16, height = 12, dpi = 300)
ggsave("figures/umap.png", plot = umap_plot, width = 16, height = 12, dpi = 300)

#umap without label box
umap_plot1 <- DimPlot(test_seu, reduction = "umap", label = F, 
                     pt.size = 0.5, repel = T, label.box = T, 
                     sizes.highlight = F, label.size = 6) + 
  #NoLegend()+
  ggtitle("Single-nucleus atlas of Paciifc oyster")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/umap1.pdf", plot = umap_plot1, width = 16, height = 12, dpi = 300)
ggsave("figures/umap1.png", plot = umap_plot1, width = 16, height = 12, dpi = 300)
########################################################################


#finally make stacked violin plot for figure
Stacked_VlnPlot(test_seu , features = c("G2457","G32588", #hyalinocytes
                                        "G12733", "G5864", "G572", "G30939", #SGC
                                        "G384", "G22986", #haemocyte 1
                                        "G13140", "G25309", #most likely gills
                                        "G22661", "G32748", "G20782", #immature haemocytes
                                        "G6983", "G2310", #macrophage-like
                                        "G8077", "G8546"),  #immune response
                x_lab_rotate = TRUE)


#17 Jan 2025
# quick check viral transcripts

# get all viral gene names
viral <- rownames(test_seu) %>% 
  as.data.frame() %>% 
  filter(str_detect(., "ORF" )) %>% 
  filter(!str_detect(., "\\." )) %>%
  rename(gene = ".")

# quick dotplot
DotPlot(test_seu, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) +
  theme(axis.text.x= element_text(angle=90, hjust=1))+
  coord_flip()

table(test_seu$sample)

#uninf control
Uninf_ctrl <- subset(test_seu, subset = sample == "Uninfected")
DotPlot(Uninf_ctrl, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


Homog_ctrl <- subset(test_seu, subset = sample == "Homogenate")
DotPlot(Homog_ctrl, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

hpi_6A <- subset(test_seu, subset = sample == "6-hpiA")
DotPlot(hpi_6A, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


hpi_6D <- subset(test_seu, subset = sample == "6-hpiD")
DotPlot(hpi_6D, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

hpi_24A <- subset(test_seu, subset = sample == "24-hpiA")
DotPlot(hpi_24A, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Gene')+
  ggtitle("viral transcripts: 24hpi-A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=10)) 

hpi_24J <- subset(test_seu, subset = sample == "24-hpiJ")
DotPlot(hpi_24J, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + 
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Gene')+
  ggtitle("viral transcripts: 24hpi-J")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=10)) 


hpi_72J <- subset(test_seu, subset = sample == "72-hpiJ")
DotPlot(hpi_72J, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Viral Gene')+
  ggtitle("Dotplot: Expression of viral transcripts in sample 72hpi-J")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=12)) 

#subset cluster 7 and check viral transcript expression
cluster7_test_seu <- subset(test_seu, idents = "Infection_response_7")

DotPlot(cluster7_test_seu, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + 
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Viral Gene')+
  ggtitle("Dotplot: Expression of viral transcripts in cluster 7")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=12))

table(cluster7_test_seu$sample)

VlnPlot(test_seu, features = "ORF124")



# no of cells per cluster
table(hpi_72J@meta.data$seurat_clusters)

hpi_96E <- subset(test_seu, subset = sample == "96-hpiE")
DotPlot(hpi_96E, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + 
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Viral Gene')+
  ggtitle("Dotplot: Expression of viral transcripts in sample 96hpi-E")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=12)) 

table(hpi_96E@meta.data$seurat_clusters)

### gimme motif analysis
#19 Jan 2025

## gimme analysis
# get all genes expressed log2fc >0.5 and adj pval < 0.05

#all_markers_test_seu generated in the section: "### find markers genes and annoatate using ORSON supl data from the French group"
all_markers_test_seu_up <- all_markers_test_seu  %>% 
  filter(avg_log2FC > 0.5) %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(!str_detect(gene, "MZ" ))

write_tsv(all_markers_test_seu_up, "all_markers_up_test_seu_for_gimme.tsv", col_names = T)

all_markers_test_seu_down <- all_markers_test_seu  %>% 
  filter(avg_log2FC < -0.5) %>% 
  filter(p_val_adj < 0.01)%>% 
  filter(!str_detect(gene, "MZ" ))

write_tsv(all_markers_test_seu_down, "all_markers_down_test_seu_for_gimme.tsv", col_names = T)

#

read_table("/external2/Pooran/seurat_thomas_scripts/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3")

gff <- as.data.frame(rtracklayer::import("/external2/Pooran/seurat_thomas_scripts/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3"))

#TSS is the first bit in 5'UTR, but we don't have gene ID for 5'UTR entry
#I can see 5'UTR and gene all start with same cooridnate, so will filter gene string
#and get TSS info from that
gff_v2 <- gff %>% 
  select(1:3, 5, 7, 10) %>% 
  filter(type == "gene") %>%
  separate(col = "ID", sep = ":", into = c("feature", "gene"), remove = FALSE) %>% 
  select(1:4, 8) %>% 
  mutate(start = start-200) %>% 
  mutate(end = start+202) %>%  #promoter is roughly 200 bp upstream of TSS
  rename(promoter = start, tss_end = end)


all_markers_test_seu_up_bed <- inner_join(gff_v2, all_markers_test_seu_up, by = "gene", relationship = "many-to-many") %>% 
  select(1:3,gene, avg_log2FC, strand, cluster) 


#save as bed file by cluster into a dir

library(dplyr)
library(purrr)

# Specify the output directory
output_dir <- "bed_clusters_up"

# Ensure the directory exists (create it if not)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Split the data by the cluster column
clusters <- split(all_markers_test_seu_up_bed, all_markers_test_seu_up_bed$cluster)

# Save each cluster as a separate file in the directory
walk(names(clusters), function(cluster_name) {
  # Construct the file path
  file_path <- file.path(output_dir, paste0("cluster_", gsub("[[:space:][:punct:]]+", "_", cluster_name), ".bed"))
  
  # Write the cluster data to a file
  write.table(
    clusters[[cluster_name]], 
    file = file_path,
    row.names = FALSE, 
    sep = "\t", 
    quote = FALSE,
    col.names = F
  )
})

# #if the batch process saving into dir doesn't work, use this
# library(dplyr)
# library(purrr)
# 
# # Split the data by the cluster column
# clusters <- split(df, df$cluster)
# 
# # Save each cluster as a separate file
# walk(names(clusters), function(cluster_name) {
#   write.table(
#     clusters[[cluster_name]], 
#     file = paste0("cluster_", gsub("\\s|\\(|\\)", "_", cluster_name), ".txt"),
#     row.names = FALSE, 
#     sep = "\t", 
#     quote = FALSE
#   )
# })

##

#write_tsv(all_markers_test_seu_up_bed, "top_markers_20_bed.bed", col_names = FALSE)

## getFasta for top markes 20 bed, bash programming
#bedtools getfasta -name -s -fi Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa -bed top_markers_20_bed.bed -fo top_markers_20.fasta


##

#3feb2025
#plot starsolo summary stats
library(readxl)

starsolo_sum <- read_xlsx("starsolo_summary.xlsx", col_names = F)
colnames(starsolo_sum) <- c("sublibrary", "description", "reads")
starsolo_sum$sublibrary[1] <- "a1"
starsolo_sum <- starsolo_sum %>% 
  drop_na()

# Assuming starsolo_sum is your DataFrame

# Create a named vector storing 'Number of Reads' for each sublibrary
num_reads_dict <- setNames(starsolo_sum$reads[starsolo_sum$description == "Number of Reads"], 
                           starsolo_sum$sublibrary[starsolo_sum$description == "Number of Reads"])

# Create a new column that divides 'reads' by 'Number of Reads' for the corresponding sublibrary
starsolo_sum$Percentage <- ifelse(
  starsolo_sum$description != "Number of Reads",
  starsolo_sum$reads / num_reads_dict[starsolo_sum$sublibrary] * 100,
  NA
)

# Display the updated DataFrame
print(starsolo_sum)




#































# quick dotplot
top10 <- top10_markers_ann %>% 
  select(1:7) 

top20_0 <- top20_markers_ann %>% 
  filter(cluster == 20)

DotPlot(seu_obj, features = top20_0$gene, 
        cols = c("blue", "red"),
         dot.scale = 8) +
  RotatedAxis()

VlnPlot(seu_obj, features = top10_0$gene)

#DimPlot(seu_obj, reduction = "umap", label=T, split.by = "Samples")


VlnPlot(seu_obj, features = "G32588", split.by = "Samples")

# quick check viral transcripts

# get all viral gene names
viral <- rownames(test_seu) %>% 
  as.data.frame() %>% 
  filter(str_detect(., "ORF" )) %>% 
  filter(!str_detect(., "\\." )) %>%
  rename(gene = ".")

# quick dotplot
DotPlot(test_seu, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) +
  theme(axis.text.x= element_text(angle=90, hjust=1))+
  coord_flip()

table(seu_obj$Samples)

#uninf control
Uninf_ctrl <- subset(seu_obj, subset = Samples == "Uninfected")
DotPlot(Uninf_ctrl, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


Homog_ctrl <- subset(seu_obj, subset = Samples == "Homogenate")
DotPlot(Homog_ctrl, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

hpi_6A <- subset(seu_obj, subset = Samples == "6-hpiA")
DotPlot(hpi_6A, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


hpi_6D <- subset(seu_obj, subset = Samples == "6-hpiD")
DotPlot(hpi_6D, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

hpi_24A <- subset(seu_obj, subset = Samples == "24-hpiA")
DotPlot(hpi_24A, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Gene')+
  ggtitle("viral transcripts: 24hpi-A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=10)) 

hpi_24J <- subset(seu_obj, subset = Samples == "24-hpiJ")
DotPlot(hpi_24J, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + 
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Gene')+
  ggtitle("viral transcripts: 24hpi-J")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=10)) 


hpi_72J <- subset(seu_obj, subset = Samples == "72-hpiJ")
DotPlot(hpi_72J, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) +
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Viral Gene')+
  ggtitle("Dotplot: Expression of viral transcripts in sample 72hpi-J")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=12)) 


# no of cells per cluster
table(hpi_72J@meta.data$seurat_clusters)

hpi_96E <- subset(seu_obj, subset = Samples == "96-hpiE")
DotPlot(hpi_96E, features = viral$gene, cols = c("blue", "red"), dot.scale = 8) + 
  coord_flip()+
  RotatedAxis()+
  ylab('Cluster') +  xlab('Viral Gene')+
  ggtitle("Dotplot: Expression of viral transcripts in sample 96hpi-E")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=12)) 

table(hpi_96E@meta.data$seurat_clusters)

# or create manually 
viral_transcripts_select <- as.data.frame(paste("ORF", 
                                                c(4,6,10,11,22,23,26,27,28,29,
                                                    41,44,47,49,51,53,56,57,66,71,76,77,
                                                    78,80,84,85,88,90,98,100,102,104,106,107,
                                                    108,110,112,113,116,117,119,121,122,124),
                                                sep = "")) %>% 
  rename_with(~ "gene")

DotPlot(seu_obj, features = viral_transcripts_select$gene, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

#DotPlot(seu_obj, features = c("ORF6", "ORF123"), cols = c("blue", "red"), dot.scale = 8)


#quick plot mito genes
mt.genes = c("ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB")

DotPlot(seu_obj, features = mt.genes, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

VlnPlot(seu_obj, features=mt.genes)


#check genes I was getting earlier in infection stages
my_genes <- c("G8077", "G8546", "G11423", "G4113", "G31799", "G15184", "G6896", "G2038")
DotPlot(seu_obj, features = my_genes, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
VlnPlot(seu_obj, features = my_genes)


my_genes2 <- c("G960", "G25637", "G1399", "G3128", "G1527")
DotPlot(seu_obj, features = my_genes2, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
VlnPlot(seu_obj, features = my_genes2)

my_genes3 <- c("G21268", "G13140", "G1208","G25110", "G19423")
DotPlot(seu_obj, features = my_genes3, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
VlnPlot(seu_obj, features = my_genes3)


#
#
#

# visualise cells by cluster to find infection stage-specific clusters, bubble plot
x <- prop.table(table(seu_obj@meta.data$Samples, seu_obj@meta.data$seurat_clusters), margin = 2)
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
  labs(title = "Distribution of cell counts by cluster across samples", x = "Sample", y = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



table(seu_obj@meta.data$seurat_clusters)
#save rds
saveRDS(seu_obj, file = "seu_obj_umap30x6.rds")
#......................................








# Get number of cells per cluster per group ID
#seu_obj <- Seurat::StashIdent(object = seu_obj, save.name = "my.clusters")
#bygroup<-table(seu_obj@meta.data$my.clusters, seu_obj@meta.data$sample)
#bysample<- table(hk@meta.data$my.clusters, seu_obj@meta.data$sample)

#to write the tables to an excel
#write.csv(bygroup, "bygroups.csv", row.names = T)
#write.csv(bysample, "M:/seu_obj/Single Cell Nuclei/bysample.csv", row.names = T)


# seu_objharmony <- seu_obj
# #with harmony  on samples. Harmony used to deal with batch effects and helps integrate the variation seen with various factors. 
# seu_objharmony <- NormalizeData(seu_objharmony, normalization.method = "LogNormalize", scale.factor = 10000)
# seu_objharmony <- FindVariableFeatures(seu_objharmony, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(seu_objharmony)
# seu_objharmony <- ScaleData(seu_objharmony, features = all.genes)
# seu_objharmony <- RunPCA(seu_objharmony, features = VariableFeatures(object = seu_objharmony))
# seu_objharmony<- RunHarmony(seu_objharmony, group.by.vars = "sample",plot_convergence = TRUE)
# ElbowPlot(seu_objharmony, ndims=50)
# seu_objharmony <- RunUMAP(seu_objharmony,reduction = "harmony",dims = 1:30)
# seu_objharmony<- FindNeighbors(seu_objharmony,reduction = "harmony", dims = 1:30)
# seu_objharmony<- FindClusters(seu_objharmony, resolution = 0.6)
# DimPlot(seu_objharmony, reduction = "umap", label=T)


# seu_objharmony <- Seurat::StashIdent(object = seu_objharmony, save.name = "my.clusters")
# bygroup<-table(seu_objharmony@meta.data$my.clusters, seu_objharmony@meta.data$groups)
# bysample<- table(seu_objharmony@meta.data$my.clusters, seu_objharmony@meta.data$sample)
# 
# write.csv(bygroup, "M:/seu_obj/Single Cell Nuclei/bygroups_harmonysample.csv", row.names = T)
# write.csv(bysample, "M:/seu_obj/Single Cell Nuclei/bysample_harmonysample.csv", row.names = T)
# 
# 
# p1<-DimPlot(seu_objharmony, reduction = "umap", label=T)
# p2<-DimPlot(seu_objharmony, reduction = "umap", label=F,group.by="sample")
# p3<-DimPlot(seu_objharmony, reduction = "umap", label=F,group.by="groups")
# ggarrange(p1,p2,p3, ncol = 3, nrow = 1)
# 
# 
# seu_obj1<-seu_objharmony
# new.cluster.ids <- c("Muscle ? 1","Cuticle 1","2","3","4","Cuticle 2","6","Cuticle 3","Muscle 1","Muscle 2","Cuticle 4","Chitins 1","12","Muscle ? 2","14","15","Chromatophores","Chitins 2","18","19","20","21")
# names(new.cluster.ids) <- levels(seu_obj)
# seu_obj <- RenameIdents(seu_obj, new.cluster.ids)
# p<-DimPlot(seu_obj, reduction = "umap", label=T,label.size = 4) + NoLegend()
# ggsave("M:/seu_obj/Single Cell Nuclei/analysis/Figures/umap.pdf", p, width=27, height=15, 
#    device = "pdf", units = "cm")  
# 
# 
# saveRDS(seu_objharmony, file = "M:/seu_obj/Single Cell Nuclei/harmony_by_sample.rds")



#Dotplot code
# moultdb<-c("LOC113804419",
# "LOC113826024",
# "LOC113820103",
# "LOC113800214",
# "LOC113814499",
# "LOC113814856",
# "LOC113827159",
# "LOC113801169",
# "LOC113828208",
# "LOC113804417",
# "LOC113800218",
# "LOC113823986",
# "LOC113820037",
# "LOC113824971")
# DotPlot(seu_obj,cols = c("white","darkorchid4"), col.min = 0.05, col.max = 3,dot.min = 0.05,  dot.scale = 8, features=moultdb)  + 
# theme(axis.text.x = element_text(angle = 45,vjust=0.8,hjust=0.8))  + coord_flip()
# 
# 
# 
# #Featureplot
# FeaturePlot(seu_obj, features = "LOC113804419")


###############................. find marker genes ################################

seu_obj <- readRDS(file = "seu_obj_umap30x6.rds")
DimPlot(seu_obj, reduction = "umap", label=T, group.by = "sample")
DimPlot(seu_obj, reduction = "umap", label=T)


#change sample names in orig.ident
seu_obj$original <- seu_obj$sample

table(seu_obj@meta.data$original)

seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre1")] <- "Uninf Ctrl"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre2")] <- "Homog Ctrl"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre3")] <- "6hpi-A"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre4")] <- "6hpi-D"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre5")] <- "24hpi-A"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre6")] <- "24hpi-J"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre7")] <- "72hpi-J"
seu_obj@meta.data$original[which(seu_obj@meta.data$original == "Ambre8")] <- "96hpi-E"

table(seu_obj@meta.data$original)

names(x=seu_obj[[]])

seu_obj$sample <- NULL
seu_obj$sample <- seu_obj$original
seu_obj$original <- NULL
names(x=seu_obj[[]])






####save
saveRDS(seu_obj, file = "seu_obj_umap30x6_sample_renamed.rds")

#plotting for cluster1
top_markers_clust_1 <- top_markers_20 %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 1)

VlnPlot(seu_obj, features = c("G13140","G1208", "G1206", "G5443", "G12478", "G12697", "G19423"), ncol = 3)
FeaturePlot(seu_obj, features = c("G13140","G1208", "G1206", "G5443", "G12478", "G12697", "G19423"), ncol = 3)



#plotting for cluster10
top_markers_clust_10 <- top_markers_20 %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 10)

VlnPlot(seu_obj, features = top_markers_clust_10$gene)

VlnPlot(seu_obj, features = c("G8077","G8546", "G11423", "G2038", "G8557", "G26708", "G13185", "G1206"), ncol = 3)
FeaturePlot(seu_obj, features = c("G8077","G8546", "G11423", "G2038", "G8557", "G26708", "G13185", "G1206"), ncol = 3)

#Finding markers for clusters and then plotting the top 5 with hierachial clustering
#all_markers <- FindAllMarkers(object = seu_obj, min.pct = 0.25, log2fc.threshold = 0.5)
all_markers <- FindAllMarkers(object = seu_obj, min.pct = 0.25, log2fc.threshold = 0.5)

# filter out expressed >1.5
all_markers_filt <- all_markers %>%
  dplyr::filter(avg_log2FC >= 1.5) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "gene_id") %>% 
  separate(gene, into = c("x", "Symbol"), remove = FALSE)


# count marker genes per cluster
all_markers_filt %>% 
  count(cluster)

# marker for cluster 11 which is hpi96 enriched
all_markers_cluster_11 <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 11) 



cluster_11_ann <- left_join(all_markers_cluster_11, mag_all_genes, by = "Symbol")



#write the marker list to file
all_markers_filt_ann <- left_join(all_markers_filt, mag_all_genes, by = "Symbol")
write_tsv(all_markers_filt_ann, "all_markers_filt_ann.tsv")

#............................................................................



#get gene symbols from ncbi
# https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/29159/
# select all columns, check box next to Gene ID to select all
# download "ncbi_dataset_taxon_29159.tsv"

#load ncbi gene id table
features_ncbi <- read_tsv("ncbi_dataset_taxon_29159.tsv", col_names = T) %>% 
  dplyr::rename(gene_id = "Nomenclature ID" )

# format chr to match gtf/gff3 files, need later to getFasta
features_ncbi$Chromosomes <- case_when(
  features_ncbi$Chromosomes == "1" ~ "LR761634.1",
  features_ncbi$Chromosomes == "2" ~ "LR761635.1",
  features_ncbi$Chromosomes == "3" ~ "LR761636.1",
  features_ncbi$Chromosomes == "4" ~ "LR761637.1",
  features_ncbi$Chromosomes == "5" ~ "LR761638.1",
  features_ncbi$Chromosomes == "6" ~ "LR761639.1",
  features_ncbi$Chromosomes == "7" ~ "LR761640.1",
  features_ncbi$Chromosomes == "8" ~ "LR761641.1",
  features_ncbi$Chromosomes == "9" ~ "LR761642.1",
  features_ncbi$Chromosomes == "10" ~ "LR761643.1",
  features_ncbi$Chromosomes == "MT" ~ "MZ497416.1"
)
# load all_genes.csv from Parse output
#genes_parse <- read_csv("DGE_unfiltered/all_genes.csv")


all_markers_filt_geneID <- left_join(all_markers_filt, features_ncbi, by = "gene_id")

#write the marker list to file
write_tsv(all_markers_filt_geneID, "all_markers_filt_geneID.tsv")

###############################............................................


all_markers_filt_geneID %>% 
  select(1, 18, 21,22)


# wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz
# gunzip Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz 






# get lit of top 10 markers for each cluster
top_markers_10 <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 10, named_vector = FALSE,
                                  make_unique = TRUE, data_frame = TRUE) %>%
  separate(gene, into = c("x", "Symbol"), remove = FALSE)

top_markers_10_ann <- left_join(top_markers_10, mag_all_genes, by = "Symbol")
write_tsv(top_markers_10_ann, "top_markers_10_ann.tsv")
  


top_markers_clust_10 <- top_markers_20 %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 10)



plot1<- VlnPlot(seu_obj, features = top_markers_clust_10$gene)
FeaturePlot(seu_obj, features = c("G20535", "G6643"))


## get top 10

top_markers_10 <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 10, named_vector = FALSE,
                                      make_unique = TRUE, data_frame = TRUE)


top_10markers_clust_10 <- top_markers_10 %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 10)



VlnPlot(seu_obj, features = top_10markers_clust_10$gene, ncol = 3)
FeaturePlot(seu_obj, features = top_10markers_clust_10$gene, ncol = 3)

VlnPlot(seu_obj, features = "gene-LOC105348737")
FeaturePlot(seu_obj, features = "gene-LOC105348737")


write_tsv(top_markers_20, "top_markers_20.tsv")

Clustered_DotPlot(seurat_object = seu_obj, features = top_markers)



### get bed file coordinates for top_markers_20
top_genes <- top_markers_20$gene

read_table("/external2/Pooran/seurat_thomas_scripts/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3")

gff <- as.data.frame(rtracklayer::import("/external2/Pooran/seurat_thomas_scripts/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3"))

#gff_v2<- gff %>% 
 # filter(end - start >3)

#TSS is the first bit in 5'UTR, but we don't have gene ID for 5'UTR entry
#I can see 5'UTR and gene all start with same cooridnate, so will filter gene string
#and get TSS info from that
gff_v2 <- gff %>% 
  select(1:3, 5, 7, 10) %>% 
  filter(type == "gene") %>%
  separate(col = "ID", sep = ":", into = c("feature", "gene"), remove = FALSE) %>% 
  select(1:4, 8)
  # select(1:6,8) %>% 
  # separate(col = "gene_ID", sep = "\\.", into = c("gene", "b"), remove = FALSE) %>% 
  # select(1:4,6,8) %>% 
  # #mutate(length_cds = end-start) %>% 
  # filter(end - start >30) 

all_markers_test_seu_up_bed <- inner_join(gff_v2, all_markers_test_seu_up, by = "gene", relationship = "many-to-many") %>% 
  select(1:3,gene, avg_log2FC, strand, cluster) 


#save by cluster
library(dplyr)
library(purrr)

# Create the directory if it doesn't exist
output_dir <- "bed_clusters_up"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Split the data by the cluster column
clusters <- split(all_markers_test_seu_up_bed, all_markers_test_seu_up_bed$cluster)

# Save each cluster as a separate file in the directory
walk(names(clusters), function(cluster_name) {
  file_path <- file.path(output_dir, paste0("cluster_", gsub("\\s|\\(|\\)", "_", cluster_name), ".txt"))
  write.table(
    clusters[[cluster_name]], 
    file = file_path, 
    row.names = FALSE, 
    sep = "\t", 
    quote = FALSE
  )
})


#write_tsv(all_markers_test_seu_up_bed, "top_markers_20_bed.bed", col_names = FALSE)

## getFasta for top markes 20 bed, bash programming
#bedtools getfasta -name -s -fi Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa -bed top_markers_20_bed.bed -fo top_markers_20.fasta


#######################....... translate all three frames using biostrings
library(Biostrings)


# import fasta file into biostring
fa_seqs  <- readDNAStringSet("top_markers_20.fasta")

# get all three position (1,2,3) as start
fa_subseqs <- lapply(1:3, function(pos) 
  subseq(fa_seqs, start=pos))

# translate
aa_seqs_list <- lapply(fa_subseqs, translate)


# Combine translated sequences into a single AAStringSet
# Initialize an empty AAStringSet
combined_aa_seqs <- AAStringSet()

# Loop to combine sequences from each frame and filter out those with stop codons
for (i in seq_along(aa_seqs_list)) {
  # Name each sequence
  names(aa_seqs_list[[i]]) <- paste0(names(fa_seqs), "_Frame_", i)
  
  # Filter sequences to remove those containing stop codons
  valid_seqs <- aa_seqs_list[[i]][!grepl("\\*", aa_seqs_list[[i]])]
  
  # Combine filtered sequences
  combined_aa_seqs <- c(combined_aa_seqs, valid_seqs)
}

# Step 5: Export the filtered amino acid sequences to a FASTA file
writeXStringSet(combined_aa_seqs, filepath = "filtered_translated_sequences.fasta")

# upload to egnog database for orhtology finding
#http://eggnog-mapper.embl.de/


##########.............................................................................

# load the top 20 marker genes file
top_markers_20 <- read_tsv ("top_markers_20.tsv")

# load eggnog output file
eggnog <- read_tsv("eggnog_out_from_aa/MM_tnx2xt46.emapper.annotations.tsv", comment = "##")


eggnog <- eggnog %>%
  separate(`#query`, sep = "\\.", into = c("a", "b"), remove = FALSE) %>% 
  separate(`a`, sep = ":", into = c("c", "gene")) 

# left join
top_markers_20_ann <- dplyr::left_join(top_markers_20, eggnog, 
                                       by = "gene", relationship = "many-to-many")

write_tsv(top_markers_20_ann, "top_markers_20_ann.tsv")

# vlnplot for each cluster

top_markers_10 <- top_markers %>% 
  filter(cluster == 10)

VlnPlot(seu_obj, features = top_markers_10$gene)
FeaturePlot(seu_obj, features = top_markers_10$gene)



# #I have my genome annotation file that I then use to add gene names and other annotations self done from papers to the gene IDs
# de_file <- read.csv("M:/seu_obj/Single Cell Nuclei/analysis/Figures/Moultdb/all_markers.csv", header = TRUE)
# ensemble_attributes <- read.csv("M:/seu_obj/Single Cell Nuclei/analysis/Gao2017/seu_obj_genome_withGao_Function.csv", header = TRUE, ",")
# final_output <- de_file %>% left_join(ensemble_attributes, by =c("gene"="Gene_ID"))
# write.csv(final_output,"M:/seu_obj/Single Cell Nuclei/analysis/Figures/Moultdb/all_markers.csv", , na="") 


## cluster17
top_markers_clust_17 <- top_markers_20 %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 17)

VlnPlot(seu_obj, features = top_markers_clust_17$gene)
FeaturePlot(seu_obj, features = top_markers_clust_17$gene)

# save gene IDs
write.table( top_markers_clust_17$gene, "ids.txt",
           ,row.names=FALSE,sep="\t", quote = FALSE,
           col.names = F)


# run extract.sh to retrieve aa sequences
# . extract.sh

# cat aa sequences
# cat output.fasta
# run protein blast or do orhtodb search

#md <- as.data.frame(seu_obj@meta.data)



########## assign cell type identity to cluster names
library(Seurat)
library(tidyverse)

seu_obj <- readRDS("seu_obj_umap30x6_sample_renamed.rds")

DimPlot(seu_obj, reduction = "umap", label=T)

new.cluster.ids <- c(new.cluster.ids <- c("0", "1", "2", "3", "4", "5", 
                                          "6","7", "8", "9", "10", "11", 
                                          "12", "13", "14", "15", "16", 
                                          "17", "18","19", "20", "21", "Excitatory neurons",
                                          "23", "24", "25"))

names(new.cluster.ids) <- levels(seu_obj)

seu_obj <- RenameIdents(seu_obj, new.cluster.ids)

DimPlot(seu_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


all_markers <- read_tsv("all_markers_filt.tsv")
library(scCustomize)
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 10, named_vector = FALSE,
                                      make_unique = TRUE, data_frame = TRUE)

top_markers_0 <- top_markers %>% 
  filter(cluster == 0)

FeaturePlot(seu_obj, features = top_markers_0$gene)
