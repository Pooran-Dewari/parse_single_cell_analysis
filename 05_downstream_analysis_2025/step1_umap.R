# working on this script: 24 March 2025


######################### STEP1: loading data and filtering ###################
#load seurat objects, merge, filter out rRNA genes, discard high mito cells


#load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
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

head(seu_obj)
# Assign sample labels to a new 'sample' column based on the 'orig.ident' column
seu_obj@meta.data$sample <- orig_to_sample[seu_obj@meta.data$orig.ident]

#add condition column to metadata for grouping later
seu_obj$condition <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "exposed",
  seu_obj$sample %in% c("24-hpiA", "24-hpiJ") ~ "early",
  seu_obj$sample %in% c("72-hpiJ", "96-hpiE") ~ "late",
  TRUE ~ "unknown"  # Optional: Handles unexpected values
)


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

#  order the conditions
seu_obj$condition <- factor(seu_obj$condition, levels=c("control",
                                                     "exposed",
                                                     "early",
                                                     "late"))
table(seu_obj@meta.data$condition)


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
library(harmony)
library(scCustomize)

# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

#standard seurat workflow until elbow plot

# load the filtered seu obj
seu_obj <- read_rds("seu_obj_filt.rds")
head(seu_obj)
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
library(harmony)
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
DimHeatmap(seu_obj, dims = 1:30, cells = 500, balanced = TRUE)

#use this to assess number of PC's to use in downstream analysis
ElbowPlot(seu_obj, ndims=50)

#looking at the elbo plot and heatmap, dims 15,20, 25,30 look reasonable to test!

#test dims=15

# going ahead with dims = 20 and res 0.6

#run harmony
seu_obj<- RunHarmony(seu_obj, group.by.vars = "sample", plot_convergence = TRUE)

seu_obj
# An object of class Seurat 
# 33466 features across 23740 samples within 1 assay 
# Active assay: RNA (33466 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, harmony

seu_obj<- FindNeighbors(seu_obj, dims = 1:18, reduction = "harmony")
seu_obj<- FindClusters(seu_obj, resolution = 0.6)


seu_obj <- RunUMAP(seu_obj, dims = 1:18, reduction = "harmony")

#plot the UMAP
DimPlot(seu_obj, reduction = "umap", label = TRUE, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = T, label.size = 6) + 
  NoLegend()+
  ggtitle("Starsolo harmony umap: 18 dim x 6 res x 3k variable features")+
  theme(plot.title = element_text(hjust = 0.5))

seu_obj #should now show pca, harmony, and umap  in reduction layer

# An object of class Seurat 
# 33466 features across 23740 samples within 1 assay 
# Active assay: RNA (33466 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, harmony, umap

########### highlight clusters

# Define the directory name
dir_name <- "18dimx6resx3kFeatres_harmony_24march"

# Source the script and pass `dir_name` as an argument
source("highlight.R")

# Call the `run_script()` function with the `dir_name`
run_script(dir_name)

saveRDS(seu_obj, file = "seu_obj_umap_18d_6r_3kRes.rds")

rm(list=ls())
######
#####
############# harmony and UMAP done !!! ####################################
####
##
#
