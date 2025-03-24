library(harmony)
#load packages
library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)


# Set the working directory
setwd("/home/pooran/Documents/parse_2025/seurat_2025/")

# load seurat object that we had created in previous step with pca done!
seu_obj <- read_rds("seu_obj_umap_18d_6r_3kRes.rds")
seu_obj 
#Note: 
# Uninfected and Homogenate: control
# 6 hours: exposed
# 24 hours: early
# 72 and 96 hours: late


######### first assign the new clusters to our original marker list
# cluster 0:  
# cluster 1:   cluster 2,7, 15 infection response of original, not sure what this one is.
# cluster 2:   Gill_ciliary
# cluster 3:   Hepatopancreas
# cluster 4:   Gill_NEC
# cluster 5:   Gill_type1
# cluster 6:   Hyalinocytes
# cluster 7:   Haemocytes_type1
# cluster 8:   cluster 10 of original
# cluster 9:   cluster 12 of original
# cluster 10:  Mantle or vesicular cells cluster 7 of original
# cluster 11:  Immature_haemocytes
# cluster 12:  Macrophage_like
# cluster 13:  Adductor_muscle
# cluster 14:  Mantle
# cluster 15:  Digestive_gland
# cluster 16:  Gill_type3 of original
# cluster 17:  Small_granule_cells




new.cluster.ids <- c(new.cluster.ids <- c("cluster 0", "cluster 1", "gill_ciliary(2)",
                                          "hepatopancreas(3)", "gill_nec(4)", "gill_type1(5)",
                                          "hyalinocytes(6)", "haemocytes_type1(7)",
                                          "cluster 8", "cluster 9", "mantle/vesicular(10)",
                                          "immature_haemocytes(11)", "macrophage_like(12)",
                                          "adductor_muscle(13)", "mantle(14)", "digestive_gland(15)",
                                          "gill_type2(16)", "small_granules_cells(17)"))


names(new.cluster.ids) <- levels(seu_obj)

seu_obj <- RenameIdents(seu_obj, new.cluster.ids)
table(seu_obj$seurat_clusters)
#add the annotations to metadata column
seu_obj$seurat_annotations <- Idents(seu_obj)

DimPlot(seu_obj, reduction = "umap", label = T, 
        pt.size = 0.5, repel = T, label.box = T, 
        sizes.highlight = F, label.size = 6) + 
  #NoLegend()+
  ggtitle("Single-nucleus atlas of Paciifc oyster: harmony UMAP")+
  theme(plot.title = element_text(hjust = 0.5))

###################################








# main idea here is to do pairwise comparison of clusters across samples/conditions
# testing one cluster to begin with
# plotting cluster 7 markers from my previous analysis, pre harmony data to select a cluster
VlnPlot(seu_obj, c("G8077", "G8546", "G31799", "G24283", 
                   "G4113", "G11423", "G8071", "G15184",
                   "G6896", "G4729", "G26708", "G6764", "G8557"))

# based on this, extract cluster 1

# subset cluster 1 from the original seurat object
cl1 <- subset(x = seu_obj, idents ="1")
head(cl1)
Idents(cl1) #this is all '1' because we subset cluster 1

DimPlot(seu_obj, split.by="sample",reduction = "umap", label=T)
DimPlot(cl1, split.by="sample",reduction = "umap", label=T)


# re-normalise and scale data for subset cl1
cl1 <- NormalizeData(cl1, normalization.method = "LogNormalize", scale.factor = 10000)
cl1 <- FindVariableFeatures(cl1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cl1)
cl1 <- ScaleData(cl1, features = all.genes)
cl1 <- RunPCA(cl1, features = VariableFeatures(object = cl1))
cl1<- RunHarmony(cl1, group.by.vars = "sample",plot_convergence = TRUE)
ElbowPlot(cl1, ndims=50)
cl1 <- RunUMAP(cl1,reduction = "harmony",dims = 1:18)          
cl1 <- FindNeighbors(cl1,reduction = "harmony", dims = 1:18)
cl1 <- FindClusters(cl1, resolution = 0.6)
DimPlot(cl1, reduction = "umap", label=T, group.by = "sample")
#saveRDS(cl1, file = "cl1.rds")



cl1_bulk <- cl1
Idents(cl1_bulk) #Idents is 0-9 clusters that were assigned during FindClusters above
# Idents will be Levels: 0 1 2 3 4 5 6 7 8 9

cl1_bulk<- FindClusters(cl1_bulk, resolution = 0.01) #low resolution to make sure only one cluster pops up

DefaultAssay(cl1_bulk) #should be RNA anyway

cl1_bulk$sample.subtypes <- paste(Idents(cl1_bulk), cl1_bulk$condition, sep = "_")
#cl1_bulk$seurat_clusters <- Idents(cl1_bulk)
Idents(cl1_bulk) <- "sample.subtypes"
Idents(cl1_bulk)
head(cl1_bulk)
DimPlot(cl1_bulk, reduction = "umap", label=T)
DimPlot(cl1_bulk, reduction = "umap", label=T)

cl1_bulk_markers <- FindAllMarkers(cl1_bulk, only.pos = TRUE, min.pct = 0.25, log2fc.threshold = 0.5)

# annotate all_markers with ORSON
ORSON <- read_csv("ORSON_french_group_2024_biorxiv_suppl2.csv")

# now annotate with cg_science paper
cg_science <- read_tsv("Cg_gene_names.tsv", col_names = F) %>% 
  as.data.frame()

colnames(cg_science)[c(1, 2)] <- c("Gene.ID", "Gene.symbol.science")

all_markers_cl1_bulk_ann <- left_join(cl1_bulk_markers, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description) %>% 
  filter(avg_log2FC > 1 & p_val_adj < 0.10)


DimPlot(cl1_bulk, reduction = "umap", label=T)

#  order the conditions
Idents(cl1_bulk) <- factor(Idents(cl1_bulk), levels = c("0_control",
                                                        "0_exposed",
                                                        "0_early",
                                                        "0_late"))

# check condition-specific expression of cluster7 markers
VlnPlot(cl1_bulk, stack=T, c("G8077", "G8546", "G31799", "G24283", "G4113",  
                    "G11423", "G8071", "G15184", "G6896", "G4729"), flip = T)

VlnPlot(cl1_bulk, c("G26708", "G6764", "G8557", "G2038", "G4723",  
                    "G5363", "G11603", "G9349", "G13185", "G14351"))

VlnPlot(cl1_bulk, c("G1230", "G5730", "G23776", "G5731", "G5361",  
                    "G17069", "G20645", "G14487", "G31523", "G20683"))

VlnPlot(cl1_bulk, c("G26870", "G18691", "G1206", "G27529", "G18199",  
                    "G5430", "G2047", "G135", "G10173", "G16354"))

VlnPlot(cl1_bulk, c("G28262", "G10014", "G14767", "G2045", "G1208",  
                    "G2120", "G6852", "G17926", "G2118", "G2952"))

########################################################################

#### can do pairwise here if need be
contrl_vs_late <- subset(x = cl1_bulk, idents =c("0_control", "0_late"))
#no need to normalise, just find markers
contrl_vs_late_markers <- FindAllMarkers(contrl_vs_late, only.pos = TRUE, min.pct = 0.25, log2fc.threshold = 0.5)
#annotation
contrl_vs_late_markers_ann <- left_join(contrl_vs_late_markers, ORSON, by = c("gene"="Gene.ID")) %>% 
  left_join(., cg_science, by = c("gene"="Gene.ID")) %>% 
  relocate(Gene.symbol.science, .after = Sequence.Description) %>% 
  filter(avg_log2FC > 1 & p_val_adj < 0.10)
####
