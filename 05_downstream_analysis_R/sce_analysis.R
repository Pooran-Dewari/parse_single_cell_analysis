# CAUTION: work in progress, very dry code.

library(BiocParallel)
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(janitor)
library(scran)
library(sctransform)
library(PCAtools)
library(pheatmap)
library(bluster)


# first set the seed to avoid shriek for all ML operations; also multicore param
set.seed(123)
bp.params <- MulticoreParam(workers = 4)

# 1= Uninfected control 
# 2= Homogenate control 
# 3= 6h post-infection (hpi) A 
# 4= 6 hpi D 
# 5= 24 hpi A 
# 6= 24 hpi J
# 7= 72 hpi J 
# 8= 96 hpi E
################################################################################





###### dirty way to add gene info into features.tsv file

# first load the original features.tsv that came from Parse outputs plus some awk
features_gigas <- read_tsv("features_gigas/features.tsv", col_names = F) %>%
  dplyr::rename(ID = X1)

# now load the gene to ID file that I downloaded from ncbi
gene_to_id_gigas <- read_tsv("features_gigas/gene_name_ncbi_dataset_c_gigas.tsv") %>%
  dplyr::select(1:3, 6, 11, 10, 14, 15) %>%
  clean_names() %>%
  dplyr::rename(ID = nomenclature_id)


# format chr to match gtf/gff3 files
gene_to_id_gigas$chromosomes <- case_when(
  gene_to_id_gigas$chromosomes == "1" ~ "LR761634.1",
  gene_to_id_gigas$chromosomes == "2" ~ "LR761635.1",
  gene_to_id_gigas$chromosomes == "3" ~ "LR761636.1",
  gene_to_id_gigas$chromosomes == "4" ~ "LR761637.1",
  gene_to_id_gigas$chromosomes == "5" ~ "LR761638.1",
  gene_to_id_gigas$chromosomes == "6" ~ "LR761639.1",
  gene_to_id_gigas$chromosomes == "7" ~ "LR761640.1",
  gene_to_id_gigas$chromosomes == "8" ~ "LR761641.1",
  gene_to_id_gigas$chromosomes == "9" ~ "LR761642.1",
  gene_to_id_gigas$chromosomes == "10" ~ "LR761643.1",
  gene_to_id_gigas$chromosomes == "MT" ~ "MZ497416.1"
) 
  
# join tables
new_features <- left_join(features_gigas, gene_to_id_gigas, by = "ID") %>%
  dplyr::select(1,5,3)

# write features.tsv, this will be used by read10xCounts
write_tsv(new_features, "features.tsv", col_names = F)

rm(features_gigas, gene_to_id_gigas, new_features)

#################################################

sample.path <- "/home/pdewari/Documents/parse/ambre_combined/DGE_filtered/for_sce"

samplesheet <- read_csv("../cell_metadata.csv") %>%
  dplyr::rename(Barcode = bc_wells) %>%
  dplyr::rename(Sample = sample) %>%
  dplyr::select(1:2) %>%
  dplyr::mutate(SampleGroup = Sample)


sce <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)
sce

dim(counts(sce))

counts(sce)[1:10, 1:10]

rowData(sce)

rownames(counts(sce))[1:5]

colData(sce)[1:5, 1]

#
genesPerCell <- colSums(counts(sce) > 0)
plot(density(genesPerCell), main="", xlab="Genes per cell")


#
sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), samplesheet, by="Barcode", sort=FALSE)
rownames(colData(sce)) <- sce$Barcode

#
plot(rowSums(counts(sce)) / rowSums(counts(sce) > 0),
     rowMeans(counts(sce) > 0),
     log = "x",
     xlab="Mean UMIs per cell",
     ylab="proportion of cells expressing the gene"
)


#
rel_expression <- t( t(counts(sce)) / colSums(counts(sce))) * 100
rownames(rel_expression) <- rowData(sce)$ID #change ID to Symbol if need be
most_expressed <- sort(rowSums( rel_expression ), decreasing = T)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)


#
colData(sce)

#
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)

#
sce <- sce[detected_genes,]

mt_genes = c("ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB")

is.mito <- which(rowData(sce)$ID %in% mt_genes)

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)
colData(sce)


plotColData(sce, x="Sample.y", y="subsets_Mito_percent") + 
  scale_y_log10() + 
  ggtitle("Mito present")


colData(sce) %>% 
  as.data.frame()


#
plotColData(sce, x="Sample.y", y="sum") + 
  scale_y_log10() + 
  ggtitle("Total count")

#
plotColData(sce, x="Sample.y", y="detected") + 
  scale_y_log10() + 
  ggtitle("Total count")

#
colData(sce) %>% 
  as.data.frame() %>% 
  arrange(subsets_Mito_percent) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mito_percent > 5)) + 
  facet_wrap(vars(SampleGroup))


###

cell_qc_results <- quickPerCellQC(colData(sce), sub.fields = TRUE)

cell_qc_results %>%
  as.data.frame() %>% 
  mutate(SampleName=colData(sce)$Sample.y) %>% 
  group_by(SampleName) %>%
  summarise(across(where(is.logical), sum))

###
sce$low_lib_size <- cell_qc_results$low_lib_size
sce$low_n_features <- cell_qc_results$low_n_features
sce$high_Mito_percent <- cell_qc_results$high_subsets_Mito_percent
sce$discard <- cell_qc_results$discard
sce.filtered <- sce[, !sce$discard]

# save normalised object
saveRDS(sce.filtered, "file1_filtered.rds")
#########################################################################




######################## use sce.filtered only from here onwards! ############

sce.filtered <- readRDS("file1_filtered.rds")

table(sce.filtered$Sample.y)

oneSamTab <- colData(sce.filtered) %>% 
  as.data.frame() %>% 
  filter( Sample.y == "Ambre1") %>% 
  dplyr::select(Sample.y,Barcode, sum) %>% 
  mutate(cell_num = 1:n())

p_before_nom <- ggplot(data=oneSamTab, aes(x=cell_num, y=sum)) +
  geom_bar(stat = 'identity') +
  ylim(0, 4000)+
  labs( x= 'Cell Index',
        y='Cell UMI counts',
        title = "Before Normalization" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )
p_before_nom


#...


clust <- quickCluster(sce.filtered, BPPARAM=bp.params)
table(clust)

sce.filtered <- computePooledFactors(sce.filtered,
                            clusters = clust,
                            min.mean = 0.1,
                            BPPARAM = bp.params)

deconv.sf <- sizeFactors(sce.filtered)
summary(deconv.sf)

assayNames(sce.filtered)
sce.filtered <- logNormCounts(sce.filtered)
assayNames(sce.filtered)

##
norm_counts <- logNormCounts(sce.filtered,transform='none' ) %>% 
  assay('normcounts') %>% 
  as.matrix() %>% 
  colSums()

norm_counts <- tibble(Barcode=names(norm_counts),
                      normCounts = log2(norm_counts)
)

norm_counts <- inner_join(norm_counts, oneSamTab, by='Barcode')


p_after_norm <- ggplot(data=norm_counts, aes(x=cell_num, y=normCounts)) +
  geom_bar(stat = 'identity') +
  labs( x= 'Cell Index',
        y='Normalized Cell UMI counts',
        title = "PBMMC_1:After Normalization" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )

p_before_nom + p_after_norm

# save normalised object
saveRDS(sce.filtered, "file2_normalised.rds")

##########


sce.filtered <- readRDS("file2_normalised.rds")


rownames(sce.filtered)

gene_var <- modelGeneVar(sce.filtered)

gene_var

hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs)
plotExpression(sce.filtered, features = hvgs[1:20], point_alpha = 0.05)

sce.filtered <- runPCA(sce.filtered, subset_row = hvgs)
sce.filtered
reducedDim(sce.filtered, "PCA")[1:10, 1:5]

percent.var <- attr(reducedDim(sce.filtered), "percentVar")
plot(percent.var, log = "y", xlab = "PC", ylab = "Variance explained (%)")

# view pca plot
plotReducedDim(sce.filtered, dimred = "PCA", 
               colour_by = "Sample.y")

#table(percent.var > 1)
#chosen_elbow <- findElbowPoint(percent.var)
#chosen_elbow
#plot(percent.var)
#abline(v=chosen_elbow, col="dodgerblue")

########

sce.filtered <- runUMAP(sce.filtered)

plotUMAP(sce.filtered, colour_by = "Sample.y")

ggcells(sce.filtered, aes(x = UMAP.1, 
                          y = UMAP.2, 
                          colour = Sample.y)) +
  geom_point()


# run with neighbors set to 50
sce.filtered <- runUMAP(sce.filtered, 
               name = "UMAP_neighbors50",
               dimred = "PCA",
               n_neighbors = 50)

ggcells(sce.filtered, aes(x = UMAP_neighbors50.1, 
                          y = UMAP_neighbors50.2, 
                          colour = Sample.y)) +
  geom_point()

# run with neighbors set to 5
sce.filtered <- runUMAP(sce.filtered, 
               name = "UMAP_neighbors5",
               dimred = "PCA",
               n_neighbors = 5)

ggcells(sce.filtered, aes(x = UMAP_neighbors5.1, y = UMAP_neighbors5.2, 
                 colour = Sample.y)) +
  geom_point() +
  labs(title = "Neighbours = 5")

# run with neighbors set to 500
sce.filtered <- runUMAP(sce.filtered, 
               name = "UMAP_neighbors500",
               dimred = "PCA",
               n_neighbors = 500)

ggcells(sce.filtered, aes(x = UMAP_neighbors500.1, y = UMAP_neighbors500.2, 
                 colour = Sample.y)) +
  geom_point() +
  labs(title = "Neighbours = 500")


# save normalised object
saveRDS(sce.filtered, "file3_dimRed.rds")



















##########

library(igraph)
sce.filtered <- readRDS("file3_dimRed.rds")




sce.filtered$leiden40 <- clusterCells(sce.filtered, 
                             use.dimred = "PCA", 
                             BLUSPARAM = SNNGraphParam(k = 40, 
                                                       cluster.fun = "leiden"))

table(sce.filtered$leiden40)

plotReducedDim(sce.filtered, 
               dimred = "UMAP",
               colour_by = "leiden40", 
               text_by = "leiden40")


#####

sil.approx <- approxSilhouette(reducedDim(sce.filtered, "PCA"),
                               clusters=sce.filtered$leiden40)
sil.approx


plotSilBeeswarm <- function(silDat){
  silTab <- silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor())
  
  plt <- silTab %>% 
    ggplot(aes(x=cluster, y=width, colour=closestCluster)) +
    ggbeeswarm::geom_quasirandom(method="smiley", alpha=0.6) +
    theme_bw()
  
  plt
}


p1 <- plotSilBeeswarm(sil.approx)
p2 <- plotReducedDim(sce.filtered, 
                     dimred = "PCA", 
                     colour_by="leiden40", 
                     text_by = "leiden40")
p1 + p2


### overlay clusters on umap

colLabels(sce.filtered) <- sce.filtered$leiden40

plotReducedDim(sce.filtered, 
               dimred = "UMAP",
               colour_by = "label", 
               text_by = "label") +
  ggtitle("Leiden k=40 clusters")


#G748 interferon regulatory factor 1
colLabels(sce.filtered) <- sce.filtered$leiden40

plotReducedDim(sce.filtered, 
               dimred = "UMAP",
               colour_by = "label" , 
               text_by = "label")


####################################

# colour by a gene expression, facet_wrap by sample names
ggcells(sce.filtered, aes(x = UMAP.1, y = UMAP.2, colour = G748)) +
  geom_point(size = 0.5) +
  facet_wrap(~ Sample.y)+
  scale_colour_gradient(low ="gray", high = "black")+
  theme_bw()


# label by cluster no.
ggcells(sce.filtered, aes(x = UMAP.1, y = UMAP.2, colour = label)) +
  geom_point(size = 1) +
  facet_wrap(~ Sample.y)

# get list of marker genes; this one is list of lists, total 20 lists
markers <- scoreMarkers(sce.filtered, 
                        groups = sce.filtered$label, 
                        block = sce.filtered$Sampl.y)

c18_markers <- as.data.frame(markers[["18"]])




################# find the gene symbol for gene ID enriched in cluster markers

# first get significatnly-enriched cluster marker genes
cluster_markers_refined <- as.data.frame(markers[["18"]]) %>% 
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 8) %>%
  filter(mean.logFC.cohen > 1) %>% 
  arrange(rank.logFC.cohen)

# get gene ID names as df for inner join later
cluster_markers_id <- as.data.frame(rownames(cluster_markers_refined)) %>% 
  rename_with(.cols = 1, ~"ID")

# option1- intersect with full gene list
intersect(cluster_markers_id$ID, gene_to_id_gigas$ID)

# option2- inner join to find the common elements to retrieve gene info
cluster_markers_coord <- inner_join(cluster_markers_id, 
                                    gene_to_id_gigas, by = "ID") %>% 
  dplyr::select(chromosomes, 
                annotation_genomic_range_start,
                annotation_genomic_range_stop,
                ID)

write_tsv(cluster_markers_coord, "/home/pdewari/Documents/parse/cluster_markers_coord.bed", col_names = F)

# get fasta sequences using bedtools (bash)
# bedtools getfasta -fi Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa -bed cluster_markers_coord.bed -fo cluster_markers_coord.bed.out


#cluster_markers_id$ID %in% gene_to_id_gigas$ID
#which(cluster_markers_id$ID %in% gene_to_id_gigas$ID)

# re-label the cells - original cluster in parenthesis
levels(colLabels(sce.filtered)) <- c("C1 (1)", "C2 (2)", 
                            "C3 (3)", "C4 (4)",
                            "C5 (5)", "C6 (6)",
                            "C7 (7)", "C8 (8)",
                            "C9 (9)","Endothelial? (10)",
                            "Muscle (11)", "Neurons? (12)",
                            "C13 (13)", "Gastric? (14)",
                            "C15 (15)", "C16 (16)",
                            "C17 (17)", "C18 (18)",
                            "C19 (19)", "C20 (20)"
                            )

plotReducedDim(sce.filtered, 
               dimred = "UMAP",
               colour_by = "label" , 
               text_by = "label")

# label by cluster no.
ggcells(sce.filtered, aes(x = UMAP.1, y = UMAP.2, colour = label)) +
  geom_point(size = 1) +
  facet_wrap(~ Sample.y)


# subset samples
sce.filtered_ambre1 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre1"]

colData(sce.filtered)

colData(sce.filtered_ambre1)

# plot UMAP for ambre1, add title
plotReducedDim(sce.filtered_ambre1, 
               dimred = "UMAP",
                     colour_by="label", 
                     text_by = "label") +
  ggplot2::ggtitle(label ="UMAP Ambre1")

# ggcells(sce.filtered_ambre1, 
#         aes(x = UMAP.1,
#                  y = UMAP.2, 
#                  colour = label)) +
#   geom_point() +
#   labs(title = "UMAP_on_PCA")
##

plotHeatmap(sce.filtered, 
            features = rownames(cluster_markers_refined),
            order_columns_by = c("label", "Sample.y"))

# scaled heatmap (z-scores)
plotHeatmap(sce.filtered, 
            features = rownames(cluster_markers_refined),
            order_columns_by = c("label", "Sample.y"),
            scale = TRUE, 
            center = TRUE,
            zlim = c(-3, 3))

plotGroupedHeatmap(sce.filtered, 
                   features = rownames(cluster_markers_refined),
                   group = "label",
                   block = "Sample.y",
                   scale = TRUE, 
                   center = TRUE,
                   zlim = c(-3, 3))


########################### grid plot


library(cowplot)

# more info about grid plotting umap here
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/scater/scater_02_dim_reduction.html

# 1= Uninfected control 
# 2= Homogenate control 
# 3= 6h post-infection (hpi) A 
# 4= 6 hpi D 
# 5= 24 hpi A 
# 6= 24 hpi J
# 7= 72 hpi J 
# 8= 96 hpi E


# compare homegenate control with 24 hr post infection
# subset samples

sce.filtered_ambre1 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre1"]
sce.filtered_ambre2 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre2"]
sce.filtered_ambre3 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre3"]
sce.filtered_ambre4 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre4"]
sce.filtered_ambre5 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre5"]
sce.filtered_ambre6 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre6"]
sce.filtered_ambre7 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre7"]
sce.filtered_ambre8 <- sce.filtered [, sce.filtered$SampleGroup =="Ambre8"]

# plot grid
# control vs 6 hpi
plot_grid(ncol = 4,
          
          plotReducedDim(sce.filtered_ambre1, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP uninfected control"),
          
          
          plotReducedDim(sce.filtered_ambre2, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP homogenate control"),
          
          plotReducedDim(sce.filtered_ambre3, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 6 hpi A"),
          
          
          plotReducedDim(sce.filtered_ambre4, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 6 hpi D")
)

###

#  6 hpi vs 24 hpi
plot_grid(ncol = 4,
          
          plotReducedDim(sce.filtered_ambre3, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 6 hpi"),
          
          
          plotReducedDim(sce.filtered_ambre4, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 6 hpi"),
          
          plotReducedDim(sce.filtered_ambre5, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 24 hpi"),
          
          
          plotReducedDim(sce.filtered_ambre6, 
                         dimred = "UMAP",
                         colour_by="label", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 24 hpi")
)


#
#  24 hpi vs 72 hpi vs 96 hpi
plot_grid(ncol = 4,
          
          plotReducedDim(sce.filtered_ambre5, 
                         dimred = "UMAP",
                         colour_by="G6643", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 24 hpi"),
          
          
          plotReducedDim(sce.filtered_ambre6, 
                         dimred = "UMAP",
                         colour_by="G6643", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 24 hpi"),
          
          plotReducedDim(sce.filtered_ambre7, 
                         dimred = "UMAP",
                         colour_by="G6643", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 72 hpi"),
          
          
          plotReducedDim(sce.filtered_ambre8, 
                         dimred = "UMAP",
                         colour_by="G6643", 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 96 hpi")
)


#################
#  24 hpi vs 72 hpi vs 96 hpi
# colour by a gene we know is highly enriched in cluster 18

markers <- scoreMarkers(sce.filtered, 
                        groups = sce.filtered$label, 
                        block = sce.filtered$Sampl.y)

cluster_markers_refined <- as.data.frame(markers[["19"]]) %>% 
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 8) %>%
  filter(mean.logFC.cohen > 1) %>% 
  arrange(rank.logFC.cohen)


# my_gene <- "G18691" # cluster 18, seen only in 96 hpi
# my_gene <- "G5032" #cluster 5
# my_gene <- "G1567" #cluster 8
# my_gene <- "G24446" #cluster 10

# my_gene <- "G22023" #cluster 12

plot_grid(ncol = 4,
          
          plotReducedDim(sce.filtered_ambre5, 
                         dimred = "UMAP",
                         colour_by = my_gene, 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 24 hpi")+
            scale_colour_gradient(low ="gray", high = "red"),
          
          
          plotReducedDim(sce.filtered_ambre6, 
                         dimred = "UMAP",
                         colour_by = my_gene, 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 24 hpi") +
            scale_colour_gradient(low ="gray", high = "red"),
          
          plotReducedDim(sce.filtered_ambre7, 
                         dimred = "UMAP",
                         colour_by = my_gene, 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 72 hpi") +
            scale_colour_gradient(low ="gray", high = "red"),
          
          
          plotReducedDim(sce.filtered_ambre8, 
                         dimred = "UMAP",
                         colour_by = my_gene, 
                         text_by = "label") +
            ggplot2::ggtitle(label ="UMAP 96 hpi") +
            scale_colour_gradient(low ="gray", high = "red")
)
