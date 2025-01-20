#aim: Extract TSS for gene hits using gff3 file

#step1: load gff3 file and convert to df for joining later
read_table("/external2/Pooran/seurat_thomas_scripts/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3")

gff <- as.data.frame(rtracklayer::import("/external2/Pooran/seurat_thomas_scripts/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3"))

#TSS is the first bit in 5'UTR, but we don't have gene ID for 5'UTR entry
#I can see 5'UTR and gene all start with same cooridnate, so will filter gene string
#and get TSS info from that
gff_v2 <- gff %>% 
  select(1:3, 5, 7, 10) %>% 
  filter(type == "gene") %>%
  separate(col = "ID", sep = ":", into = c("feature", "gene"), remove = FALSE) %>% 
  select(1:4, 8)

#step2: prepare gene hits df
all_markers_test_seu <- FindAllMarkers(object = test_seu, min.pct = 0.20, log2fc.threshold = 0.25)

all_markers_test_seu_up <- all_markers_test_seu  %>% 
  filter(avg_log2FC > 0.5) %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(!str_detect(gene, "MZ" ))

#step3: join tables to extract TSS for gene hits
all_markers_test_seu_up_bed <- inner_join(gff_v2, all_markers_test_seu_up, by = "gene", relationship = "many-to-many") %>% 
  select(1:3,gene, avg_log2FC, strand, cluster) 


#step4: save as bed file by cluster into a dir

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
    quote = FALSE
  )
})
