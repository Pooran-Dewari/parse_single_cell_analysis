library(Seurat)

setwd("/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_Resc_collated/fix_hex_oligdT") 

source("merge_parse_hexamers_polyA_captures.R") 

mat1 <- Read10X(data.dir ="/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_Resc_collated/a1_subdir")

merge1 <- merge_hexamer_polyA_columns(mtx = mat1, kit = "WT_mega", version = "v2", bc_directory = "/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_Resc_collated/fix_hex_oligdT")

seu_obj1 <- CreateSeuratObject(counts = merge1, project = "oyster1")

saveRDS(seu_obj1, file = "seu_obj1.rds")
