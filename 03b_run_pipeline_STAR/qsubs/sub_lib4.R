library(Seurat)

setwd("/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_Resc_collated/fix_hex_oligdT") 

source("merge_parse_hexamers_polyA_captures.R") 

#update subdir
mat <- Read10X(data.dir ="/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_Resc_collated/a4_subdir")

merge <- merge_hexamer_polyA_columns(mtx = mat, kit = "WT_mega", version = "v2", bc_directory = "/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_Resc_collated/fix_hex_oligdT")

#update project
seu_obj <- CreateSeuratObject(counts = merge, project = "oyster4")

#update file
saveRDS(seu_obj, file = "seu_obj4.rds")
