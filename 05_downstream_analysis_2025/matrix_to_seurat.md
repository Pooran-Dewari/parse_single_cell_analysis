**Until now:** Single-nucleus Parse datasets were mapped to C gigas genome using starsolo. We will be using EM matrix output files  for downstream analysis.  

# 1: Prepare files for seurat

## 1a: Directory structure of starsolo outputs

tree -d -L 1 /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict  
├── a1_resultsSolo.out  
├── a2_resultsSolo.out  
├── a3_resultsSolo.out  
├── a4_resultsSolo.out  
├── a5_resultsSolo.out  
├── a6_resultsSolo.out  
├── a7_resultsSolo.out  
├── a8_resultsSolo.out  

Directory structure of one of the Solo.out is as follows:  

tree /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/a2_resultsSolo.out/GeneFull/raw  

├── **barcodes.tsv -> use this**  
├── **features.tsv -> use this**  
├── matrix.mtx  
├── **UniqueAndMult-EM.mtx -> use this**  
├── UniqueAndMult-Rescue.mtx  
└── UniqueAndMult-Uniform.mtx  

we will use barcodes.tsv,  features.tsv, and  UniqueAndMult-EM.mtx  for downstream analysis.   

## 1b: Collate the above three files for each sub-library into one dir
The attached script 'star_outputs_unique_multi_EM_collate_gz.sh' will look into each of the eight Solo.out > GeneFull/raw dirs and copy the above three files into a new dir.    
Note: To comply with seurat requirements, we will rename UniqueAndMult-EM.mtx --> matrix.mtx, and gzip all three files after they have been copied.  
Run the bash script to collate these files, make sure you have all 8 sub-libs Solo.out dirs in the working dir.  

Directory structure of collated files is as follows:

tree /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_EM_collated  
/exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star_strict/star_outputs_Uni_Mult_EM_collated  
├── a1_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
├── a2_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
├── a3_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
├── a4_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
├── a5_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
├── a6_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
├── a7_subdir  
│   ├── barcodes.tsv.gz  
│   ├── features.tsv.gz  
│   └── matrix.mtx.gz  
└── a8_subdir  
    ├── barcodes.tsv.gz  
    ├── features.tsv.gz  
    └── matrix.mtx.gz  
    

## 1c: Merge Parse hexamers and polyA captures
Parse uses a mix of hexamers and poly(A) to capture mRNAs. We need to assign each nucleus by merging reads that come up from these two sources but the same mRNA molecule.  


