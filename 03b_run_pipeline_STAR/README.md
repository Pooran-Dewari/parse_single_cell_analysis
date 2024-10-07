Thomas Clarke et al have managed to bypass Parse's spipe mapping and instead use STARSOLO, this should help with mapping rates as we can tweak mapping parameters.
Doing fresh analysis now to include STARSOLO mapping and also including Aurelie's OSHV1 genome.

***


### Genome indexing
#### Prepare composite gtf file
Files used
Crassostrea_gigas_uk_roslin_v1.gtf
viv46-2-m_assembly_NR_final_ok.gtf  
First create gtf files if using gff files, then merge the oyster gtf with the viral genome!  
 ```
#convert to gtf
agat_convert_sp_gff2gtf.pl --gff viv46-2-m_assembly_NR_final_ok.gff -o viv46-2-m_assembly_NR_final_ok.gtf
#merge host and pathogen gtf
agat_sp_merge_annotations.pl --gff ../Crassostrea_gigas_uk_roslin_v1.gtf --gff viv46-2-m_assembly_NR_final_ok.gtf --out Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf
```
add this bit of info about chr length to the merged gtf file Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf file using nano command manually
```
##sequence-region   viv46-2-m_assembly_NR_final 1 186279
```

#### Prepare composite fasta file
Files used
Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa
viv46-2-m_assembly_NR_2017.fasta

```
cat ../Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa viv46-2-m_assembly_NR_2017.fasta > Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa
```
Create a dir and move the composite genome and gtf files  
- Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa
- Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf

#### Index genome
While in composite_genome dir, run the command below to index genome
```
STAR  --runMode genomeGenerate --runThreadN 4 \
 --genomeDir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_star \
 --genomeFastaFiles Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa \
 --sjdbGTFfile Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf
```

***


### Run STARSOLO

#### prepare whitelist files
Need to gather barcode whitelist files from Parse spipe output (combine output > process > barcode_data.csv).  
The barcode 1 is sample type, I have copied barcode 1 for my samples wells into whitelist1.txt  (don't just include all 96 wells unless you have samples in the entire plate!!)  
The barcode 2 and 3 (for cell identity) are same really, I have copied barcode 2 into whitelist2.txt

```
$ head whitelist1.txt
CATTCCTA
CTTCATCA
CCTATATC
ACATTTAC
ACTTAGCT
CCAATTCT
```
#### run starsolo
## Make sure no spaces or empty lines in the STAR code, it just doesn't like it and will through some error! Do not create temp dir beforehand!

```
#!/bin/bash

#$ -V -cwd
#$ -l rl9=false
#$ -l h_rt=16:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 12

. /etc/profile.d/modules.sh
module load roslin/star/2.7.10b

read1="/exports/eddie/scratch/pdewari/newvolume/A3/split-fq/A3_EKDL240002475-1A_223M7CLT4_1_group_oyster_R1.fastq.gz"
read2="/exports/eddie/scratch/pdewari/newvolume/A3/split-fq/A3_EKDL240002475-1A_223M7CLT4_1_group_oyster_R2.fastq.gz"
sample="3"
temp="temp"$sample
out="a"$sample"_results"

STAR \
--runThreadN 32 \
--genomeDir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_star \
--readFilesCommand zcat \
--readFilesIn $read1 $read2 \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
--runDirPerm All_RWX \
--soloType CB_UMI_Complex \
--soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
--soloUMIposition 0_0_0_9 \
--soloCBwhitelist /exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/composite_genome/whitelist2.txt \
/exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/composite_genome/whitelist2.txt \
/exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/composite_genome/whitelist1.txt \
--soloCBmatchWLtype EditDist_2 \
--soloBarcodeReadLength 0 \
--soloFeatures GeneFull Gene \
--soloCellReadStats Standard \
--outTmpDir /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star/$temp \
--limitBAMsortRAM 230854492160 \
--outFileNamePrefix /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_star/$out \
--soloMultiMappers Rescue EM Uniform \
--outFilterScoreMinOverLread 0.50 \
--outFilterMatchNminOverLread 0.50

```

***

## prepare files for Seurat loading
Seurat needs the following three files:  
- matrix.mtx,
- genes.tsv (or features.tsv),
- and barcodes.tsv

Starsolo outputs these files:  
1- barcodes.tsv    
2- features.tsv   
3- matrix.mtx    
4- UniqueAndMult-EM.mtx    
5- UniqueAndMult-Rescue.mtx    
6- UniqueAndMult-Uniform.mtx  


### collate outputs for sub-libs

The script below will collate all relevant output files from the 8 sub-libraries into star_outputs_collated dir.

 
```
#!/bin/bash

# Logfile to store actions
logfile="logfile.txt"

# Clear previous log entries
> $logfile

# Create the final output directory 'star_outputs_collated' in the current working directory
final_dir="./star_outputs_collated"
echo "Creating final output directory: $final_dir"
echo "Creating final output directory: $final_dir" >> $logfile
mkdir -p "$final_dir"

# Check if the final directory was created successfully
if [ $? -eq 0 ]; then
    echo "Directory $final_dir created successfully."
    echo "Directory $final_dir created successfully." >> $logfile
else
    echo "Failed to create directory $final_dir."
    echo "Failed to create directory $final_dir." >> $logfile
    exit 1
fi

# Loop through the directories a1_resultsSolo.out to a8_resultsSolo.out
for dir in a*_resultsSolo.out; do
    # Navigate to the GeneFull/raw subdirectory
    sub_dir="$dir/GeneFull/raw"

    # Extract the prefix (e.g., a1, a2, a3) from the directory name
    prefix=$(echo $dir | cut -d'_' -f1)

    # Create the new subdirectory inside the final output directory with the extracted prefix
    new_dir="$final_dir/${prefix}_subdir"

    echo "Creating directory: $new_dir"
    echo "Creating directory: $new_dir" >> $logfile
    mkdir -p "$new_dir"

    # Check if the directory was created successfully
    if [ $? -eq 0 ]; then
        echo "Directory $new_dir created successfully."
        echo "Directory $new_dir created successfully." >> $logfile
    else
        echo "Failed to create directory $new_dir."
        echo "Failed to create directory $new_dir." >> $logfile
        continue
    fi

    # Files to copy
    files_to_copy=("barcodes.tsv" "features.tsv")

    # Copy each file (barcodes.tsv, features.tsv)
    for file in "${files_to_copy[@]}"; do
        src_file="$sub_dir/$file"
        dest_file="$new_dir/$file"

        echo "Copying $src_file to $dest_file"
        echo "Copying $src_file to $dest_file" >> $logfile
        cp "$src_file" "$dest_file"

        # Check if the file copy was successful
        if [ $? -eq 0 ]; then
            echo "Copied $file successfully."
            echo "Copied $file successfully." >> $logfile
        else
            echo "Failed to copy $file."
            echo "Failed to copy $file." >> $logfile
        fi
    done

    # Copy and rename UniqueAndMult-Rescue.mtx to matrix.mtx
    src_matrix="$sub_dir/UniqueAndMult-Rescue.mtx"
    dest_matrix="$new_dir/matrix.mtx"

    echo "Copying and renaming $src_matrix to $dest_matrix"
    echo "Copying and renaming $src_matrix to $dest_matrix" >> $logfile
    cp "$src_matrix" "$dest_matrix"

    # Check if the matrix file copy was successful
    if [ $? -eq 0 ]; then
        echo "Copied and renamed UniqueAndMult-Rescue.mtx to matrix.mtx successfully."
        echo "Copied and renamed UniqueAndMult-Rescue.mtx to matrix.mtx successfully." >> $logfile
    else
        echo "Failed to copy and rename UniqueAndMult-Rescue.mtx."
        echo "Failed to copy and rename UniqueAndMult-Rescue.mtx." >> $logfile
    fi

    echo "--------------------------------" >> $logfile
done

echo "Script execution complete. Check $logfile for details."
```

### move files to server and prepare final barcode file for each sub-lib

To prepare the final barcodes.tsv file, we need to extract well information for each sample barcode (i.e. barcode 1) from Parse output (combine output > process > barcode_data.csv), and then merge them with starsolo output file barcode.tsv, see the R script below:  

Note: this is for sub-lib 1, we will have to create barcode.tsv for each sub-lib; we will create individual Seurat objects for these sub-lib and can later merge them.  

```
library(tidyverse)
setwd("/home/pooran/Documents/seurat_thomas_scripts/barcode")

# load parse barcodes
parse <- read_csv("barcode_data.csv", comment = "#")


# tweak Parse barcode df to your experimental design
parse_experiment <- parse %>% 
  mutate(barcode = case_when(
    row_number() <= 192 ~ "barcode1",    # Rows 1-192
    row_number() <= 288 ~ "barcode2",    # Rows 193-288
    TRUE               ~ "barcode3"      # Rows 289-384
  )) %>% 
  mutate(sample = case_when(
    barcode == "barcode1" & well %in% c("C9", "C10")  ~ "Uninfected",
    barcode == "barcode1" & well %in% c("C11", "C12") ~ "Homogenate",
    barcode == "barcode1" & well %in% c("D1", "D2")   ~ "6h-hpiA",
    barcode == "barcode1" & well %in% c("D3", "D4")   ~ "6-hpiD",
    barcode == "barcode1" & well %in% c("D5", "D6")   ~ "24-hpiA",
    barcode == "barcode1" & well %in% c("D7", "D8")   ~ "24-hpiJ",
    barcode == "barcode1" & well %in% c("D9", "D10")  ~ "72-hpiJ",
    barcode == "barcode1" & well %in% c("D11", "D12") ~ "96-hpiE",
    TRUE ~ NA_character_  # If none of the conditions match, assign NA
  )) %>% 
  # my samples are in well C9 to D12
  filter(barcode == "barcode1") %>% 
  #filter(grepl("C(9|10|11|12)|D(1[0-2]?|[1-9])", well)) %>% 
  filter(!is.na(sample))


# load starsolo barcodes
starsolo <- read_delim("barcodes_a1.tsv", col_names = F, delim = "_") %>% 
  rename(
    bc2 = X1,
    bc3 = X2,
    bc1 = X3
  )


# intersect starsolo df with parse to get sample metadata
merged_bc <- inner_join(starsolo, parse_experiment, by = c("bc1" = "sequence"), keep = T) %>% 
  select(barcode=sequence, well, sample) %>% 
  distinct()

# write file
write_tsv(merged_bc, "barcodes_a1_seurat.tsv", col_names = F)
```

final barcode.tsv will have this format:  
```
ACATTCAT D10   72-hpiJ   
ACTGTGGG C12   Homogenate
ATCCTTAC C10   Uninfected
ATTCTGTC D5    24-hpiA   
CACCTTTA D2    6h-hpiA   
CACTTTCA D1    6h-hpiA
```
We have four matrix files from starsolo run in our case, could use any, but Thomas advised should try UniqueAndMult-Rescue.mtx first. We will have to rename this file as matrix.mtx.

