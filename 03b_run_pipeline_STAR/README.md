Thomas Clarke et al have managed to bypass Parse's spipe mapping and instead use STARSOLO, this should help with mapping rates as we can tweak mapping parameters.
Doing fresh analysis now to include STARSOLO mapping and also including Aurelie's OSHV1 genome.

***
### Notes:   

Multi-gene reads
Multi-gene reads are concordant with (i.e. align equally well to) transcripts of two or more genes. One class of multi-gene read are those that map uniquely to a genomic region where two or more genes overlap. Another class are those reads that map to multiple loci in the genome, with each locus annotated to a different gene.    

Including multi-gene reads allows for more accurate gene quantification and, more importantly, enables detection of gene expression from certain classes of genes that are supported only by multi-gene reads, such as overlapping genes and highly similar paralog families.  

The multi-gene read recovery options are specified with --soloMultiMappers. Several algorithms are implemented:  

- **--soloMultiMappers Uniform** uniformly distributes the multi-gene UMIs to all genes in its gene set. Each gene gets a fractional count of 1/N_genes, where N_genes is the number of genes in the set. This is the simplest possible option, and it offers higher sensitivity for gene detection at the expense of lower precision.  --soloMultiMappers Uniform distributes the multi-gene UMIs proportionally to the number of unqiue UMIs per gene. UMIs that map to genes that are not supported by unique UMIs are distributed uniformly.  

- **--soloMultiMappers EM** uses Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their genes, taking into account other UMIs (both unique- and multi-gene) from the same cell (i.e. with the same CB). Expectation-Maximization (EM) algorithm is used to find the gene expression values that maximize the likelihood function. Recovering multi-gene reads via MLE-EM model was previously used to quantify transposable elements in bulk RNA-seq {TEtranscripts} and in scRNA-seq {Alevin; Kallisto-bustools}.  
 
- **--soloMultiMappers Rescue** distributes multi-gene UMIs to their gene set proportionally to the sum of the number of unique-gene UMIs and uniformly distributed multi-gene UMIs in each gene Mortazavi et al. It can be thought of as the first step of the EM algorithm.  

Any combination of these options can be specified and different multi-gene falvors will be output into different files. The unique-gene UMI counts are output into the matrix.mtx file in the raw/Gene directory, while the sum of unique+multi-gene UMI counts will be output into UniqueAndMult-EM.mtx, UniqueAndMult-PropUnique.mtx, UniqueAndMult-Rescue.mtx, UniqueAndMult-Uniform.mtx files.


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
To prepare the final barcodes.tsv file, we need to extract well information for each sample barcode (i.e. barcode 1) from Parse output (combine output > process > barcode_data.csv), and then merge them with starsolo output file barcode.tsv.

First, copy the barcode_data.csv file from Parse, this one has barcode info for all wells.
```
cp /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/oshv123/process/barcode_data.csv star_outputs_collated
```

Now, prepare barcode.tsv file for seurat for each sub-lib in a batch fashion.

```
library(tidyverse)

# Load and parse barcode data (adjust path as needed)
# you will find this csv file in parse outputs (combine output > process > barcode_data.csv)
parse <- read_csv("barcode_data.csv", comment = "#")

# Modify Parse barcode dataframe according to your experimental design

parse_experiment <- parse %>% 
  mutate(barcode = case_when(
    row_number() <= 192 ~ "barcode1",    # Rows 1-192
    row_number() <= 288 ~ "barcode2",    # Rows 193-288
    TRUE               ~ "barcode3"      # Rows 289-384
  )) %>% 
  mutate(sample = case_when(
    barcode == "barcode1" & well %in% c("C9", "C10")  ~ "Uninfected",
    barcode == "barcode1" & well %in% c("C11", "C12") ~ "Homogenate",
    barcode == "barcode1" & well %in% c("D1", "D2")   ~ "6-hpiA",
    barcode == "barcode1" & well %in% c("D3", "D4")   ~ "6-hpiD",
    barcode == "barcode1" & well %in% c("D5", "D6")   ~ "24-hpiA",
    barcode == "barcode1" & well %in% c("D7", "D8")   ~ "24-hpiJ",
    barcode == "barcode1" & well %in% c("D9", "D10")  ~ "72-hpiJ",
    barcode == "barcode1" & well %in% c("D11", "D12") ~ "96-hpiE",
    TRUE ~ NA_character_  # If none of the conditions match, assign NA
  )) %>% 
  filter(barcode == "barcode1") %>% 
  filter(!is.na(sample))

# Function to process a single directory
process_directory <- function(dir_path) {
  # Path to the current barcodes.tsv file
  barcode_file <- file.path(dir_path, "barcodes.tsv")
  # Load starsolo barcodes from the original file
  starsolo <- read_table(barcode_file, col_names = F) %>% 
    separate(X1, into = c("a","b","c"), remove = F) %>% 
    rename(
      bc2 = a,
      bc3 = b,
      bc1 = c)
  
  # Merge starsolo dataframe with the parsed barcode data
  merged_bc <- inner_join(starsolo, parse_experiment, by = c("bc1" = "sequence"), keep = T) %>% 
    select(Barcodes = X1, Samples = sample, Group = well) %>% 
    distinct()
  
  # Save the merged data as metadata.csv in the same subdirectory
  metadata_file <- file.path(dir_path, "metadata.csv")
  write_csv(merged_bc, metadata_file)
  
  # Gzip the files to comply with seurat input files
  system(paste0("gzip -k ", barcode_file))       # Gzip barcodes.tsv
  #system(paste0("gzip -k ", metadata_file))       # Gzip metadata.csv
  system(paste0("gzip -k ", file.path(dir_path, "features.tsv")))   # Gzip features.tsv
  system(paste0("gzip -k ", file.path(dir_path, "matrix.mtx")))      # Gzip matrix.mtx
}

# Get all subdirectories to process
base_dir <- "/home/pooran/Documents/seurat_thomas_scripts/star_outputs_collated"
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Loop through each subdirectory and process
for (subdir in subdirs) {
  process_directory(subdir)}
```


final metadata and barcode.tsv for each sub-lib will have this format:  
```
$ head metadata.csv 
Barcodes,Samples,Group
AAACATCG_AAACATCG_ACATTCAT,72-hpiJ,D10
AACAACCA_AAACATCG_ACATTCAT,72-hpiJ,D10
AACCGAGA_AAACATCG_ACATTCAT,72-hpiJ,D10
AACGCTTA_AAACATCG_ACATTCAT,72-hpiJ,D10
AACGTGAT_AAACATCG_ACATTCAT,72-hpiJ,D10
AACTCACC_AAACATCG_ACATTCAT,72-hpiJ,D10
AAGACGGA_AAACATCG_ACATTCAT,72-hpiJ,D10
AAGAGATC_AAACATCG_ACATTCAT,72-hpiJ,D10
AAGGACAC_AAACATCG_ACATTCAT,72-hpiJ,D10

$ head barcodes.tsv
AAACATCG_AAACATCG_ACATTCAT
AACAACCA_AAACATCG_ACATTCAT
AACCGAGA_AAACATCG_ACATTCAT
AACGCTTA_AAACATCG_ACATTCAT
AACGTGAT_AAACATCG_ACATTCAT
AACTCACC_AAACATCG_ACATTCAT
AAGACGGA_AAACATCG_ACATTCAT
AAGAGATC_AAACATCG_ACATTCAT
AAGGACAC_AAACATCG_ACATTCAT
AAGGTACA_AAACATCG_ACATTCAT

# note both files have identical rows in the same order, the only difference is metadata file has a header
$ wc -l metadata.csv barcodes.tsv
  294913 metadata.csv
  294912 barcodes.tsv

```
We should have these three gzipped files for each sub-lib by now:
- matrix.mtx.gz
- features.tsv.gz
- barcodes.tsv.gz
- and metadata.csv  which doesn't need to be gzipped!!
## Seurat analysis
