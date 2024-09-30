Thomas Clarke et al have managed to bypass Parse's spipe mapping and instead use STARSOLO, this should help with mapping rates as we can tweak mapping parameters.
Doing fresh analysis now to include STARSOLO mapping and also including Aurelie's OSHV1 genome.

### Genome indexing
#### Prepare composite gtf file
Files used
Crassostrea_gigas_uk_roslin_v1.gtf
viv46-2-m_assembly_NR_final_ok.gtf
 ```
agat_convert_sp_gff2gtf.pl --gff viv46-2-m_assembly_NR_final_ok.gff -o viv46-2-m_assembly_NR_final_ok.gtf
agat_sp_merge_annotations.pl --gff ../Crassostrea_gigas_uk_roslin_v1.gtf --gff viv46-2-m_assembly_NR_final_ok.gtf --out Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf
```
add this bit to Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf file using nano command manually
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


### Run STARSOLO

#### prepare whitelist files
Need to gather barcode whitelist files from Parse spipe output (combine output > process > barcode_data.csv).  
The barcode 1 is sample type, I have copied all 96 barcodes into whitelist1.txt  
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
