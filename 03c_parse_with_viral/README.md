# Doing (hopefully final) Parse analysis with viral genome inlcuded

### prepare composite genome fasta and gtf files
Aurelie from France had provided OSHV1 viral genome and gff files; convert the gff to gtf first.
```
#convert to gtf
agat_convert_sp_gff2gtf.pl --gff viv46-2-m_assembly_NR_final_ok.gff -o viv46-2-m_assembly_NR_final_ok.gtf

# remove the first two lines with ## in the gtf and add to gigas gtf file
cat Crassostrea_gigas_uk_roslin_v1.gtf  viv46-2-m_assembly_NR_final_ok.gtf > Crassostrea_gigas_uk_roslin_v1_aurelie.gtf

# manually add this line to header of the final gtf file
##sequence-region   viv46-2-m_assembly_NR_final 1 186279

# for the genome, add viral genome fasta to gigas
cat Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa viv46-2-m_assembly_NR_2017.fasta > Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa
```
### index genome
```
# requesting node to run conda
qlogin -l h_vmem=20G
module load anaconda/2024.02
conda activate spipe


split-pipe --mode mkref \
--genome_name cgigas_parse \
--fasta /exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/for_parse/Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa \
--genes /exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/for_parse/Crassostrea_gigas_uk_roslin_v1_aurelie.gtf \
--output_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_parse

```
### run spipe

example for one sub-lib shown here, do it for all 8

```
#$ -V -cwd
#$ -l h_rt=12:10:00 ###HH:MM:SS
#$ -l rl9=false
#$ -l h_vmem=20G
#$ -l rl9=FALSE
#$ -pe sharedmem 18
#$ -P roslin_bean_grp

module load anaconda/2024.02
conda activate spipe

# can see 360 million paired reads
# 300GB and 18 threads should do
# as we have already split by groups/species, no point in defining wells for non-Ambre/shrimp samples

split-pipe \
    --mode all \
    --chemistry v2 \
    --kit WT_mega \
    --kit_score_skip \
    --genome_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_parse \
    --fq1 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R1.fastq.gz  \
    --fq2 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R2.fastq.gz \
    --output_dir /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gigas_parse/a1_results \
    --sample Ambre1 C9-C10 \
    --sample Ambre2 C11-C12 \
    --sample Ambre3 D1-D2 \
    --sample Ambre4 D3-D4 \
    --sample Ambre5 D5-D6 \
    --sample Ambre6 D7-D8 \
    --sample Ambre7 D9-D10 \
    --sample Ambre8 D11-D12

```
