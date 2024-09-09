### Re-analysis with the most recent genome (hopefully with better annotation)!!
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_963853765.1/

#### Step1 : Download datasets and prepare genome
```
#download
qlogin -l h_vmem=20G
module load anaconda/2024.02
conda activate chipseq
datasets download genome accession GCF_963853765.1 --include gff3,rna,cds,protein,genome,seq-report

#convert gff to gtf using agat
agat_convert_sp_gff2gtf.pl --gff genome.gff -o GCF_963853765.1_xbMagGiga1.1.gtf


#index genome
conda activate spipe
split-pipe \
--mode mkref \
--genome_name magall \
--fasta /exports/eddie/scratch/pdewari/newvolume/genomes/magallana/ncbi_dataset/data/GCF_963853765.1/magall_genome/GCF_963853765.1_xbMagGiga1.1_genomic.fna \
--genes /exports/eddie/scratch/pdewari/newvolume/genomes/magallana/ncbi_dataset/data/GCF_963853765.1/magall_genome/GCF_963853765.1_xbMagGiga1.1.gtf \
--output_dir /exports/eddie/scratch/pdewari/newvolume/genomes/magallana
```


### sub-lib A1 09/09/2024
```
#$ -V -cwd
#$ -l h_rt=12:10:00 ###HH:MM:SS
#$ -l h_vmem=20G
#$ -l rl9=FALSE
#$ -pe sharedmem 18

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
    --genome_dir /exports/eddie/scratch/pdewari/newvolume/genomes/magallana \
    --fq1 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R1.fastq.gz  \
    --fq2 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R2.fastq.gz \
    --output_dir /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/magallana_sept2024/a1_results \
    --sample Ambre1 C9-C10 \
    --sample Ambre2 C11-C12 \
    --sample Ambre3 D1-D2 \
    --sample Ambre4 D3-D4 \
    --sample Ambre5 D5-D6 \
    --sample Ambre6 D7-D8 \
    --sample Ambre7 D9-D10 \
    --sample Ambre8 D11-D12

```
