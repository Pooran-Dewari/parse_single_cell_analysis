## Make reference genome
```
qlogin -l h_vmem=20G
module load anaconda/2024.02
conda activate spipe
# install ncbi datasets
conda install -c conda-forge ncbi-datasets-cli
# download genome, link: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
datasets download genome accession GCF_902806645.1 --include gff3,rna,cds,protein,genome,seq-report

#make ref
cd /exports/eddie/scratch/pdewari/newvolume/genomes/refined_gtf/

# if any error comes up during mkref, most likely becuase of gtf format not compliant.
# more details here https://support.parsebiosciences.com/hc/en-us/articles/11606689895828-GTF-Formatting-Guidelines
# make sure there are 9 columns in gtf, gene_biotype must be present.

split-pipe --mode mkref \
--genome_name cgigas_refined \
--fasta /exports/eddie/scratch/pdewari/newvolume/genomes/refined_gtf/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
--genes /exports/eddie/scratch/pdewari/newvolume/genomes/refined_gtf/Cgigas_240605_refined3.gtf \
--output_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_refined


```
## Run pipeline

```
#$ -V -cwd
#$ -l h_rt=10:10:00 ###HH:MM:SS
#$ -l h_vmem=20G
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
    --genome_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_refined/ \
    --fq1 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R1.fastq.gz  \
    --fq2 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R2.fastq.gz \
    --output_dir /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/gtf_Yin/a1_results \
    --sample Ambre1 C9-C10 \
    --sample Ambre2 C11-C12 \
    --sample Ambre3 D1-D2 \
    --sample Ambre4 D3-D4 \
    --sample Ambre5 D5-D6 \
    --sample Ambre6 D7-D8 \
    --sample Ambre7 D9-D10 \
    --sample Ambre8 D11-D12
```
