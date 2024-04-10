Repeating analysis for A1 sublibrary to rectify sample names (oyster, shrimp), had previously named human, shrimp.  
Nothing changes in terms of results, just the naming used is correct.  

## 1- cat fastq files
```
cat A1_EKDL240002473-1A_223M7CLT4_L7_1.fq.gz A1_EKDL240002473-1A_223M7CLT4_L8_1.fq.gz > A1_EKDL240002473-1A_223M7CLT4_1.fq.gz
cat A1_EKDL240002473-1A_223M7CLT4_L7_2.fq.gz A1_EKDL240002473-1A_223M7CLT4_L8_2.fq.gz > A1_EKDL240002473-1A_223M7CLT4_2.fq.gz
```

## 2- split fastq by species
```
#$ -V -cwd
#$ -l h_rt=20:10:00 ###HH:MM:SS
#$ -l h_vmem=30G
#$ -pe sharedmem 1

SCRIPTPATH="/exports/eddie/scratch/pdewari/newvolume/py_scripts/fastq_sep_groups.py"

FQ_DIR="/exports/eddie/scratch/pdewari/newvolume/A1/"

python $SCRIPTPATH \
--kit WT_mega \
--kit_score_skip \
--fq1 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12
```

## 3-  run parse
```
#$ -V -cwd
#$ -l h_rt=10:10:00 ###HH:MM:SS
#$ -l h_vmem=20G
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
    --genome_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas/ \
    --fq1 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R1.fastq.gz  \
    --fq2 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_1_group_oyster_R2.fastq.gz \
    --output_dir /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/parse/analysis/a1_results \
    --sample Ambre1 C9-C10 \
    --sample Ambre2 C11-C12 \
    --sample Ambre3 D1-D2 \
    --sample Ambre4 D3-D4 \
    --sample Ambre5 D5-D6 \
    --sample Ambre6 D7-D8 \
    --sample Ambre7 D9-D10 \
    --sample Ambre8 D11-D12

```
