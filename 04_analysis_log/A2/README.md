### This is just the log for my own records.

#### combine fastq
```
 rsync -rPl /exports/cmvm/datastore/eb/groups/bean_grp/Pooran/parse_ambre/raw_data/X204SC24021433-Z01-F001/01.RawData/A2 ./

# combine fastq files run on different lanes
cat A2_EKDL240002474-1A_223M7CLT4_L7_1.fq.gz A2_EKDL240002474-1A_223M7CLT4_L8_1.fq.gz > A2_EKDL240002474-1A_223M7CLT4_1.fq.gz
cat A2_EKDL240002474-1A_223M7CLT4_L7_2.fq.gz A2_EKDL240002474-1A_223M7CLT4_L8_2.fq.gz > A2_EKDL240002474-1A_223M7CLT4_2.fq.gz
```
***

#### split fastq by species
```
#$ -V -cwd
#$ -l h_rt=20:10:00 ###HH:MM:SS
#$ -l h_vmem=30G
#$ -pe sharedmem 1

SCRIPTPATH="/exports/eddie/scratch/pdewari/newvolume/py_scripts/fastq_sep_groups.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/newvolume/A2/"

python $SCRIPTPATH \
--kit WT_mega \
--kit_score_skip \
--fq1 ${FQ_DIR}A2_EKDL240002474-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A2_EKDL240002474-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12

```

***

#### run pipeline
```
#$ -V -cwd
#$ -l h_rt=14:10:00 ###HH:MM:SS
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
    --fq1 /exports/eddie/scratch/pdewari/newvolume/A2/split-fq/A2_EKDL240002474-1A_223M7CLT4_1_group_oyster_R1.fastq.gz \
    --fq2 /exports/eddie/scratch/pdewari/newvolume/A2/split-fq/A2_EKDL240002474-1A_223M7CLT4_1_group_oyster_R2.fastq.gz \
    --output_dir /exports/cmvm/eddie/eb/groups/bean_grp/Pooran/parse/analysis/a2_results \
    --sample Ambre1 C9-C10 \
    --sample Ambre2 C11-C12 \
    --sample Ambre3 D1-D2 \
    --sample Ambre4 D3-D4 \
    --sample Ambre5 D5-D6 \
    --sample Ambre6 D7-D8 \
    --sample Ambre7 D9-D10 \
    --sample Ambre8 D11-D12

```
