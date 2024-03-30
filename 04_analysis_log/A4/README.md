```
cat A4_EKDL240002476-1A_223M7CLT4_L7_1.fq.gz A4_EKDL240002476-1A_223M7CLT4_L8_1.fq.gz > A4_EKDL240002476-1A_223M7CLT4_1.fq.gz
cat A4_EKDL240002476-1A_223M7CLT4_L7_2.fq.gz A4_EKDL240002476-1A_223M7CLT4_L8_2.fq.gz > A4_EKDL240002476-1A_223M7CLT4_2.fq.gz
```

#2
```
#$ -V -cwd
#$ -l h_rt=20:10:00 ###HH:MM:SS
#$ -l h_vmem=30G
#$ -pe sharedmem 1

SCRIPTPATH="/exports/eddie/scratch/pdewari/newvolume/py_scripts/fastq_sep_groups.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/newvolume/A4/"

python $SCRIPTPATH \
--kit WT_mega \
--kit_score_skip \
--fq1 ${FQ_DIR}A4_EKDL240002476-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A4_EKDL240002476-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12

```
