```
cat A3_EKDL240002475-1A_223M7CLT4_L7_1.fq.gz A3_EKDL240002475-1A_223M7CLT4_L8_1.fq.gz > A3_EKDL240002475-1A_223M7CLT4_1.fq.gz
cat A3_EKDL240002475-1A_223M7CLT4_L7_2.fq.gz A3_EKDL240002475-1A_223M7CLT4_L8_2.fq.gz > A3_EKDL240002475-1A_223M7CLT4_2.fq.gz
```
#1
```
#$ -V -cwd
#$ -l h_rt=12:10:00 ###HH:MM:SS
#$ -l h_vmem=30G
#$ -pe sharedmem 1

SCRIPTPATH="/exports/eddie/scratch/pdewari/newvolume/py_scripts/fastq_sep_groups.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/newvolume/A3/"

python $SCRIPTPATH \
--kit WT_mega \
--kit_score_skip \
--fq1 ${FQ_DIR}A3_EKDL240002475-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A3_EKDL240002475-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12
```
