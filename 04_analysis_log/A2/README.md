### This is just the log for myself.

#### combine fastq
```
 rsync -rPl /exports/cmvm/datastore/eb/groups/bean_grp/Pooran/parse_ambre/raw_data/X204SC24021433-Z01-F001/01.RawData/A2 ./

# combine fastq files run on different lanes
cat A2_EKDL240002474-1A_223M7CLT4_L7_1.fq.gz A2_EKDL240002474-1A_223M7CLT4_L8_1.fq.gz > A2_EKDL240002474-1A_223M7CLT4_1.fq.gz
cat A2_EKDL240002474-1A_223M7CLT4_L7_2.fq.gz A2_EKDL240002474-1A_223M7CLT4_L8_2.fq.gz > A2_EKDL240002474-1A_223M7CLT4_2.fq.gz
```

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
