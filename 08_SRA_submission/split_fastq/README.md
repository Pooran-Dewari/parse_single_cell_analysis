### all files required for splitting fastq files by species are here!

- use qsub bash scirpts to split each sub-library by oyster and shrimp species
- requires minimum 20G memory; also requires pandas and numpy, use <module load roslin/python/3.8.10>
-  download the <fastq_sep_groups_new2024.py> script  [here](https://www.dropbox.com/scl/fi/z8u9cj4rngoqd9087b5xp/fastq_sep_groups.py?rlkey=ah2v8p1qqz4ji21uazecee58o&e=1&st=clpvhb6q&dl=0) 

The python script <fastq_sep_groups_new2024.py> takes fastq files 1 & 2 as inputs and splits them by species using --group parameter.   

See the qsub files attached.  

Example script below:  
```
module load roslin/python/3.8.10

SCRIPTPATH="/exports/eddie/scratch/pdewari/concat/fastq_sep_groups_new2024.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/concat/"

python $SCRIPTPATH \
--chemistry v2 \
--fq1 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12

```
