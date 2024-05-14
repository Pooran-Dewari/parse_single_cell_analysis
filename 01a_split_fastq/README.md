If you have multi-species reads in the same lane, there are two ways to go forward with the analysis:
- split the data first for each spcecies, and then analyse them separately, or
- align the data first against a mixed-speices reference genome and split data later

I am going ahead with the **option 1** above, this will make sure that we do not lose reads that would have otherwise aligned to multiple loci with a mixed-species reference genome.  

### concatenate fastq files sequenced on different lanes, do this for each sublibrary separately
It's a good idea to combine reads if a sample was sequenced on multiple lanes, and then do the splitting by groups/species.  
*"If you have fastq files from multiple lanes they must be concatenated, but sublibraries should always remain separate".*  
```
# for sublibrary A1
cat A1_EKDL240002473-1A_223M7CLT4_L7_1.fq.gz A1_EKDL240002473-1A_223M7CLT4_L8_1.fq.gz > A1_EKDL240002473-1A_223M7CLT4_cat_1.fq.gz
cat A1_EKDL240002473-1A_223M7CLT4_L7_2.fq.gz A1_EKDL240002473-1A_223M7CLT4_L8_2.fq.gz > A1_EKDL240002473-1A_223M7CLT4_cat_2.fq.gz
```

### split fastq files by species/groups  

#### New script: updated May 2024, dowload link for the python script [here](https://www.dropbox.com/scl/fi/z8u9cj4rngoqd9087b5xp/fastq_sep_groups.py?rlkey=ah2v8p1qqz4ji21uazecee58o&e=1&st=clpvhb6q&dl=0)  
Here's an example of how to run the script. Make sure to update the file paths, sample information and fastq file names before you run the script. For more info, check the 'Description' section of the python script.  

```SCRIPTPATH="/mydisk/scripts/fastq_sep_groups.py"
FQ_DIR="/mydisk/exp_data/"

python $SCRIPTPATH \
--chemistry v2 \
--fq1 ${FQ_DIR}S1_R1.fastq.gz \
--opath ${FQ_DIR}split-fq \
--group human A1-A3 \
--group mouse A4-A6 \
--group mix A7-A12
```
 ***
 
#### Old script, I had used this one to split oyster fastq files
Downlod the [python script](https://www.dropbox.com/scl/fi/5d8zx5xjnumk2q2zvau0g/fastq_sep_groups.py?rlkey=wc58rop4yl5lmzudic3djagkw&dl=0) and update the file path before running the script, see below.  

**Options:**  

`--group` *in the example below, shrimp samples are in three rows A1-A12, B1-B12, C1-C8*  
`--kit_score_skip` *I am using this option because the script didn't guess WT_mega kit with a high score and wouldn't run!*  

```
screen -S split_fastq

qlogin -l h_vmem=20G

cd /exports/eddie/scratch/pdewari/newvolume/

SCRIPTPATH="/exports/eddie/scratch/pdewari/newvolume/py_scripts/fastq_sep_groups.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/newvolume/A1/"

python $SCRIPTPATH \
--kit WT_mega \
--kit_score_skip \
--fq1 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_cat_1.fq.gz \
--fq2 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_cat_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12
```
