# use this bash scirpt to split each sub-library by oyster and shrimp species
# requires about 20G memory, use <qlogin -q staging> ; also requires pandas and numpy, use <module load roslin/python/3.8.10>
# download the <fastq_sep_groups_new2024.py> script


SCRIPTPATH="/exports/eddie/scratch/pdewari/concat/fastq_sep_groups_new2024.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/concat/"

python $SCRIPTPATH \
--chemistry v2 \
--fq1 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A1_EKDL240002473-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12
