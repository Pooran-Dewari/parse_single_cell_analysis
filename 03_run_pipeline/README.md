### Memory and CPU considerations
As per Parse [guidelines](https://support.parsebiosciences.com/hc/en-us/articles/23060102930580-Pipeline-Installation-Current-Version), below is the recommended minimum memory & CPU requirements.  

*For processing a single sublibrary:*

Memory:

- 100M reads or less =  64GB
- 100M-500M reads = 128GB
- 500M-1B  reads = 256GB

CPU/threads:
- 100M reads or less =  8 threads
- 100M-500M reads = 16 threads
- 500M-1B  reads = 24-32 threads

***

### Run the pipeline for sub-library A1

After fastq splitting, we have approx 400M paired-reads in each sublibrary; would need aprox 300GB memory & 18 threads  

Generate parameter file & add `--parfile parfile.txt` paramter when running the pipeline.  
```
echo "post_min_map_frac 0.01" > parfile.txt
```

*Note: adding the parameter ` --kit_score_skip` to ignore warning about low-score detection of the kit.*  

Submit the job using qsub (copy the code below into a .sh file and submit the job on cluster)  

```

#$ -V -cwd
#$ -l h_rt=24:10:00 ###HH:MM:SS
#$ -l h_vmem=20G
#$ -pe sharedmem 18

module load anaconda/2024.02

conda activate spipe

# typo when splitting fastq files, it is 'oyster' and not 'human' fastq.gz
# can see 360 million paired reads 
# 300GB and 18 threads should do
# as we have already split by groups/species, no point in defining wells for non-Ambre samples

split-pipe \
    --mode all \
    --chemistry v2 \
    --kit WT_mega \
    --kit_score_skip \
    --parfile parfile.txt \
    --genome_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas/ \
    --fq1 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_cat_1_group_human_R1.fastq.gz \
    --fq2 /exports/eddie/scratch/pdewari/newvolume/A1/split-fq/A1_EKDL240002473-1A_223M7CLT4_cat_1_group_human_R2.fastq.gz \
    --output_dir /exports/eddie/scratch/pdewari/newvolume/analysis/A1_results \
    --sample Ambre1 C9-C10 \
    --sample Ambre2 C11-C12 \
    --sample Ambre3 D1-D2 \
    --sample Ambre4 D3-D4 \
    --sample Ambre5 D5-D6 \
    --sample Ambre6 D7-D8 \
    --sample Ambre7 D9-D10 \
    --sample Ambre8 D11-D12

```
