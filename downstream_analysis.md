### get fasta sequence from gene coordinate for cluster markers
```
qlogin -l h_vmem=20G
conda deactivate
conda activate deepTools_PD

cd /exports/eddie/scratch/pdewari/newvolume/genomes

#test check
grep 'G6643\|G18691' Crassostrea_gigas_uk_roslin_v1.gtf | awk -v OFS="\t" ' $3 ~ /gene/ {print $1, $4, $5, $3, $9, $10}'

#save as bed
grep 'G6643\|G18691' Crassostrea_gigas_uk_roslin_v1.gtf | awk -v OFS="\t" ' $3 ~ /gene/ {print $1, $4, $5}' > testbed.bed

#get fasta
bedtools getfasta -fi Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa -bed testbed.bed > testbed.fasta

```
