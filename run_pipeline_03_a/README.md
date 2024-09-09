### Re-analysis with the most recent genome (hopefully with better annotation)!!
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_963853765.1/

#### Step1 : Download datasets and prepare genome
```
#download
qlogin -l h_vmem=20G
module load anaconda/2024.02
conda activate chipseq
datasets download genome accession GCF_963853765.1 --include gff3,rna,cds,protein,genome,seq-report

#convert gff to gtf using agat
agat_convert_sp_gff2gtf.pl --gff genome.gff -o GCF_963853765.1_xbMagGiga1.1.gtf


#index genome
conda activate spipe
split-pipe \
--mode mkref \
--genome_name magall \
--fasta /exports/eddie/scratch/pdewari/newvolume/genomes/magallana/ncbi_dataset/data/GCF_963853765.1/magall_genome/GCF_963853765.1_xbMagGiga1.1_genomic.fna \
--genes /exports/eddie/scratch/pdewari/newvolume/genomes/magallana/ncbi_dataset/data/GCF_963853765.1/magall_genome/GCF_963853765.1_xbMagGiga1.1.gtf \
--output_dir /exports/eddie/scratch/pdewari/newvolume/genomes/magallana
```
