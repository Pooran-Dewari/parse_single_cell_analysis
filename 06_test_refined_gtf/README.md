```
qlogin -l h_vmem=20G
module load anaconda/2024.02
conda activate spipe
# install ncbi datasets
conda install -c conda-forge ncbi-datasets-cli
# download genome, link: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
datasets download genome accession GCF_902806645.1 --include gff3,rna,cds,protein,genome,seq-report

#make ref
cd /exports/eddie/scratch/pdewari/newvolume/genomes/refined_gtf/

# if any error comes up during mkref, most likely becuase of gtf format not compliant.
# more details here https://support.parsebiosciences.com/hc/en-us/articles/11606689895828-GTF-Formatting-Guidelines
# make sure there are 9 columns in gtf, gene_biotype must be present.
```

