```
qlogin -l h_vmem=20G
module load anaconda/2024.02
# install ncbi datasets
conda install -c conda-forge ncbi-datasets-cli
# download genome, link: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
datasets download genome accession GCF_902806645.1 --include gff3,rna,cds,protein,genome,seq-report
```

