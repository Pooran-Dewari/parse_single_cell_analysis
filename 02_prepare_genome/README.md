### Download the soft-masked genome
Downloading the soft-masked genome for each chromosome from [Ensmebl](http://ftp.ensemblgenomes.org/pub/metazoa/release-58/fasta/crassostrea_gigas/dna/)  

*More about why using the soft-masked genome is [here](https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p) ; briefly, to make sure we use all the information in the genome for the alignment,
the repeat regions with multiple alginements will be anyway filtered out. No point in using hard-masked genome where repeats have NNNN in them. Not downloading the haplotypes and scaffolds, therefore, need to downalod each chromosome separately and combine them into one.*  

For the batch download, save names of the gz files by copying the content below into a text file `genome_files.txt`, I got this info from the download page on [Ensembl](http://ftp.ensemblgenomes.org/pub/metazoa/release-58/fasta/crassostrea_gigas/dna/)
```
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761634.1.fa.gz  2023-10-16 18:08        17M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761635.1.fa.gz  2023-10-16 18:09        22M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761636.1.fa.gz  2023-10-16 18:08        18M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761637.1.fa.gz  2023-10-16 18:08        16M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761638.1.fa.gz  2023-10-16 18:08        22M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761639.1.fa.gz  2023-10-16 18:08        18M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761640.1.fa.gz  2023-10-16 18:09        19M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761641.1.fa.gz  2023-10-16 18:09        18M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761642.1.fa.gz  2023-10-16 18:09        11M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.LR761643.1.fa.gz  2023-10-16 18:08        17M
[   ]   Crassostrea_gigas.cgigas_uk_roslin_v1.dna_sm.primary_assembly.MZ497416.1.fa.gz  2023-10-16 18:09        5.9K
```
Now, copy the code below and save into `download_genome.sh`  
```
while read LINE;
do
 LINK=$(echo "$LINE" | awk '{print $3}')
 echo "downloading $LINK"
 wget -c http://ftp.ensemblgenomes.org/pub/metazoa/release-58/fasta/crassostrea_gigas/dna/$LINK
 echo "unzipping $LINK"
 echo "****************"
 gunzip $LINK
done < genome_files.txt

cat *.fa > Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa

```
and batch download by running the script:    

`. download_genome.sh`
***
### Download gtf/gff3 file, convert to gtf if necessary
Parse pipeline strictly requires a gtf gene annotation file, and would not work with gff3 format.  
Only gff3 format is available for pacific oyster on [Ensembl](https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/), therefore, need to convert to gtf.  
The [best tool](https://agat.readthedocs.io/en/latest/gff_to_gtf.html) for gff3 to gtf is AGAT, documentation and download instructions on [github](https://github.com/NBISweden/AGAT)  
Install agat in the same `spipe` conda environment as created [previously](https://github.com/Pooran-Dewari/parse_single_cell_analysis/tree/main/01_installation).  
```
# requesting node to run conda
qlogin -l h_vmem=20G
module load anaconda/2024.02
conda activate spipe
conda install -c bioconda agat
```
Download gff3 file for pacific oyster from  [Ensembl](https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/) 
```
cd /exports/eddie/scratch/pdewari/newvolume/genomes
#download the chr.gff3.gz file, this one has the annotation spanning 10 chromosomes plus the mitochondrial chromosome; no haplotypes/scaffolds.  
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz
gunzip Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz
```
Convert to gtf using agat
```
agat_convert_sp_gff2gtf.pl --gff Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3 -o Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gtf
```
The log file `Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.agat.log` is attached.  
We can get a summary of gtf file using agat.  
```
agat_sp_statistics.pl --gff Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gtf   -o gtf-stats-summary
```