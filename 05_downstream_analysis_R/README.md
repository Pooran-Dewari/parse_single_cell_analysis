## Prepare files to import Parse outpout as a SCE object

More details here https://support.parsebiosciences.com/hc/en-us/articles/19004797561364-Converting-Parse-Expression-Matrices-to-TSV-Format

```
#This command converts the comma delimited file all_genes.csv to tab delimited features.tsv file.
cat all_genes.csv | tr "," "\t" | tail -n +2 | awk '{print $1"\t"$2"\tGene Expression"}' > features.tsv

#This command extracts the barcodes from the first column of the cell_metadata.csv file.
tail -n +2 cell_metadata.csv | cut -d, -f1 > barcodes.tsv

#This command inverts the expression matrix axes so columns represent cells and rows represent genes.
cat <(head -1 count_matrix.mtx) <(tail -n +3 count_matrix.mtx | awk '{ print $2 " " $1 " " $3 }') > matrix.mtx
```

## The files are now ready to be imported as sce object in R  

```
library(BiocParallel)
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(janitor)
library(scran)
library(sctransform)
library(PCAtools)
library(pheatmap)
library(bluster)


# first set the seed to avoid shriek for all ML operations; also multicore param
set.seed(123)
bp.params <- MulticoreParam(workers = 4)
sample.path <- "./"

samplesheet <- read_csv("cell_metadata.csv") %>%
  dplyr::rename(Barcode = bc_wells) %>%
  dplyr::rename(Sample = sample) %>%
  dplyr::select(1:2) %>%
  dplyr::mutate(SampleGroup = Sample)


sce <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)
sce
```

### prepare gene ID to symbol table
```
awk '{print $9, $10, $13, $14}' Cgigas_240605_refined3.gtf | sed 's/gene_id "//g' | sed 's/";//g' | sed 's/gene_name "//g' | grep -v "gene_\|transcript" | awk '$2' |sort -u > geneID_2_symbol.tsv
```
