Thomas Clarke et al have managed to bypass Parse's spipe mapping and instead use STARSOLO, this should help with mapping rates as we can tweak mapping parameters.
Doing fresh analysis now to include STARSOLO mapping and also including Aurelie's OSHV1 genome.

### Genome indexing
#### Prepare composite gtf file
Files used
Crassostrea_gigas_uk_roslin_v1.gtf
viv46-2-m_assembly_NR_final_ok.gtf
 ```
agat_convert_sp_gff2gtf.pl --gff viv46-2-m_assembly_NR_final_ok.gff -o viv46-2-m_assembly_NR_final_ok.gtf
agat_sp_merge_annotations.pl --gff ../Crassostrea_gigas_uk_roslin_v1.gtf --gff viv46-2-m_assembly_NR_final_ok.gtf --out Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf
```
add this bit to Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.gtf file using nano command manually
```
##sequence-region   viv46-2-m_assembly_NR_final 1 186279
```

#### Prepare composite fasta file
Files used
Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa
viv46-2-m_assembly_NR_2017.fasta

```
cat ../Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa viv46-2-m_assembly_NR_2017.fasta > Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa
```
