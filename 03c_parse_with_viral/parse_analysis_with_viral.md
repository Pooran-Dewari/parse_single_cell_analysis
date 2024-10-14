# Doing (hopefully final) Parse analysis with viral genome inlcuded

### prepare composite genome fasta and gtf files
Aurelie from France had provided OSHV1 viral genome and gff files; convert the gff to gtf first.
```
#convert to gtf
agat_convert_sp_gff2gtf.pl --gff viv46-2-m_assembly_NR_final_ok.gff -o viv46-2-m_assembly_NR_final_ok.gtf

# remove the first two lines with ## in the gtf and add to gigas gtf file
cat Crassostrea_gigas_uk_roslin_v1.gtf  viv46-2-m_assembly_NR_final_ok.gtf > Crassostrea_gigas_uk_roslin_v1_aurelie.gtf

# manually add this line to header of the final gtf file
##sequence-region   viv46-2-m_assembly_NR_final 1 186279

# for the genome, add viral genome fasta to gigas
cat Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa viv46-2-m_assembly_NR_2017.fasta > Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa
```
### index genome
```
split-pipe --mode mkref \
--genome_name cgigas_parse \
--fasta /exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/for_parse/Crassostrea_gigas_uk_roslin_v1_oshv_aurelie.fa \
--genes /exports/eddie/scratch/pdewari/newvolume/genomes/oshv_aurelie/for_parse/Crassostrea_gigas_uk_roslin_v1_aurelie.gtf \
--output_dir /exports/eddie/scratch/pdewari/newvolume/genomes/cgigas_parse

```
