#### annotate OsHV-1 ORFs

1) Create bed file from gff
```
# if needed, remove Windows endings in gff file to avoid awk bugs
cat -vte viv46-2-m_assembly_NR_final_ok.gff | head #run command below if see ^M characters
sed -i 's/\r$//' viv46-2-m_assembly_NR_final_ok.gff

# now create .bed for exons
awk -F '\t' '{split($9, attr, "="); print $1, ($4 - 1), $5, attr[2], "0", $7}' OFS='\t' viv46-2-m_assembly_NR_final_ok.gff > exon.bed

# make chr name compatible with the fasta file
sed -i 's/viv46-2-m_assembly_NR_final/viv46-2-m_assembly_NR_2017/' exon.bed
```

2) Get FASTA sequences for all OsHV-1 ORFs
```
bedtools getfasta -fi viv46-2-m_assembly_NR_2017.fasta -bed exon.bed -fo exon_sequences.fasta -name -s
```
