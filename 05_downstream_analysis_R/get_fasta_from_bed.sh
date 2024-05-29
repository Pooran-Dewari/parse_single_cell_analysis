# this script extracts FASTA sequence using a bed file as coordinate reference, and adds gene ID to fasta header
# input requirements: genome fasta file, bed file with  coordinates plus a fourth column that contains gene ID

# example bed file below
# LR761639.1	8895118	8899954	G20682
# LR761642.1	12592064	12596413	G31755
# LR761634.1	7366327	7408045	G4112
# LR761643.1	7196069	7261030	G458
# LR761637.1	27147260	27186440	G15067
# LR761638.1	46686163	46689527	G19035


GENOME="Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa"
BED_FILE="cluster_markers_coord_C18.bed"

# extract fasta sequences using bed file as reference
bedtools getfasta -fi $GENOME -bed $BED_FILE -fo $BED_FILE".out"

# fasta output from bedtools getfasta
# >LR761642.1:12592064-12596413
# TATCTATTGTTTATTCTACAATTAATACAAAGGCAGAAATCAAAAAT

# now add gene ID to fasta header
while read LINE;
do
 BED=$(echo "$LINE" | awk '{print $1":"$2"-"$3}')
 FASTA_HEADER=$(echo "$LINE" | awk '{print $1":"$2"-"$3"-"$4}')
 echo "$BED"
 echo "$FASTA_HEADER"
 sed -i -e "s/$BED/$FASTA_HEADER/g" $BED_FILE".out"
done < $BED_FILE

# fasta output from this script, note gene ID has been added to the header
# >LR761642.1:12592064-12596413-G31755
# TATCTATTGTTTATTCTACAATTAATACAAAGGCAGAAATCAAAAAT
