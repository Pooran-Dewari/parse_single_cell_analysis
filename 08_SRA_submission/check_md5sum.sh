# this script requires sub-libraries in different directories, each dir should have fastq files and a corresponding MD5.txt
# for example, A1 A2 A3 A4.... sub-dirs

# tree structure of A1 shown below
# A1_EKDL240002473-1A_223M7CLT4_L7_1.fq.gz
# A1_EKDL240002473-1A_223M7CLT4_L7_2.fq.gz
# A1_EKDL240002473-1A_223M7CLT4_L8_1.fq.gz
# A1_EKDL240002473-1A_223M7CLT4_L8_2.fq.gz
# MD5.txt

# the script will ouput md5sum check results in a .txt file called "results_md5sum.txt"

RES="results_md5sum.txt"
touch $RES

ls -d */ | while read DIR;
do
 cd $DIR
 echo "working on $DIR" >> ../$RES 2>&1
 md5sum -c MD5.txt >> ../$RES 2>&1
 echo "done $DIR !" >> ../$RES 2>&1
 cd ..
done
