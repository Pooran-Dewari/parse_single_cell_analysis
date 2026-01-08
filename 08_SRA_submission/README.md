## Prepare files for SRA submission
- Copy raw fastq files from datastore to eddie scratch
- Check md5sum for fastq file integrity
- Combine fastq files (if same sub-library sequenced on a different lane) using <concat_fastq.sh>; keep log file <concat.log>
- Split fastq files by species, see split_fastq dir
- Generate md5sum hashkeys and upload oyster files to SRA

Transfer files from local to server
```
scp -r fastq_sep_groups_new2024.py  username@eddie.ecdf.ed.ac.uk:<file path>
```
