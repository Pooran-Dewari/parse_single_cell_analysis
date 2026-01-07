## Prepare files for SRA submission
- Copy files from datastore to scratch
- check md5sum for fastq file integrity
- combine fastq files (if same sub-library sequenced on a different lane) using concat_fastq.sh; keep log
- split groups (into oyster and shrimp) using python script fastq_sep_groups_new2024.py & split_fastq_2026.sh

```
scp -r fastq_sep_groups_new2024.py  username@eddie.ecdf.ed.ac.uk:<file path>
```
