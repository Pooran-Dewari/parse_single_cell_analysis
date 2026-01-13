## Prepare data for SRA submission
- Step 1 : Copy raw fastq files from datastore to eddie scratch and check md5sum for fastq file integrity
- Step 2: All our sub-libraries were sequenced on multiple lanes to gain more coverage. Combine these fastq files using <concat_fastq.sh>; keep log file <concat.log>
- Step 3: Ambre, Alex, and James had prepared Parse libraries. Nuclei from both oyster and shrimp were processed together. We need to split fastq files by species prior to mapping to oyster genome. For splitting fastq by species, see split_fastq dir.
- Step 4: Upload oyster files to SRA using sftp username@sftp-private.ncbi.nlm.nih.gov; use put *.fastq.gz; keep a log of md5sum for these files.

Transfer files from local to server
```
scp -r fastq_sep_groups_new2024.py  username@eddie.ecdf.ed.ac.uk:<file path>
```
