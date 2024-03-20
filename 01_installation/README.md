Official documentation from Parse Biosciences [here](https://support.parsebiosciences.com/hc/en-us/articles/23060102930580-Pipeline-Installation-Current-Version), login required.

### create a conda environment for installation
```
# requesting node to run conda
qlogin -l h_vmem=20G
module load anaconda/2024.02
#create a new conda environment 'spipe', and install all required dependencies in it
conda create -n spipe conda-forge::python==3.10
#the command above will take a while
```

### create directory structure, download parse pipeline, and install
```
conda activate spipe
cd /exports/eddie/scratch/pdewari/
#create directory structure for download
mkdir newvolume
cd newvolume/
mkdir analysis expdata genomes
wget -c https://support.parsebiosciences.com/hc/en-us/article_attachments/23931973337876
unzip ParseBiosciences-Pipeline.1.2.0.zip
cd ParseBiosciences-Pipeline.1.2.0/
bash ./install_dependencies_conda.sh -i -y
pip install --no-cache-dir ./
```
### directory structure after the install
```
tree -d
```
It should output something like as below
```
.
├── analysis
├── expdata
├── genomes
└── ParseBiosciences-Pipeline.1.2.0
    ├── build
    │   ├── bdist.linux-x86_64
    │   ├── lib
    │   │   └── splitpipe
    │   │       ├── barcodes
    │   │       ├── config
    │   │       ├── scripts
    │   │       │   ├── config
    │   │       │   └── db
    │   │       └── templates
    │   └── scripts-3.10
    ├── splitpipe
    │   ├── barcodes
    │   ├── config
    │   ├── scripts
    │   │   ├── config
    │   │   └── db
    │   └── templates
    └── splitpipe.egg-info

23 directories

```
### check if installed correctly
```
split-pipe --help
```
If installed correctly, it should print this help message on the screen

```
usage: split-pipe [-h] [-m MODE] [-c CHEMISTRY] [--kit KIT] [-p PARFILE] [--run_name RUN_NAME] [--fq1 FQ1] [--fq2 FQ2] [--output_dir OUTPUT_DIR] [--genome_dir GENOME_DIR]
                  [--parent_dir PARENT_DIR] [--targeted_list TARGETED_LIST] [--sample SAMPLE_NAME WELLS] [--samp_list SAMP_LIST] [--samp_sltab SAMP_SLTAB] [--genome_name [GENOME_NAME ...]]
                  [--genes [GENES ...]] [--fasta [FASTA ...]] [--gfasta GENOME_NAME FASTA] [--sublibraries [SUBLIBRARIES ...]] [--sublib_list SUBLIB_LIST] [--sublib_pref SUBLIB_PREF]
                  [--sublib_suff SUBLIB_SUFF] [--tscp_use TSCP_USE] [--tscp_min TSCP_MIN] [--tscp_max TSCP_MAX] [--cell_use CELL_USE] [--cell_est CELL_EST] [--cell_xf CELL_XF]
                  [--cell_min CELL_MIN] [--cell_max CELL_MAX] [--cell_list CELL_LIST] [--focal_barcoding] [--crispr] [--crsp_guides CRSP_GUIDES] [--crsp_read_thresh CRSP_READ_THRESH]
                  [--crsp_tscp_thresh CRSP_TSCP_THRESH] [--crsp_max_mm] [--crsp_use_star] [--tcr_analysis] [--tcr_genome TCR_GENOME] [--tcr_check] [--use_imgt_db] [--no_save_anndata]
                  [--kit_list] [--chem_list] [--bc_list] [--bc_round_set ROUND NAME] [--rseed RSEED] [--nthreads NTHREADS] [--no_keep_going] [--reuse] [--keep_temps] [--one_step]
                  [--until_step UNTIL_STEP] [--clear_runproc] [--start_timeout START_TIMEOUT] [--kit_score_skip] [--dryrun] [-e] [-V]

SplitPipe data processing pipeline v1.2.0
```

