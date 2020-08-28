[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# MegaPath-Nano

## Introduction

The ultra-long ONT sequencing technology benefits metagenomic profiling with high alignment specificity. Yet, its high sequencing error per read remains a hurdle to distinguish among closely related pathogens at lower taxonomic ranks, and for refined drug-level antimicrobial resistance prediction. In this study, we present MegaPath-Nano, successor to the NGS-based MegaPath, an accurate compositional analysis software with drug-level AMR identification for ONT metagenomic sequencing data. MegaPath-Nano takes ONT raw reads as input, and performs  data cleansing, taxonomic profiling, and drug-level AMR detection within a single workflow. The major output of our tool includes 1) a taxonomic profiling report down to strain level with abundance estimated; and 2) an integrated class and drug level AMR report in tabular format with supportive information from different detection tools. As a key feature for taxonomic profiling, MegaPath-Nano performs a global-optimization on multiple alignments and reassigns predictably misplaced reads to a single most likely species. To perform a consistent and comprehensive AMR detection analysis, MegaPath-Nano uses a novel consensus-based approach to detect AMR, incorporating a collection of AMR software and databases. We benchmarked against other state-of-the-art software, including WIMP, Kraken 2, MetaMaps, ARMA and ARGpore using real sequencing data, and we achieved the best performance in both tasks. MegaPath-Nano is therefore a well rounded ONT metagenomic tool for clinical use in practice.

## Prerequisites

## Conda Virtual Environment Setup
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mpn python=3.6
conda activate mpn

# composition analysis module's dependencies
conda install pandas==0.23 psutil pybedtools qcat bioconvert

# AMR detection module's dependencies
conda install bcftools samtools cgecore pysam conda-forge tabulate
```

## Git clone MegaPath-Nano
```
git clone --depth 1 https://github.com/HKU-BAL/MegaPath-Nano
cd MegaPath-Nano
```

## Database Installation
To use MegaPath-Nano, users need to download RefSeq database and build index first. Script for database preparation is under db_preparation/.
```
cd db_preparation

# download RefSeq:
./refseq_download.sh ${DB_DIR}

# download taxnomy to get names.dmp and nodes.dmp:
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar xzf taxdump.tar.gz

# build assembly metadata:
./updateAssemblyMetadata.sh ${DB_DIR} nodes.dmp ${ASSEMBLY_DIR}

# generate config files:
./updateConfigFile.sh ${DB_DIR} ${CONFIG_DIR}

# prepare SQL db data:
./updateDB.sh ${DB_DIR} ${SQL_DIR} names.dmp nodes.dmp
# and then follows the sqlite script in updateDB.sh to import data to SQL tables.
```
## Usage
For all available options, please check [Usage.md](docs/Usage.md)
```
python megapath_nano.py --query ${FASTQ}

Required Arguments:
    --query
        Query file (fastq or fasta)

Optional Arguments:
    --aligner
        Path to aligner, default 'minimap2'.
    --max_aligner_thread INT
        Maximum number of threads used by aligner, default 64.
    --output_prefix
        Output Prefix, query file name will be used for output prefix by default.
    --output_folder
        Output folder, default ./.
    -h, --help
        Show help message and exit
```

## Demo data

Demo data for AMR detection of five patient isolates are available at:
http://www.bio8.cs.hku.hk/dataset/MegaPath-Nano/.
Samples were prepared using ONT Rapid Sequencing Kit,and sequenced using ONT R9.4.1 flowcells.
