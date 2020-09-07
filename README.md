[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# MegaPath-Nano

## Introduction

The ultra-long ONT sequencing technology benefits metagenomic profiling with high alignment specificity. Yet, its high sequencing error per read remains a hurdle to distinguish among closely related pathogens at lower taxonomic ranks, and for refined drug-level antimicrobial resistance prediction. In this study, we present MegaPath-Nano, successor to the NGS-based MegaPath, an accurate compositional analysis software with drug-level AMR identification for ONT metagenomic sequencing data. MegaPath-Nano takes ONT raw reads as input, and performs  data cleansing, taxonomic profiling, and drug-level AMR detection within a single workflow. The major output of our tool includes 1) a taxonomic profiling report down to strain level with abundance estimated; and 2) an integrated class and drug level AMR report in tabular format with supportive information from different detection tools. As a key feature for taxonomic profiling, MegaPath-Nano performs a global-optimization on multiple alignments and reassigns predictably misplaced reads to a single most likely species. To perform a consistent and comprehensive AMR detection analysis, MegaPath-Nano uses a novel consensus-based approach to detect AMR, incorporating a collection of AMR software and databases. We benchmarked against other state-of-the-art software, including WIMP, Kraken 2, MetaMaps, ARMA and ARGpore using real sequencing data, and we achieved the best performance in both tasks. MegaPath-Nano is therefore a well rounded ONT metagenomic tool for clinical use in practice.

## Prerequisites
Tested OS environment: Ubuntu 16.04

## Conda Virtual Environment Setup
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mpn python=3.6
conda activate mpn

# installing all dependencies for both modules
conda install pandas==0.23 psutil pybedtools qcat bioconvert seqtk minimap2 bcftools samtools cgecore pysam tabulate rgi ncbi-amrfinderplus
pip install pandarallel

```

## Git clone MegaPath-Nano
```
git clone --depth 1 https://github.com/HKU-BAL/MegaPath-Nano
cd MegaPath-Nano
```

## Database Installation
To use MegaPath-Nano, users need to download RefSeq database and build index first. Script for database preparation is under db_preparation/. 
```
# Taxon
# download RefSeq:
./refseq_download.sh [${DB_DIR}=MegaPath-Nano/genomes/refseq/]

# build assembly metadata:
./updateAssemblyMetadata.sh [${DB_DIR}=MegaPath-Nano/genomes/refseq/] [${ASSEMBLY_DIR}=MegaPath-Nano/genomes/]

# generate config files:
./updateConfigFile.sh [${DB_DIR}=MegaPath-Nano/genomes/refseq/] [${CONFIG_DIR}=MegaPath-Nano/config/]

# prepare SQL db data:
./updateDB.sh [${DB_DIR}=MegaPath-Nano/genomes/refseq/] [${SQL_DIR}=MegaPath-Nano/db/]
# and then follows the sqlite script in updateDB.sh to import data to SQL tables.


# AMR
# prepare AMR databases:
./prepareAMR_DB.sh

```

## Basic usage
(1) Run taxonomic analysis and AMR deteciton module
```
python MegaPath-Nano_taxon.py --query ${fq/fa} [options]

required arguments:
  --query
                              Query file (fastq or fasta)

optional arguments:
  --max_aligner_thread INT    Maximum number of threads used by aligner, default 64. Actual number of threads is min( available num of cores, threads specified).
  --output_prefix             Output Prefix, query file name will be used for output prefix by default.
  --output_folder             Output folder, default ./.
```
For all available options, please check [Usage.md](docs/Usage.md)

(2) Run taxonomic analysis module only
```
python MegaPath-Nano_taxon.py --query ${fq/fa} [options] --no_AMR_module

required arguments:
  --query
                              Query file (fastq or fasta)

```

(3) Run AMR deteciton module only
```
python MegaPath-Nano_AMR.py --query_bam ${bam} --output_folder ${dir} [options]

required arguments:
  --query_bam QUERY_BAM
                              Input bam
  --output_folder OUTPUT_FOLDER
                              Output directory

optional arguments:
  --taxon TAXON               Taxon-specific options for AMRFinder [e.g. --taxon Escherichia], see usage for the full list of curated organisms
  --threads THREADS           Max num of threads
```

(4) Filter FQ/FA only: Adaptor trimming, read filtering and trimming, human or decoy filtering
```
python MegaPath-Nano_taxon.py --query ${fq/fa} --filter_fq_only [options]
```

## Demo data

The demo data for AMR detection of five patient isolates are available for download on http://www.bio8.cs.hku.hk/dataset/MegaPath-Nano/. Samples were prepared using ONT Rapid Sequencing Kit, and sequenced using ONT R9.4.1 flowcells.
  
The experimental validation results of these AMR demo datasets are listed on [Supplementary_info_AMR](docs/Supplementary_info_demo_AMR_data.md).
