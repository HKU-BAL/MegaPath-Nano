[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/HKU-BAL/MegaPath-Nano/master)

# MegaPath-Nano

## Introduction

The ultra-long ONT sequencing technology benefits metagenomic profiling with high alignment specificity. Yet, its high sequencing error per read remains a hurdle to distinguish among closely related pathogens at lower taxonomic ranks, and for refined drug-level antimicrobial resistance prediction. In this study, we present MegaPath-Nano, successor to the NGS-based MegaPath, an accurate compositional analysis software with drug-level AMR identification for ONT metagenomic sequencing data. MegaPath-Nano takes ONT raw reads as input, and performs  data cleansing, taxonomic profiling, and drug-level AMR detection within a single workflow. The major output of our tool includes 1) a taxonomic profiling report down to strain level with abundance estimated; and 2) an integrated class and drug level AMR report in tabular format with supportive information from different detection tools. As a key feature for taxonomic profiling, MegaPath-Nano performs a global-optimization on multiple alignments and reassigns predictably misplaced reads to a single most likely species. To perform a consistent and comprehensive AMR detection analysis, MegaPath-Nano uses a novel consensus-based approach to detect AMR, incorporating a collection of AMR software and databases. We benchmarked against other state-of-the-art software, including WIMP, Kraken 2, MetaMaps, ARMA and ARGpore using real sequencing data, and we achieved the best performance in both tasks. MegaPath-Nano is therefore a well rounded ONT metagenomic tool for clinical use in practice.

## Prerequisites

## Option 1: Bioconda
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mpn -c bioconda megapath-nano
conda activate mpn
```


## Option 2: Conda Virtual Environment Setup
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mpn python=3.6
conda activate mpn

# installing all dependencies for both modules
conda install pandas==1.1.5 psutil==5.6.5 pybedtools==0.8.0 porechop==0.2.4 bioconvert==0.4.3 seqtk==1.3 minimap2==2.21 bcftools==1.9 samtools==1.9 pysam==0.16.0 tabulate==0.8.9 cgecore==1.5.6 ncbi-amrfinderplus==3.10.5 rgi==5.2.0 biopython==1.72 pyahocorasick==1.1.7 filetype==1.0.7libdeflate==1.6

# git clone MegaPath-Nano
git clone --depth 1 https://github.com/HKU-BAL/MegaPath-Nano
```

## Option 3: Docker
```
sudo docker build -f ./Dockerfile -t mpn_image . 
sudo docker run -it mpn_image /bin/bash
```



## Preprocessed Database Download
```
# Option 1, Bioconda: cd ${CONDA_PREFIX}/MegaPath-Nano
# conda info --env can show the ${CONDA_PREFIX} in the current environment.
# Option 2, Conda Virtual Env: cd ./MegaPath-Nano (the git clone)
# Option 3, Docker: cd /opt/MegaPath-Nano
cd ${MEGAPATH_NANO_DIR}

# Taxon
wget http://www.bio8.cs.hku.hk/dataset/MegaPath-Nano/MegaPath-Nano_db.v1.0.tar.gz
tar -xvzf MegaPath-Nano_db.v1.0.tar.gz

# AMR
rgi load --card_json bin/amr_db/card/card.json
amrfinder -u
```

## Alternative: Online Database Installation
The latest RefSeq database can be downloaded with the scripts under db_preparation/. 
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

# (optional) add custom FASTA sequences to the decoy database 
python addDecoyDB.py --decoy_fasta ${fasta}

# AMR
# prepare AMR databases:
./prepareAMR_DB.sh

```

## Basic usage
(1) Run taxonomic analysis and AMR deteciton module
```
python megapath_nano.py --query ${fq/fa} [options]

required arguments:
  --query
                              Query file (fastq or fasta)

optional arguments:
  --max_aligner_thread INT    Maximum number of threads used by aligner, default: 64. Actual number of threads is min( available num of cores, threads specified)
  --output_prefix             Output Prefix, default: query file name
  --output_folder             Output folder, default: current working directory 
```

(2) Run taxonomic analysis module only
```
python megapath_nano.py --query ${fq/fa} --taxon_module_only [options]

```

(3) Run AMR deteciton module only with **FASTQ/FASTA**
```
python megapath_nano.py --query ${fq/fa} --AMR_module_only [options]

```

(4) Filter FQ/FA only: Adaptor trimming, read filtering and trimming, human or decoy filtering
```
python megapath_nano.py --query ${fq/fa} --filter_fq_only [options]
```
For all available options, please check [Usage.md](docs/Usage.md)


(5) Run AMR deteciton module only with **BAM**
```
python megapath_nano_amr.py --query_bam ${bam} --output_folder ${dir} [options]

required arguments:
  --query_bam QUERY_BAM
                              Input bam
  --output_folder OUTPUT_FOLDER
                              Output directory

optional arguments:
  --taxon TAXON               Taxon-specific options for AMRFinder [e.g. --taxon Escherichia], see usage for the full list of curated organisms
  --threads THREADS           Max num of threads, default: available num of cores
```



## Demo data

The demo data for AMR detection of five patient isolates are available for download on http://www.bio8.cs.hku.hk/dataset/MegaPath-Nano/. Samples were prepared using ONT Rapid Sequencing Kit, and sequenced using ONT R9.4.1 flowcells.
  
The experimental validation results of these AMR demo datasets are listed on [Supplementary_info_AMR](docs/Supplementary_info_demo_AMR_data.md).

## Demo run

```
wget http://www.bio8.cs.hku.hk/dataset/MegaPath-Nano/Escherichia_coli_isolate2_HKUBAL_20200103.fastq
python megapath_nano.py --query Escherichia_coli_isolate2_HKUBAL_20200103.fastq
```
