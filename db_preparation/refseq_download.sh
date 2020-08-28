#!/bin/bash

SCRIPT=$(readlink -f $0)
ROOT_PATH=$(dirname $(dirname ${SCRIPT}))
DB_DIR=${ROOT_PATH}/genomes/refseq/
if [ $# -eq 1 ]
then
    DB_DIR=$1
fi

if [ ! -d ${DB_DIR} ]; then
    mkdir -p ${DB_DIR}
fi

python3 refseq_download.py --bacteria --get_summary --db_dir ${DB_DIR}
python3 refseq_download.py --fungi --get_summary --db_dir ${DB_DIR}
python3 refseq_download.py --protozoa --get_summary --db_dir ${DB_DIR}
python3 refseq_download.py --viral --get_summary --db_dir ${DB_DIR}
python3 refseq_download.py --archaea --get_summary --db_dir ${DB_DIR}
python3 refseq_download.py --vertebrate_mammalian --get_summary --db_dir ${DB_DIR}
python3 refseq_download.py --plasmid --num 8 --db_dir ${DB_DIR} #8 = number of plasmid file can be downloaded on ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/
