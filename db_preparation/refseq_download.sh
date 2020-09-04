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

nohup python3 refseq_download.py --bacteria --get_summary --db_dir ${DB_DIR} > nohup.bacteria_refseq_download &
nohup python3 refseq_download.py --fungi --get_summary --db_dir ${DB_DIR} > nohup.fungi_refseq_download &
nohup python3 refseq_download.py --protozoa --get_summary --db_dir ${DB_DIR} > nohup.protozoa_refseq_download &
nohup python3 refseq_download.py --viral --get_summary --db_dir ${DB_DIR} > nohup.viral_refseq_download &
nohup python3 refseq_download.py --archaea --get_summary --db_dir ${DB_DIR} > nohup.archaea_refseq_download &
nohup python3 refseq_download.py --vertebrate_mammalian --get_summary > nohup.vertebrate_mammalian_refseq_download &
#8 = number of plasmid file can be downloaded on ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/--db_dir ${DB_DIR}
nohup python3 refseq_download.py --plasmid --num 8 --db_dir ${DB_DIR} > nohup.plasmid_refseq_download & 
