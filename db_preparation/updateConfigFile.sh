#!/bin/bash
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
ROOT_PATH=$(dirname $(dirname ${SCRIPT}))
DB_DIR=${ROOT_PATH}/genomes/refseq/
CONFIG_DIR=${ROOT_PATH}/config/

if [ $# -eq 2 ]
then
    DB_DIR=$1
    CONFIG_DIR=$2
fi

if [ ! -d ${DB_DIR} ]; then
    echo "DB not exist. Please use refseq_download.sh to download RefSeq DB first."
    exit
fi

if [ ! -d ${CONFIG_DIR} ]; then
    mkdir -p ${CONFIG_DIR}
fi

cd ${CONFIG_DIR}

python3 ${SCRIPT_PATH}/genConfigFile.py --assemblySummary ${DB_DIR}/archaea/assembly_summary.txt --outputFile species_id.archaea.genome_set --outputFile2 assembly_id.archaea.genome_set --function 1
python3 ${SCRIPT_PATH}/genConfigFile.py --assemblySummary ${DB_DIR}/bacteria/assembly_summary.txt --outputFile species_id.bacteria.genome_set --outputFile2 assembly_id.bacteria.genome_set --function 1
python3 ${SCRIPT_PATH}/genConfigFile.py --assemblySummary ${DB_DIR}/protozoa/assembly_summary.txt --outputFile species_id.protozoa.genome_set --outputFile2 assembly_id.protozoa.genome_set --function 1
python3 ${SCRIPT_PATH}/genConfigFile.py --assemblySummary ${DB_DIR}/fungi/assembly_summary.txt --outputFile species_id.fungi.genome_set --outputFile2 assembly_id.fungi.genome_set --function 1
python3 ${SCRIPT_PATH}/genConfigFile.py --assemblySummary ${DB_DIR}/viral/assembly_summary.txt --outputFile species_id.viral.genome_set --outputFile2 assembly_id.viral.genome_set --function 1

python3 ${SCRIPT_PATH}/genConfigFile.py --assemblySummary ${DB_DIR}/vertebrate_mammalian/assembly_summary.txt --outputFile human.genome_set --function 2 #only the first column is used
python3 ${SCRIPT_PATH}/genConfigFile.py --num 8 --outputFile plasmid.genome_set --function 3

#same set of genomes for global selection (species) and specific selection (assembly)
cat species_id.archaea.genome_set species_id.bacteria.genome_set species_id.protozoa.genome_set species_id.fungi.genome_set species_id.viral.genome_set assembly_id.archaea.genome_set assembly_id.bacteria.genome_set assembly_id.protozoa.genome_set assembly_id.fungi.genome_set assembly_id.viral.genome_set > species_id.genome_set && rm species_id.archaea.genome_set species_id.bacteria.genome_set species_id.protozoa.genome_set species_id.fungi.genome_set species_id.viral.genome_set assembly_id.archaea.genome_set assembly_id.bacteria.genome_set assembly_id.protozoa.genome_set assembly_id.fungi.genome_set assembly_id.viral.genome_set
cp species_id.genome_set assembly_id.genome_set
