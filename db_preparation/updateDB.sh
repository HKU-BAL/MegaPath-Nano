#!/bin/bash
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
ROOT_PATH=$(dirname $(dirname ${SCRIPT}))
DB_DIR=${ROOT_PATH}/genomes/refseq/
SQL_DB_DIR=${ROOT_PATH}/db/

if [ $# -eq 2 ]
then
    DB_DIR=$1
    SQL_DB_DIR=$2
fi

if [ ! -d ${DB_DIR} ]; then
    echo "DB not exist. Please use refseq_download.sh to download RefSeq DB first."
    exit
fi

if [ ! -d ${SQL_DB_DIR} ]; then
    mkdir -p ${SQL_DB_DIR}
fi

cd ${SQL_DB_DIR}

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar xzf taxdump.tar.gz && rm citations.dmp delnodes.dmp division.dmp gencode.dmp merged.dmp readme.txt

#sequence_name table
python ${SCRIPT_PATH}/genSequenceName.py --function 1 --sequenceName abhvfp.sequence_name --db_dir ${DB_DIR}
python ${SCRIPT_PATH}/genSequenceName.py --function 2 --sequenceName plasmid.sequence_name --db_dir ${DB_DIR} --num 8
cat abhvfp.sequence_name plasmid.sequence_name > sequence_name.csv && rm abhvfp.sequence_name plasmid.sequence_name

#assembly_summary table
python ${SCRIPT_PATH}/genAssemblySummary.py --db_dir ${DB_DIR} --assemblySummary assembly_summary

#ranks
python ${SCRIPT_PATH}/genRank.py --ranks ranks.csv

#names and nodes
python ${SCRIPT_PATH}/parseDml.py --dmp names.dmp --outputFile abhvfp.names --function 1
python ${SCRIPT_PATH}/parseDml.py ----outputFile plasmid.names --function 3 --num 8 
cat abhvfp.names plasmid.names > names.csv && rm abhvfp.names plasmid.names

python ${SCRIPT_PATH}/parseDml.py --dmp nodes.dmp --outputFile abhvfp.nodes --function 2
python ${SCRIPT_PATH}/parseDml.py ----outputFile plasmid.nodes --function 4 --num 8
cat abhvfp.nodes plasmid.nodes > nodes.csv && rm abhvfp.nodes plasmid.nodes

cp ${SCRIPT_PATH}/source.csv source.csv

#source.csv table is provided

#sqlite3 script for importing to db
#DELETE FROM table #remove all records from current table

sqlite3 ncbi_taxonomy.db << 'END_SQL'
.mode csv
.separator "\t"
.import assembly_summary.csv assembly_summary
.import sequence_name.csv sequence_name
.import ranks.csv ranks
.import names.csv names
.import nodes.csv nodes
.import source.csv source
END_SQL
