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
python ${SCRIPT_PATH}/genAssemblySummary.py --db_dir ${DB_DIR} --assemblySummary assembly_summary.csv

#ranks
python ${SCRIPT_PATH}/genRank.py --ranks ranks.csv

#names and nodes
python ${SCRIPT_PATH}/parseDml.py --dmp names.dmp --outputFile abhvfp.names --function 1
python ${SCRIPT_PATH}/parseDml.py --outputFile plasmid.names --function 3 --num 8 
cat abhvfp.names plasmid.names > names.csv && rm abhvfp.names plasmid.names

python ${SCRIPT_PATH}/parseDml.py --dmp nodes.dmp --outputFile abhvfp.nodes --function 2
python ${SCRIPT_PATH}/parseDml.py --outputFile plasmid.nodes --function 4 --num 8
cat abhvfp.nodes plasmid.nodes > nodes.csv && rm abhvfp.nodes plasmid.nodes

cp ${SCRIPT_PATH}/source.csv source.csv

#source.csv table is provided

#sqlite3 script for importing to db
#DELETE FROM table #remove all records from current table

sqlite3 ncbi_taxonomy.db << 'END_SQL'
.mode csv
.separator "\t"
CREATE TABLE assembly_summary(
assembly_id char(20) not null,
bioproject char(20),
biosample char(20),
wgs_master char(20),
refseq_category char(30),
taxid int not null,
species_taxid int not null,
organism_name char(150),
infraspecific_name char(150),
isolate char(150),
version_status char(15),
assembly_level char(20),
release_type char(15),
genome_rep char(15),
seq_rel_date char(10),
asm_name char(150),
submitter char(255),
gbrs_paired_asm char(20),
paired_asm_comp char(20),
ftp_path char(250),
excluded_from_refseq char(100),
relation_to_type_material char(100));
.import assembly_summary.csv assembly_summary
CREATE TABLE sequence_name (sequence_id char(20), sequence_name char(100));
CREATE UNIQUE INDEX idx_sequence_name_sequence_id on sequence_name (sequence_id);
.import sequence_name.csv sequence_name
CREATE TABLE ranks (
        rank VARCHAR NOT NULL,
        height INTEGER NOT NULL,
        PRIMARY KEY (rank),
        UNIQUE (height)
);
.import ranks.csv ranks
create table names(tax_id,tax_name,unique_name,name_class,source_id,is_primary,is_classified);
.import names.csv names
CREATE TABLE nodes (
        tax_id VARCHAR NOT NULL,
        parent_id VARCHAR,
        rank VARCHAR,
        embl_code VARCHAR,
        division_id VARCHAR,
        source_id INTEGER,
        is_valid BOOLEAN,
        PRIMARY KEY (tax_id),
        FOREIGN KEY(rank) REFERENCES ranks (rank),
        FOREIGN KEY(source_id) REFERENCES source (id),
        CHECK (is_valid IN (0, 1))
);
.import nodes.csv nodes
CREATE TABLE source (
        id INTEGER NOT NULL,
        name VARCHAR,
        description VARCHAR,
        PRIMARY KEY (id),
        UNIQUE (name)
);
.import source.csv source
END_SQL
