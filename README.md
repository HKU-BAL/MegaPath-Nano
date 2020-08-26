# MegaPath-Nano
MegaPath-Nano is an Accurate Compositional Analysis and Drug-level Antimicrobial Resistance Detection Software for Oxford Nanopore Long-read Metagenomics

## Database Installation
To use MegaPath-Nano, users need to download RefSeq database and build index first. Script for database preparation is under db_preparation/.
```
#download refseq:
./refseq_download.sh ${DB_DIR}

#download taxnomy to get names.dmp and nodes.dmp:
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar xzf taxdump.tar.gz

#build assembly metadata:
./updateAssemblyMetadata.sh ${DB_DIR} nodes.dmp ${ASSEMBLY_DIR}

#prepare SQL db data:
./updateDB.sh ${DB_DIR} ${SQL_DIR} names.dmp nodes.dmp
#and then follows the sqlite script in updateDB.sh to import data to SQL tables.
```
