#!/bin/bash
if [[ $# -ne 4 ]]; then
    echo "$0 db_dir sql_dir names.dmp nodes.dmp"
    exit
fi

if [ ! -d ${2} ]
then
    mkdir -p ${2}
fi
cd ${2}

#sequence_name table
python genSequenceName.py --function 1 --sequenceName abhvfp.sequence_name --db_dir ${1}
python genSequenceName.py --function 2 --sequenceName plasmid.sequence_name --db_dir ${1} --num 8
cat abhvfp.sequence_name plasmid.sequence_name > sequence_name.csv

#assembly_summary table
python genAssemblySummary.py --db_dir ${1} --assemblySummary assembly_summary

#ranks
python genRank.py --ranks ranks.csv

#names and nodes
python parseDml.py --dmp ${3} --outputFile abhvfp.names --function 1
python parseDml.py ----outputFile plasmid.names --function 3 --num 8 
cat abhvfp.names plasmid.names > names.csv

python parseDml.py --dmp ${4} --outputFile abhvfp.nodes --function 2
python parseDml.py ----outputFile plasmid.nodes --function 4 --num 8
cat abhvfp.nodes plasmid.nodes > nodes.csv

#source.csv table is provided

#sqlite3 script
#import to db
#DELETE FROM table #remove all records from current table

#sqlite>.mode csv
#sqlite>.separator "\t"
#sqlite>.import assembly_summary.csv assembly_summary
#sqlite>.import sequence_name.csv sequence_name
#sqlite>.import ranks.csv ranks
#sqlite>.import names.csv names
#sqlite>.import nodes.csv nodes
#sqlite>.import source.csv source

