#!/bin/bash
if [[ $# -ne 2 ]]; then
    echo "$0 db_dir config_dir"
    exit
fi

if [ ! -d ${2} ]
then
    mkdir -p ${2}
fi
cd ${2}

python3 genConfigFile.py --assemblySummary ${1}/archaea/assembly_summary.txt --outputFile species_id.archaea.genome_set --outputFile2 assembly_id.archaea.genome_set --function 1
python3 genConfigFile.py --assemblySummary ${1}/bacteria/assembly_summary.txt --outputFile species_id.bacteria.genome_set --outputFile2 assembly_id.bacteria.genome_set --function 1
python3 genConfigFile.py --assemblySummary ${1}/protozoa/assembly_summary.txt --outputFile species_id.protozoa.genome_set --outputFile2 assembly_id.protozoa.genome_set --function 1
python3 genConfigFile.py --assemblySummary ${1}/fungi/assembly_summary.txt --outputFile species_id.fungi.genome_set --outputFile2 assembly_id.fungi.genome_set --function 1
python3 genConfigFile.py --assemblySummary ${1}/viral/assembly_summary.txt --outputFile species_id.viral.genome_set --outputFile2 assembly_id.viral.genome_set --function 1

python3 genConfigFile.py --assemblySummary ${1}/vertebrate_mammalian/assembly_summary.txt --outputFile human.genome_set --function 2 #only first column is used
python3 genConfigFile.py --num 8 --outputFile plasmid.genome_set --function 3

cat species_id.archaea.genome_set species_id.bacteria.genome_set species_id.protozoa.genome_set species_id.fungi.genome_set species_id.viral.genome_set > species_id.genome_set
cat assembly_id.archaea.genome_set assembly_id.bacteria.genome_set assembly_id.protozoa.genome_set assembly_id.fungi.genome_set assembly_id.viral.genome_set > assembly_id.genome_set
