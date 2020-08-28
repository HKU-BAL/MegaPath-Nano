#!/bin/bash
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
ROOT_PATH=$(dirname $(dirname ${SCRIPT}))
DB_DIR=${ROOT_PATH}/genomes/refseq/
ASSEMBLY_DIR=${ROOT_PATH}/genomes/

if [ $# -eq 2 ]
then
    DB_DIR=$1
    ASSEMBLY_DIR=$2
fi

if [ ! -d ${DB_DIR} ]; then
    echo "DB not exist. Please use refseq_download.sh to download RefSeq DB first."
    exit
fi

if [ ! -d ${ASSEMBLY_DIR} ]; then
    mkdir -p ${ASSEMBLY_DIR}
fi

cd ${ASSEMBLY_DIR}

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar xzf taxdump.tar.gz && rm citations.dmp delnodes.dmp division.dmp gencode.dmp merged.dmp readme.txt

python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --archaea --nodesDmp nodes.dmp --assemblyLength archaea_assembly_length --assemblyPath archaea_assembly_path --assemblyTaxid archaea_assembly_tax_id --sequenceSummary archaea_sequence_summary --db_dir ${1}
python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --bacteria --nodesDmp nodes.dmp --assemblyLength bacteria_assembly_length --assemblyPath bacteria_assembly_path --assemblyTaxid bacteria_assembly_tax_id --sequenceSummary bacteria_sequence_summary --db_dir ${1}
python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --fungi --nodesDmp nodes.dmp --assemblyLength fungi_assembly_length --assemblyPath fungi_assembly_path --assemblyTaxid fungi_assembly_tax_id --sequenceSummary fungi_sequence_summary --db_dir ${1}
python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --protozoa --nodesDmp nodes.dmp --assemblyLength protozoa_assembly_length --assemblyPath protozoa_assembly_path --assemblyTaxid protozoa_assembly_tax_id --sequenceSummary protozoa_sequence_summary --db_dir ${1}
python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --vertebrate_mammalian --nodesDmp nodes.dmp --assemblyLength human_assembly_length --assemblyPath human_assembly_path --assemblyTaxid human_assembly_tax_id --sequenceSummary human_sequence_summary --db_dir ${1}
python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --viral --nodesDmp nodes.dmp --assemblyLength viral_assembly_length --assemblyPath viral_assembly_path --assemblyTaxid viral_assembly_tax_id --sequenceSummary viral_sequence_summary --db_dir ${1}
python3 ${SCRIPT_PATH}/genAssemblyMetadata.py --plasmid --assemblyLength plasmid_assembly_length --assemblyPath plasmid_assembly_path --sequenceSummary plasmid_sequence_summary --assemblyTaxid plasmid_assembly_tax_id --num 8 --db_dir ${1}

#then cat all together
cat archaea_assembly_length bacteria_assembly_length human_assembly_length viral_assembly_length fungi_assembly_length protozoa_assembly_length > assembly_length && rm archaea_assembly_length bacteria_assembly_length human_assembly_length viral_assembly_length fungi_assembly_length protozoa_assembly_length
cat archaea_assembly_path bacteria_assembly_path human_assembly_path viral_assembly_path fungi_assembly_path protozoa_assembly_path > assembly_path && rm archaea_assembly_path bacteria_assembly_path human_assembly_path viral_assembly_path fungi_assembly_path protozoa_assembly_path
cat archaea_assembly_tax_id bacteria_assembly_tax_id human_assembly_tax_id viral_assembly_tax_id fungi_assembly_tax_id protozoa_assembly_tax_id > assembly_tax_id && rm archaea_assembly_tax_id bacteria_assembly_tax_id human_assembly_tax_id viral_assembly_tax_id fungi_assembly_tax_id protozoa_assembly_tax_id
cat archaea_sequence_summary bacteria_sequence_summary human_sequence_summary viral_sequence_summary fungi_sequence_summary protozoa_sequence_summary > sequence_summary && rm archaea_sequence_summary bacteria_sequence_summary human_sequence_summary viral_sequence_summary fungi_sequence_summary protozoa_sequence_summary

