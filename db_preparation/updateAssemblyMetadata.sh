#!/bin/bash
if [[ $# -ne 3 ]]; then
    echo "$0 db_dir nodes.dmp assembly_dir"
    exit
fi

if [ ! -d ${3} ]
then
    mkdir -p ${3}
fi
cd ${3}

python3 genAssemblyMetadata.py --archaea --nodesDmp ${2} --assemblyLength archaea_assembly_length --assemblyPath archaea_assembly_path --assemblyTaxid archaea_assembly_tax_id --sequenceSummary archaea_sequence_summary --db_dir ${1}
python3 genAssemblyMetadata.py --bacteria --nodesDmp ${2} --assemblyLength bacteria_assembly_length --assemblyPath bacteria_assembly_path --assemblyTaxid bacteria_assembly_tax_id --sequenceSummary bacteria_sequence_summary --db_dir ${1}
python3 genAssemblyMetadata.py --fungi --nodesDmp ${2} --assemblyLength fungi_assembly_length --assemblyPath fungi_assembly_path --assemblyTaxid fungi_assembly_tax_id --sequenceSummary fungi_sequence_summary --db_dir ${1}
python3 genAssemblyMetadata.py --protozoa --nodesDmp ${2} --assemblyLength protozoa_assembly_length --assemblyPath protozoa_assembly_path --assemblyTaxid protozoa_assembly_tax_id --sequenceSummary protozoa_sequence_summary --db_dir ${1}
python3 genAssemblyMetadata.py --vertebrate_mammalian --nodesDmp ${2} --assemblyLength human_assembly_length --assemblyPath human_assembly_path --assemblyTaxid human_assembly_tax_id --sequenceSummary human_sequence_summary --db_dir ${1}
python3 genAssemblyMetadata.py --viral --nodesDmp ${2} --assemblyLength viral_assembly_length --assemblyPath viral_assembly_path --assemblyTaxid viral_assembly_tax_id --sequenceSummary viral_sequence_summary --db_dir ${1}
python3 genAssemblyMetadata.py --plasmid --assemblyLength plasmid_assembly_length --assemblyPath plasmid_assembly_path --sequenceSummary plasmid_sequence_summary --assemblyTaxid plasmid_assembly_tax_id --num 8 --db_dir ${1}

#then cat all together
cat archaea_assembly_length bacteria_assembly_length human_assembly_length viral_assembly_length fungi_assembly_length protozoa_assembly_length > assembly_length
cat archaea_assembly_path bacteria_assembly_path human_assembly_path viral_assembly_path fungi_assembly_path protozoa_assembly_path > assembly_path
cat archaea_assembly_tax_id bacteria_assembly_tax_id human_assembly_tax_id viral_assembly_tax_id fungi_assembly_tax_id protozoa_assembly_tax_id > assembly_tax_id
cat archaea_sequence_summary bacteria_sequence_summary human_sequence_summary viral_sequence_summary fungi_sequence_summary protozoa_sequence_summary > sequence_summary

