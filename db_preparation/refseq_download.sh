if [[ $# -ne 1 ]]; then
    echo "$0 db_dir"
    exit
fi

python3 refseq_download.py --bacteria --get_summary --db_dir ${1}
python3 refseq_download.py --fungi --get_summary --db_dir ${1}
python3 refseq_download.py --protozoa --get_summary --db_dir ${1}
python3 refseq_download.py --viral --get_summary --db_dir ${1}
python3 refseq_download.py --archaea --get_summary --db_dir ${1}
python3 refseq_download.py --vertebrate_mammalian --get_summary --db_dir ${1}
python3 refseq_download.py --plasmid --num 8 --db_dir ${1} #8 = number of plasmid file can be downloaded on ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/
