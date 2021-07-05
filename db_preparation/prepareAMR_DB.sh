#!/bin/bash
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
ROOT_PATH=$(dirname $(dirname ${SCRIPT}))

cd ${ROOT_PATH}/bin
#TODO update to 4.0
git clone -b 3.2.1 https://git@bitbucket.org/genomicepidemiology/resfinder.git

cd amr_db

mkdir card
wget https://card.mcmaster.ca/latest/data -P card
tar -C card -xvf card/data ./card.json
rgi load --card_json card/card.json

git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git resfinder

mkdir megares
wget https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta -O megares/megares_full_database_v2.00.fsa

wget -r -np -nd -R "*.html*" http://proteininformatics.org/mkumar/lactamasedb/downloadnuc/ -P cbmar
awk 'FNR==1{print ""}{print}' cbmar/*tide.fasta > cbmar/cbmar_nucl.fsa && rm -f cbmar/*tide.fasta
wget -r -np -nd -R "*.html*" http://proteininformatics.org/mkumar/lactamasedb/download_protein/ -P cbmar
awk 'FNR==1{print ""}{print}' cbmar/*.fasta > cbmar/cbmar_prot.fsa && rm -f cbmar/*.fasta
makeblastdb -in cbmar/cbmar_prot.fsa   -title "cbmar_prot" -dbtype prot

#download AMRFinder database
amrfinder -u
