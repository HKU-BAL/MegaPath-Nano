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

## Usage
For all available options, please check [Usage.md](docs/Usage.md)
```
python megapath_nano.py --query ${FASTQ}

General options:
    --query
        Query file (fastq or fasta)
    --human
        Human genome set in config folder, default human.genome_set.
    --decoy
        Decoy genome set in config folder, default plasmid.genome_set.
    --species
        Genome set for species identification in config folder, default species_id.genome_set.
    --assembly
        Genome set for assembly identification in config folder, default assembly_id.genome_set.
    --taxonomy_db
        Taxonomy database, default 'db/ncbi_taxonomy.db'.
    --tool_folder
        Tool folder, default 'tools/'.
    --config_folder
        Config folder, default 'config/'.
    --assembly_folder
        Assembly folder, default 'genomes/'.
    --aligner
        Path to aligner, default 'minima2'.
    --max_aligner_thread INT
        Maximum number of threads used by aligner, default 64.
    --output_prefix
        Output Prefix, query file name will be used for output prefix by default.
    --output_folder
        Output folder, default ./.
```
