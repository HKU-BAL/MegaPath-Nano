import argparse
from shutil import copyfile
import os
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add fasta to decoy db')
    parser.add_argument('--decoy_fasta', required=True)

    cwd=os.path.dirname(os.path.realpath(__file__))
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--assembly_dir', default=f'{NANO_DIR}/genomes' )
    parser.add_argument('--config_folder', help='Config file folder', default=f'{NANO_DIR}/config' )

    FLAGS = parser.parse_args()

    decoy_name=os.path.basename(os.path.splitext(FLAGS.decoy_fasta)[0])
    path= f'{FLAGS.assembly_dir}/refseq/plasmid/{decoy_name}.fa'
    copyfile(FLAGS.decoy_fasta,path)
    with open(f'{FLAGS.config_folder}/plasmid.genome_set' ,'a') as config_plasmid:
        config_plasmid.write(f'{decoy_name}\n')
    assemblyLengthWriter = open(f'{FLAGS.assembly_dir}/assembly_length' , 'a')
    assemblyPathWriter = open(f'{FLAGS.assembly_dir}/assembly_path' , 'a')
    assemblyTaxidWriter = open(f'{FLAGS.assembly_dir}/assembly_tax_id' , 'a')
    sequenceSummaryWriter = open(f'{FLAGS.assembly_dir}/sequence_summary' , 'a')
    totalLength=0
    with open(path, 'rt') as fi:
        for record in SeqIO.parse(fi, 'fasta'):
            totalLength += len(record)
            sequenceSummaryWriter.write(f"{record.id}\t{len(record)}\t{decoy_name}\n")
    assemblyLengthWriter.write(f"{decoy_name}\t{totalLength}\n")
    assemblyPathWriter.write(f"{decoy_name}\t{path}\n")
    arbitrary_taxid='35'
    assemblyTaxidWriter.write(f"{decoy_name}\t1000000099\t1000000099\t1000000001\t{arbitrary_taxid}\n")
