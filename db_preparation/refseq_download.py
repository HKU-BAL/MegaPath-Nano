from __future__ import print_function, division

import argparse
import pandas
import os
import sys
import subprocess
import gzip
from Bio import SeqIO

FLAGS = None
class Assembly:
    #assemblyAccession, refseqCategory, taxid, speciesTaxid, assemblyLevel, organismName
    def __init__(self, assemblyAccession, refseqCategory, taxid, speciesTaxid, assemblyLevel, organismName):
        self.assemblyAccession = assemblyAccession
        self.refseqCategory = refseqCategory
        self.taxid = taxid
        self.speciesTaxid = speciesTaxid
        self.assemblyLevel = assemblyLevel
        self.organismName = organismName

    def __str__(self):
        return "%s\t%s\t%s\t%s"%(self.assemblyAccession, self.speciesTaxid, self.taxid, self.organismName)

def getValidAssemble(assembly_summary,assemblySummary_path):
    prevSpeciesTaxid = None
    assemblyDict = {}
    referenceFound = False
    representativeFound = False
    targetAssembly = []
    with open(assemblySummary_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            eles = line.strip().split('\t')
            speciesTaxid = eles[6]
            if prevSpeciesTaxid is None:
                prevSpeciesTaxid = speciesTaxid
            if prevSpeciesTaxid != speciesTaxid:
                for assemblyKey in assemblyDict:
                    assemblyInfo = assemblyDict[assemblyKey]
                    if assemblyInfo.refseqCategory == 'reference genome':
                        targetAssembly.append(assemblyKey)
                    elif assemblyInfo.refseqCategory == 'representative genome':
                        if referenceFound:
                            if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                                continue
                            targetAssembly.append(assemblyKey)
                    elif assemblyInfo.refseqCategory == 'na':
                        if referenceFound or representativeFound:
                            if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                                continue
                            targetAssembly.append(assemblyKey)
                assemblyDict.clear()
                referenceFound = False
                representativeFound = False
                prevSpeciesTaxid = speciesTaxid
            assemblyInfo = Assembly(eles[0], eles[4], eles[5], eles[6], eles[11], eles[7])
            if eles[4] == 'reference genome':
                referenceFound = True
            elif eles[4] == 'representative genome':
                representativeFound = True
            assemblyDict[eles[0]] = assemblyInfo
    if prevSpeciesTaxid != speciesTaxid:
        for assemblyKey in assemblyDict:
            assemblyInfo = assemblyDict[assemblyKey]
            if assemblyInfo.refseqCategory == 'reference genome':
                targetAssembly.append(assemblyKey)
            elif assemblyInfo.refseqCategory == 'representative genome':
               	if referenceFound:
                    if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                        continue
                    targetAssembly.append(assemblyKey)
            elif assemblyInfo.refseqCategory == 'na':
                if referenceFound or representativeFound:
                    if assemblyInfo.assemblyLevel == 'Scaffold' or assemblyInfo.assemblyLevel == 'Contig':
                        continue
                    targetAssembly.append(assemblyKey)
    return targetAssembly

def download(genome, db_dir):
    genome_dir = os.path.join(db_dir, genome)
    if os.path.exists(genome_dir) != True:
        status = os.system('mkdir -p %s'%(genome_dir))
        if status != 0:
            sys.exit()

    if FLAGS.numBatch == 0:
        if FLAGS.get_summary == True:
            print('Download assembly_summary', flush=True)
            if genome == "vertebrate_mammalian":
                status = os.system(
                    'wget -a wget.log -N -P %s ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/%s/Homo_sapiens/assembly_summary.txt'%(genome_dir, genome))
            else:
                status = os.system(
                    'wget -a wget.log -N -P %s ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/%s/assembly_summary.txt'%(genome_dir, genome))
            if status != 0:
                print("assembly_summary failed", flush=True)
                sys.exit()

    assembly_summary = pandas.read_csv('%s/assembly_summary.txt'%(genome_dir), dtype=str, sep='\t',header=1)
    assemblySummary_path = '%s/assembly_summary.txt'%(genome_dir)
    targetAssembly = getValidAssemble(assembly_summary,assemblySummary_path)
    
    for i in range(assembly_summary.shape[0]):

    	if FLAGS.numBatch == 0 or i % FLAGS.numBatch == FLAGS.batch:

            prefix = assembly_summary['ftp_path'][i][assembly_summary['ftp_path'][i].find(assembly_summary['# assembly_accession'][i]):]
            
            print('   ', genome, i + 1, 'of', assembly_summary.shape[0], ': ', assembly_summary['# assembly_accession'][i], end='', flush=True)
            assembly_path=os.path.join(genome_dir, assembly_summary['# assembly_accession'][i])
            if os.path.exists(assembly_path) != True:
                status = os.system('mkdir %s'%(assembly_path))
                if status != 0:
                    print("mkdir failed", flush=True)
                    sys.exit()

            if os.path.exists(os.path.join(assembly_path, "%s_assembly_report.txt"%(prefix))) != True:
                status = os.system("wget -a wget.log -N -P " + assembly_path  + " " +
                                assembly_summary['ftp_path'][i] + "/" + prefix + "_assembly_report.txt")
                if status != 0:
                    print("assembly_report failed", end='', flush=True)
            
            if os.path.exists(os.path.join(assembly_path, "%s_genomic.fna.gz"%(prefix))) != True:
                status = os.system("wget -a wget.log -N -P " + assembly_path + " " +
                                assembly_summary['ftp_path'][i] + "/" + prefix + "_genomic.fna.gz")
                if status != 0:
                    print("fasta failed", end='', flush=True)

            if os.path.exists(os.path.join(assembly_path, "md5checksums.txt")) != True:
                status = os.system("wget -a wget.log -N -P " + assembly_path + " " +
                                assembly_summary['ftp_path'][i] + "/md5checksums.txt")
                if status != 0:
                    print("checksum failed", end='', flush=True)

            print("", flush=True)

def downloadPlasmid(genome, num, db_dir):
    genome_dir = os.path.join(db_dir, genome)
    if os.path.exists(genome_dir) != True:
        status = os.system('mkdir -p %s'%(genome_dir))
        if status != 0:
            sys.exit()

    assemblySummaryWriter = open(os.path.join(genome_dir,"assembly_summary.txt"), 'w')
    assemblySummaryWriter.write("# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material\n")
    for i in range(1, num+1):
        assemblySummaryWriter.write("PLA_00000000%d.1\t\t\t\t\t100000000%d\t100000000%d\tPLA_00000000%d.1\t\t\t\t\t\t\t\t\t\t\t\tftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.%d.1.genomic.fna.gz\t\t\n"%(i, i, i, i, i))
    assemblySummaryWriter.close()

    for i in range(1, num+1):
        assembly_path = os.path.join(genome_dir, "PLA_00000000%d.1"%(i))
        if os.path.exists(assembly_path) != True:
            status = os.system('mkdir -p %s'%(assembly_path))
            if status != 0:
                print("mkdir failed", flush=True)
                sys.exit()

        if os.path.exists(os.path.join(assembly_path, 'plasmid.%d.1.genomic.fna.gz'%(i))) != True:
            status = os.system("wget -a wget.log -N -P " + assembly_path + " ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.%d.1.genomic.fna.gz"%(i))
            if status != 0:
                print("fasta failed", end='', flush=True)
            status = os.system("mv %s %s"%(os.path.join(assembly_path, 'plasmid.%d.1.genomic.fna.gz'%(i)), os.path.join(assembly_path, 'plasmid.%d.1.genomic.original.fna.gz'%(i))))
            if status != 0:
                print("rename fasta failed", end='', flush=True)
            
            path=os.path.join(assembly_path, 'plasmid.%d.1.genomic.original.fna.gz'%(i))
            outputPath = os.path.join(assembly_path, 'plasmid.%d.1.genomic.fna.gz'%(i))
            fileWriter = open(outputPath, 'w')
            with gzip.open(path, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    description = ' '.join(record.description.split(' ')[1:])
                    record.id = record.id+"_PLA"
                    fileWriter.write(">%s %s\n%s\n"%(record.id, description, record.seq))

            fileWriter.close()

def main():
    if FLAGS.viral:
        print("Download viral", flush=True)
        download("viral", FLAGS.db_dir)
    if FLAGS.bacteria:
        print("Download bacteria", flush=True)
        download("bacteria", FLAGS.db_dir)
    if FLAGS.fungi:
        print("Download fungi", flush=True)
        download("fungi", FLAGS.db_dir)
    if FLAGS.protozoa:
        print("Download protozoa", flush=True)
        download("protozoa", FLAGS.db_dir)
    if FLAGS.archaea:
        print("Download archaea", flush=True)
        download("archaea", FLAGS.db_dir)
    if FLAGS.plasmid:
        print("Download plasmid", flush=True)
        downloadPlasmid("plasmid", FLAGS.num, FLAGS.db_dir)
    if FLAGS.vertebrate_mammalian:
        print("Download vertebrate_mammalian", flush=True)
        download("vertebrate_mammalian", FLAGS.db_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refseq download')

    parser.add_argument('--viral', dest='viral', default=False, action='store_true')
    parser.add_argument('--bacteria', dest='bacteria', default=False, action='store_true')
    parser.add_argument('--fungi', dest='fungi', default=False, action='store_true')
    parser.add_argument('--protozoa', dest='protozoa', default=False, action='store_true')
    parser.add_argument('--archaea', dest='archaea', default=False, action='store_true')
    parser.add_argument('--plasmid', dest='plasmid', default=False, action='store_true')
    parser.add_argument('--num', type=int, default=8, required=False)
    parser.add_argument('--vertebrate_mammalian', dest='vertebrate_mammalian', default=False, action='store_true')
    
    parser.add_argument('--get_summary', dest='get_summary', default=False, action='store_true')

    parser.add_argument('--batch', default=0, type=int)
    parser.add_argument('--numBatch', default=0, type=int)
    parser.add_argument('--db_dir', default='./')

    FLAGS, UNPARSED = parser.parse_known_args()

    main()
