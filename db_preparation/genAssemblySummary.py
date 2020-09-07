import argparse
import os  
import pandas

FLAGS = None

def summary(refseqGenome, assemblySummary):

    assembly_summary = pandas.read_csv(refseqGenome+'/archaea/assembly_summary.txt', dtype=str, sep='\t', header=1)
    assembly_summary_2 = pandas.read_csv(refseqGenome+'/bacteria/assembly_summary.txt', dtype=str, sep='\t', header=1)
    assembly_summary = assembly_summary.append(assembly_summary_2)
    assembly_summary_2 = pandas.read_csv(refseqGenome+'/vertebrate_mammalian/assembly_summary.txt', dtype=str, sep='\t', header=1)
    assembly_summary = assembly_summary.append(assembly_summary_2)
    assembly_summary_2 = pandas.read_csv(refseqGenome+'/viral/assembly_summary.txt', dtype=str, sep='\t', header=1)
    assembly_summary = assembly_summary.append(assembly_summary_2)
    assembly_summary_2 = pandas.read_csv(refseqGenome+'/fungi/assembly_summary.txt', dtype=str, sep='\t', header=1)
    assembly_summary = assembly_summary.append(assembly_summary_2)
    assembly_summary_2 = pandas.read_csv(refseqGenome+'/protozoa/assembly_summary.txt', dtype=str, sep='\t', header=1)
    assembly_summary = assembly_summary.append(assembly_summary_2)
    assembly_summary_2 = pandas.read_csv(refseqGenome+'/plasmid/assembly_summary.txt', dtype=str, sep='\t')
    assembly_summary = assembly_summary.append(assembly_summary_2)

    assembly_summary.to_csv(assemblySummary+'.original', sep='\t', header=False, index=False)
    
    fileWriter = open(assemblySummary, 'w')
    with open(assemblySummary+'.original', 'r') as f:
        for line in f:
            eles = line.split('\t')
            for i in range(len(eles)):
                eles[i] = eles[i].replace('"', '')
            fileWriter.write("\t".join(eles))
    fileWriter.close()
    os.remove(assemblySummary+'.original')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refseq assembly info update')
    parser.add_argument('--db_dir', default="genomes/refseq/", required=False)
    parser.add_argument('--assemblySummary', default=None, required=False)

    FLAGS, UNPARSED = parser.parse_known_args()

    summary(assemblySummary = FLAGS.assemblySummary, refseqGenome = FLAGS.db_dir)
