#!/usr/bin/env python3
import sys
import argparse 
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from shutil import copyfile
import psutil
import pysam
import pandas
import numpy



def process_tax_id(bam_path,FLAGS):
    """
    """
    acc_id_list=set()
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        for read in samfile:
            acc_id_list.add(read.reference_name)

    if acc_id_list==False:
        print("No accession ID is mapped.")
        sys.exit()

    with pysam.FastxFile(FLAGS.REFSEQ_PATH) as fin, open(f"ref.fa", mode='w') as fout, open(f"ref.bed", mode='w') as bed_out:
        for entry in fin:
            if entry.name in acc_id_list:
                fout.write(f'{entry}\n')
                bed_out.write(f'{entry.name}\t0\t{len(entry.sequence)}\n')
    process_accession_num(bam_path,FLAGS)



def process_accession_num(bam_path,FLAGS):
    """
    """
    bam_basename=os.path.basename(bam_path)
    consensus_subprocess=subprocess.Popen(
            f"bedtools bamtobed -i {bam_path} |bedops -m - > {bam_basename}.merged.bed;"
            f"bedops -d ref.bed {bam_basename}.merged.bed > {bam_basename}.0cov.bed;"
            f"bcftools mpileup --threads {FLAGS.threads} -R {bam_basename}.merged.bed -Ou -f ref.fa {bam_path} | bcftools call --threads {FLAGS.threads} -Oz -mv -o calls.vcf.gz;"
            "tabix calls.vcf.gz;"
            f"cat ref.fa | bcftools consensus -m {bam_basename}.0cov.bed calls.vcf.gz |seqtk cutN - > cns.fa",shell=True)
    consensus_subprocess.communicate()

    os.makedirs("results/resf", exist_ok=True)
    r_subprocess = subprocess.Popen(f'python {FLAGS.NANO_DIR_PATH}/bin/resfinder/resfinder.py -p {FLAGS.NANO_DIR_PATH}/bin/amr_db/resfinder/ -i cns.fa -o results/resf -x',shell=True)
    os.makedirs("results/card", exist_ok=True)
    c_subprocess = subprocess.Popen(f'rgi main --input_sequence cns.fa --output_file results/card/results_tab.tsv --input_type contig -n {FLAGS.threads}',shell=True)
    os.makedirs("results/amrf", exist_ok=True )
    add_flag=''
    if FLAGS.taxon != None:
        add_flag='-O '+FLAGS.taxon
    a_subprocess = subprocess.Popen(f'amrfinder {add_flag} -n cns.fa > results/amrf/results_tab.tsv',shell=True)
    os.makedirs("results/megares", exist_ok=True )
    m_subprocess = subprocess.Popen(f'blastn -subject {FLAGS.NANO_DIR_PATH}/bin/amr_db/megares/megares_full_database_v2.00.fsa -query cns.fa -out results/megares/results_tab.tsv -outfmt "6 qseqid sseqid pident qcovs" -qcov_hsp_perc {FLAGS.blast_qcov_hsp_perc} -perc_identity {FLAGS.blast_perc_identity}',shell=True)
    os.makedirs("results/cbmar", exist_ok=True )
    cbmar_nucl_subprocess = subprocess.Popen(f'blastn -subject {FLAGS.NANO_DIR_PATH}/bin/amr_db/cbmar/cbmar_nucl.fsa -query cns.fa -out results/cbmar/results_tab.tsv -outfmt "6 qseqid sseqid pident qcovs" -qcov_hsp_perc {FLAGS.blast_qcov_hsp_perc} -perc_identity {FLAGS.blast_perc_identity}',shell=True)
    cbmar_prot_subprocess = subprocess.Popen('prodigal -m -a cns.prot.fa -i cns.fa && '
            f'blastp -query cns.prot.fa -out cns.prot.fa.blast.txt -db {FLAGS.CBMAR_PROT_DB_PATH} -outfmt "6 qseqid sseqid pident qcovs"', shell=True)

    r_subprocess.communicate()
    c_subprocess.communicate()
    a_subprocess.communicate()
    m_subprocess.communicate()
    cbmar_nucl_subprocess.communicate()
    #  TODO separate protein subprocess
    cbmar_prot_subprocess.communicate()

def format_blast_output(in_path):
    out_path=in_path.replace('_blast.txt','.txt')
    with open(in_path) as fi, open(out_path,'w') as fo:
        if 'megares' in in_path:
            drug_class=["Aminocoumarins", "Aminoglycosides", "Bacitracin", "betalactams", "Cationic_antimicrobial_peptides", "Elfamycins", "Fluoroquinolones", "Fosfomycin", "Fusidic_acid", "Glycopeptides", "Lipopeptides", "Metronidazole", "MLS", "Multi-drug_resistance", "Mycobacterium_tuberculosis-specific_Drug", "Phenicol", "Rifampin", "Sulfonamides", "Tetracyclines", "Thiostrepton", "Trimethoprim", "Tunicamycin"]
            drug_dict=dict()
            for line in fi:
                header,ID=line.split('\t')[1:3]
                for drug in drug_class:
                    if drug in header:
                        if drug in drug_dict:
                            drug_dict[drug]=[f'{drug_dict[drug][0]};{header}',f'{drug_dict[drug][1]};{ID}']
                        else:
                            drug_dict[drug]=[header,ID]
                        break
            for drug, header_id_list in drug_dict.items():
                fo.write(f'{k}\t{header_id_list[0]}\t{header_id_list}[1]\n')
        elif 'cbmar' in in_path:
            drug="Betalactamase"
            header_str=''
            ID_str=''
            for line in fi:
                header,ID=line.split('\t')[1:3]
                header_str=f'{header_str};{header}'
                ID_str=f'{ID_str};{ID}'
            fo.write(f'{drug}\t{header_str}\t{ID_str}\n')

def canonicalize(name):
    if name.endswith('s'):
        name=name[:-1]
    return name.replace('-','').upper()

def parse_output(db_acc_id,db_gene_idscore,db_path):
    if os.path.isfile(f"results/{db_path}") and os.path.getsize(f"results/{db_path}"):
        with open(f"results/{db_path}") as f:
            if any(db in db_path for db in ('card','resf','amrf')):
                next(f)
            for line in f:
                if 'card' in db_path:
                    col=line.split('\t')
                    gene_idscore=f'{col[8]}[{col[9]}]'
                    acc_id=col[1]
                    drug=col[14]
                    #acc_id=acc_id.split(':')[0]
                elif 'resf' in db_path:
                    col=line.split('\t')
                    gene_idscore=f'{col[1]}[{col[2]}]'
                    acc_id=col[4]
                    drug=col[0]
                elif 'amrf' in db_path:
                    col=line.split('\t')
                    gene_idscore=f'{col[5]}[{col[16]}]'
                    acc_id=col[1]
                    drug=col[11]
                    #acc_id=acc_id.split(':')[0]
                elif 'megares' in db_path:
                    #put arg
                    drug_class=["AMINOCOUMARINS", "AMINOGLYCOSIDES", "BACITRACIN", "BETALACTAMS", "CATIONIC_ANTIMICROBIAL_PEPTIDES", "ELFAMYCINS", "FLUOROQUINOLONES", "FOSFOMYCIN", "FUSIDIC_ACID", "GLYCOPEPTIDES", "LIPOPEPTIDES", "METRONIDAZOLE", "MLS", "MULTI-DRUG_RESISTANCE", "MYCOBACTERIUM_TUBERCULOSIS-SPECIFIC_DRUG", "PHENICOL", "RIFAMPIN", "SULFONAMIDES", "TETRACYCLINES", "THIOSTREPTON", "TRIMETHOPRIM", "TUNICAMYCIN"]
                    drug_dict=dict()
                    acc_id,gene,id_score=line.split('\t')[0:3]
                    gene_idscore=f'{gene}[{id_score}]'
                    #acc_id=acc_id.split(':')[0]
                    for pheno in drug_class:
                        if pheno in gene_idscore:
                            drug=pheno
                            break
                elif 'cbmar' in db_path:
                    drug="BETALACTAMASE"
                    acc_id,gene,id_score=line.split('\t')[0:3]
                    gene_idscore=f'{gene}[{id_score}]'
                    #acc_id=acc_id.split(':')[0]
                drug=canonicalize(drug)
                if drug not in db_acc_id:
                    db_acc_id[drug] = acc_id
                    db_gene_idscore[drug] = gene_idscore.strip()
                elif drug in db_acc_id:
                    if acc_id not in db_acc_id[drug]:
                        db_acc_id[drug] = f'{db_acc_id[drug]}:{acc_id}'
                        db_gene_idscore[drug] = f'{db_gene_idscore[drug]}:'

                    if gene_idscore.strip() != "":
                        if db_gene_idscore[drug][-1] != ":":
                            db_gene_idscore[drug] = f'{db_gene_idscore[drug]};{gene_idscore}'
                        elif db_gene_idscore[drug][-1] == ":":
                            db_gene_idscore[drug] = f'{db_gene_idscore[drug]}{gene_idscore}'

def merge_results():
    #TODO threading
    """
    consolidating the results of all accession IDs
    :return:
    """
    card_acc_id = {}
    card_gene_idscore = {}

    resf_acc_id = {}
    resf_gene_idscore = {}

    amrf_acc_id = {}
    amrf_gene_idscore = {}

    mega_acc_id = {}
    mega_gene_idscore = {}

    cbmar_acc_id = {}
    cbmar_gene_idscore = {}

    for db_data in ((card_acc_id,card_gene_idscore,'card/results_tab.tsv.txt'),(resf_acc_id,resf_gene_idscore,'resf/results_tab.tsv'),(amrf_acc_id,amrf_gene_idscore,'amrf/results_tab.tsv'),(mega_acc_id,mega_gene_idscore,'megares/results_tab.tsv'),(cbmar_acc_id,cbmar_gene_idscore,'cbmar/results_tab.tsv')):
        parse_output(*db_data)

    df_dict=dict()
    for antibiotic in set().union(card_acc_id,resf_acc_id,mega_acc_id,cbmar_acc_id,amrf_acc_id):
        if antibiotic not in df_dict:
            df_dict[antibiotic]=dict()
        if antibiotic in card_acc_id:
            df_dict[antibiotic]['card_acc_id']=card_acc_id[antibiotic]
            df_dict[antibiotic]['card_gene[idscore]']=card_gene_idscore[antibiotic]
        if antibiotic in resf_acc_id:
            df_dict[antibiotic]['resfinder_acc_id']=resf_acc_id[antibiotic]
            df_dict[antibiotic]['resfinder_gene[idscore]']=resf_gene_idscore[antibiotic]
        if antibiotic in mega_acc_id:
            df_dict[antibiotic]['megares_acc_id']=mega_acc_id[antibiotic]
            df_dict[antibiotic]['megares_gene[idscore]']=mega_gene_idscore[antibiotic]
        if antibiotic in amrf_acc_id:
            df_dict[antibiotic]['amrfinder_acc_id']=amrf_acc_id[antibiotic]
            df_dict[antibiotic]['amrfinder_gene[idscore]']=amrf_gene_idscore[antibiotic]
        if antibiotic in cbmar_acc_id:
            df_dict[antibiotic]['cbmar_acc_id']=cbmar_acc_id[antibiotic]
            df_dict[antibiotic]['cbmar_gene[idscore]']=cbmar_gene_idscore[antibiotic]

    df=pandas.DataFrame(df_dict).T
    df.to_csv('results.csv')

    #  output cbmar protein results
    family_details_df=pandas.read_csv(FLAGS.FAMILY_DETAILS_PATH,encoding = "utf-8")
    drug_set=set()
    df = pandas.read_csv(f'cns.prot.fa.blast.txt',sep='\t',header=None)
    df[2]=df[2].astype(int)
    df[3]=df[3].astype(int)
    df=df[(df[2]>=FLAGS.blast_perc_identity) & (df[3]>=FLAGS.blast_qcov_hsp_perc)]
    for index, row in family_details_df.iterrows():
        if type(row['Uniprot ID'])==str:
            if df[1].str.contains(row['Uniprot ID']).any():
                drug_set.add(row['Hydrolytic profile'])
    with open('cbmar_protein_blasted_hydrolytic_profile.txt','w') as fo:
        for drug in drug_set:
            fo.write(str(drug)+'\n')


def main(FLAGS):
    bam_path = os.path.abspath(FLAGS.query_bam)
    os.makedirs(FLAGS.output_folder,exist_ok=True)
    os.chdir(FLAGS.output_folder)
    os.makedirs("results",exist_ok=True)

    process_tax_id(bam_path,FLAGS)

    print("All results have been generated")

    print("Merging results")
    merge_results()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MegaPath-Nano: AMR Detection') 
    parser.add_argument('--query_bam', required=True,help='Input sorted and indexed bam')
    parser.add_argument('--output_folder', required=True,help='Output directory')
    parser.add_argument('--taxon', help='Taxon-specific options for AMRFinder, curated organisms: Campylobacter, Enterococcus_faecalis, Enterococcus_faecium, Escherichia, Klebsiella, Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius, Vibrio_cholerae')
    parser.add_argument('--threads', default=int(psutil.cpu_count(logical=True)/2), help='Num of threads')
    parser.add_argument('--blast_perc_identity', default=90, help='The threshold of percentage of identical matches in blast')
    parser.add_argument('--blast_qcov_hsp_perc', default=60, help='The threshold of percentage of query coverage in blast')
    CWD=os.path.dirname(os.path.realpath(__file__))
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--REFSEQ_PATH', default=f'{NANO_DIR}/genomes/refseq/refseq.fna.gz', help='The path of reference files. RefSeq by default')
    parser.add_argument('--NANO_DIR_PATH', default=NANO_DIR, help='The path of root directory of MegaPath-Nano')
    parser.add_argument('--CBMAR_PROT_DB_PATH', default=f'{NANO_DIR}/bin/amr_db/cbmar/cbmar_prot.fsa', help='The path of betalactamase family details in protein, collected from http://proteininformatics.org/mkumar/lactamasedb/cllasification.html.')
    parser.add_argument('--FAMILY_DETAILS_PATH', default=f'{NANO_DIR}/bin/amr_db/cbmar/cbmar_family_details.csv', help='The path of betalactamase family details in protein, collected from http://proteininformatics.org/mkumar/lactamasedb/cllasification.html.')
    FLAGS = parser.parse_args()
    main(FLAGS)
