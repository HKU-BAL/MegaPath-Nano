#!/usr/bin/python
import sys
import argparse 
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from shutil import copyfile
import psutil
import pysam
import pandas



FLAGS=None
def process_tax_id(bam_path):
    """
    """
    try:
        acc_id_list = subprocess.check_output(
            "samtools view -@ {threads} {bam}|cut -f3|sort|uniq".format(threads=FLAGS.threads,bam=bam_path) ,encoding='utf-8',
            shell=True).strip().split("\n")
    except subprocess.CalledProcessError:
        print("No accession ID is mapped with the taxonomy ID.")
        sys.exit()
    with pysam.FastxFile(FLAGS.REFSEQ_PATH) as fin:
        for entry in fin:
            if str(entry.name) in acc_id_list:
                with open( "ref_{acc_id}.fa".format(acc_id=entry.name), mode='a') as fout, open("ref_{acc_id}.bed".format(acc_id=entry.name), mode='a') as bed_out:
                    fout.write(str(entry)+'\n')
                    bed_out.write("{name}\t0\t{end}\n".format(name=entry.name,end=len(entry.sequence)))
    with ThreadPoolExecutor(int(FLAGS.threads)) as executor:
        for acc_id in acc_id_list:
            executor.submit(process_accession_num,acc_id)

def process_accession_num(acc_id):
    """
    """
    
    copyfile("header.sam","sample_{acc_id}.sam".format(acc_id=acc_id))
    os.system("samtools view sample.sorted.bam | awk '$3 ~ /{acc_id}/' >> sample_{acc_id}.sam;" 
    "samtools view -b sample_{acc_id}.sam > sample_{acc_id}.bam;samtools index sample_{acc_id}.bam;"
    "bedtools bamtobed -i  sample_{acc_id}.bam |bedops -m - > sample_{acc_id}.merged.bed;"
    "bedops -d ref_{acc_id}.bed sample_{acc_id}.merged.bed > sample_{acc_id}.0cov.bed;"
    "bcftools mpileup -R sample_{acc_id}.merged.bed -Ou -f ref_{acc_id}.fa sample_{acc_id}.bam | bcftools call -Oz -mv -o calls_{acc_id}.vcf.gz;"
    "tabix calls_{acc_id}.vcf.gz;"
    "cat ref_{acc_id}.fa | bcftools consensus -m sample_{acc_id}.0cov.bed calls_{acc_id}.vcf.gz > cns_{acc_id}.fa".format(acc_id=acc_id))
    os.makedirs("results/{acc_id}".format(acc_id=acc_id), exist_ok=True)

    reformat_command = "gawk -F'\\t' 'NR>1 {if($1 in g){if($2>s[$8][$1]){s[$8][$1]=$2}}else {g[$8][$1]=$1; s[$8][$1]=$2};}END {for (i in g) {for (j in g[i]){fg[i]=sprintf(\"%s;%s\", fg[i], g[i][j]);fs[i]=sprintf(\"%s;%s\", fs[i],s[i][j]);}}for (i in g){if(i==\"\"){print \"Unknown \\t \" fg[i] \"\\t\" fs[i]}else {print i \"\\t \" fg[i] \"\\t\" fs[i]}}}'"

    os.makedirs("results/{acc_id}/resf".format(acc_id=acc_id), exist_ok=True)
    r_subprocess = subprocess.Popen('python resfinder/resfinder.py -p amr_db/resfinder/ -i cns_{acc_id}.fa -o results/{acc_id}/resf -t 0.90 -l 0.60  && '.format(acc_id=acc_id) + \
    reformat_command + ' results/{acc_id}/resf/results_tab.txt > results/{acc_id}/resf/resf_temp.txt && '.format(acc_id=acc_id) + \
    "awk -F'\\t' -vOFS=, '{gsub(\" resistance\", \"\");split($1,a,\", \"); for (i in a) print a[i]\"\\t\"$2\"\\t\"$3}'  results/%s/resf/resf_temp.txt >  results/%s/resf/resf_temp2.txt" % (acc_id, acc_id)+ ' && ' + \
    "awk -F'\\t' -vOFS=, '{split($1,a,\"and \"); for (i in a) if(a[i]==\"\") {} else { print a[i]\"\\t\"$2\"\\t\"$3} ;}'  results/%s/resf/resf_temp2.txt >  results/%s/resf/resf.txt" % (acc_id, acc_id), shell=True)



    os.makedirs("results/{acc_id}/card".format(acc_id=acc_id), exist_ok=True)
    c_subprocess = subprocess.Popen('rgi main --input_sequence cns_{acc_id}.fa --output_file results/{acc_id}/card/results --input_type contig -n {threads} && '
        "awk -F'\t' 'NR>1{{print $15 \"\\t\" $9 \"\\t\" $10}}' results/{acc_id}/card/results.txt > results/{acc_id}/card/card_temp.txt && "
        "awk -F'\t' -vOFS=, '{{split($1,a,\"; \"); for (i in a) if(a[i] !~ /antibiotic/){{print a[i]\"\\t\"$2\"\\t\"$3}} else {{split(a[i],b,\" \"); print b[1]\"\\t\"$2\"\\t\"$3}} ;}}' results/{acc_id}/card/card_temp.txt > results/{acc_id}/card/card.txt".format(acc_id=acc_id,threads=FLAGS.threads), shell=True)


    os.makedirs("results/{acc_id}/amrf".format(acc_id=acc_id), exist_ok=True )
    add_flag=''
    if FLAGS.taxon != None:
        add_flag='-O '+FLAGS.taxon
    a_subprocess = subprocess.Popen('amrfinder {add_flag} -n cns_{acc_id}.fa > results/{acc_id}/amrf/results.txt && ' 
                "awk -F'\\t' 'NR>1{{print $12 \"\\t\" $6 \"\\t\" $17}}' results/{acc_id}/amrf/results.txt > results/{acc_id}/amrf/amrf_temp.txt && " 
                "awk -F'\\t' -vOFS=, '{{split($1,a,\"; \"); for (i in a) if(a[i] !~ /antibiotic/){{print a[i]\"\\t\"$2\"\\t\"$3}} else {{split(a[i],b,\" \"); print b[1]\"\\t\"$2\"\\t\"$3}} ;}}' results/{acc_id}/amrf/amrf_temp.txt > results/{acc_id}/amrf/amrf.txt".format(add_flag=add_flag,acc_id=acc_id), shell=True)
    os.makedirs("results/{acc_id}/megares".format(acc_id=acc_id), exist_ok=True )
    m_subprocess = subprocess.Popen('python {bin_dir}/blast_amr.py -i cns_{acc_id}.fa -o results/{acc_id}/megares/ -d megares_full_database_v2.00 -p {bin_dir}/amr_db/megares -l 0.9 -t 0.6  && ' 
            '{reformat_command} results/{acc_id}/megares/results_tab.txt > results/{acc_id}/megares/megares.txt'.format(acc_id=acc_id,reformat_command=reformat_command,bin_dir=FLAGS.NANO_DIR_PATH+"/bin"), shell=True)

    os.makedirs("results/{acc_id}/cbmar".format(acc_id=acc_id), exist_ok=True )
    cbmar_nucl_subprocess = subprocess.Popen('python {bin_dir}/blast_amr.py -i cns_{acc_id}.fa -o results/{acc_id}/cbmar/ -d cbmar_nucl -p {bin_dir}/amr_db/cbmar -l 0.9 -t 0.6  && ' 
            '{reformat_command} results/{acc_id}/cbmar/results_tab.txt > results/{acc_id}/cbmar/cbmar.txt'.format(acc_id=acc_id,reformat_command=reformat_command,bin_dir=FLAGS.NANO_DIR_PATH+"/bin"), shell=True)
    cbmar_prot_subprocess = subprocess.Popen('prodigal -m -a cns_{acc_id}.prot.fa -i cns_{acc_id}.fa && '
            'blastp -query cns_{acc_id}.prot.fa -out cns_{acc_id}.prot.fa.blast.txt -db {CBMAR_PROT_DB_PATH} -outfmt "6 qseqid sseqid pident qcovs"'.format(acc_id=acc_id,CBMAR_PROT_DB_PATH=FLAGS.CBMAR_PROT_DB_PATH), shell=True)

    r_subprocess.communicate()
    c_subprocess.communicate()
    a_subprocess.communicate()
    m_subprocess.communicate()
    cbmar_nucl_subprocess.communicate()
    #  TODO separate protein subprocess
    cbmar_prot_subprocess.communicate()


def parse_output(acc_id,db_acc_id,db_gene,db_score,db_name):
    if os.path.isfile("{acc_id}/{db_name}/{db_name}.txt".format(acc_id=acc_id,db_name=db_name)) and os.path.getsize("{acc_id}/{db_name}/{db_name}.txt".format(acc_id=acc_id,db_name=db_name)):
        with open("{acc_id}/{db_name}/{db_name}.txt".format(acc_id=acc_id,db_name=db_name)) as f:
            for line in f:
                col = line.split("\t")
                drug=col[0].lower().strip()
                if drug.endswith('s'):
                    drug = drug[:-1]
                if drug not in db_acc_id:
                    db_acc_id[drug] = acc_id
                    db_gene[drug] = col[1]
                    db_score[drug] = col[2].strip()
                else:
                    if acc_id not in db_acc_id[drug]:
                        db_acc_id[drug] = db_acc_id[drug] + "|" + acc_id
                        db_gene[drug] = db_gene[drug] + "|"
                        db_score[drug] = db_score[drug] + "|"

                    if col[1].strip() != "":
                        if db_gene[drug][-1] != "|":
                            db_gene[drug] = db_gene[drug] + ";" + col[1]
                            db_score[drug] = db_score[drug] + ";" + col[2].strip()
                        elif db_gene[drug][-1] == "|":
                            db_gene[drug] = db_gene[drug] + col[1]
                            db_score[drug] = db_score[drug] + col[2].strip()

def merge_results(dir_arr):
    #TODO threading
    """
    consolidating the results of all accession IDs
    :param dir_arr:
    :return:
    """
    card_acc_id = {}
    card_gene = {}
    card_score = {}

    resf_acc_id = {}
    resf_gene = {}
    resf_score = {}

    amrf_acc_id = {}
    amrf_gene = {}
    amrf_score = {}

    mega_acc_id = {}
    mega_gene = {}
    mega_score = {}

    cbmar_acc_id = {}
    cbmar_gene = {}
    cbmar_score = {}

    for acc_id in dir_arr:
        for db_data in ((card_acc_id,card_gene,card_score,'card'),(resf_acc_id,resf_gene,resf_score,'resf'),(amrf_acc_id,amrf_gene,amrf_score,'amrf'),(mega_acc_id,mega_gene,mega_score,'mega'),(cbmar_acc_id,cbmar_gene,cbmar_score,'cbmar')):
            parse_output(acc_id,*db_data)
    
    df_dict=dict()
    for antibiotic in set().union(card_acc_id,resf_acc_id,mega_acc_id,cbmar_acc_id,amrf_acc_id):
        if antibiotic not in df_dict:
            df_dict[antibiotic]=dict()
        if antibiotic in card_acc_id:
            df_dict[antibiotic]['card_acc_id']=card_acc_id[antibiotic]
            df_dict[antibiotic]['card_gene']=card_gene[antibiotic]
            df_dict[antibiotic]['card_score']=card_score[antibiotic]
        if antibiotic in resf_acc_id:
            df_dict[antibiotic]['resfinder_acc_id']=resf_acc_id[antibiotic]
            df_dict[antibiotic]['resfinder_gene']=resf_gene[antibiotic]
            df_dict[antibiotic]['resfinder_score']=resf_score[antibiotic]
        if antibiotic in mega_acc_id:
            df_dict[antibiotic]['megares_acc_id']=mega_acc_id[antibiotic]
            df_dict[antibiotic]['megares_gene']=mega_gene[antibiotic]
            df_dict[antibiotic]['megares_score']=mega_score[antibiotic]
        if antibiotic in amrf_acc_id:
            df_dict[antibiotic]['amrfinder_acc_id']=amrf_acc_id[antibiotic]
            df_dict[antibiotic]['amrfinder_gene']=amrf_gene[antibiotic]
            df_dict[antibiotic]['amrfinder_score']=amrf_score[antibiotic]
        if antibiotic in cbmar_acc_id:
            df_dict[antibiotic]['cbmar_acc_id']=cbmar_acc_id[antibiotic]
            df_dict[antibiotic]['cbmar_gene']=cbmar_gene[antibiotic]
            df_dict[antibiotic]['cbmar_score']=cbmar_score[antibiotic]

    df=pandas.DataFrame(df_dict).T.sort_index(axis=0).reindex(['card_acc_id','card_gene','card_score','amrfinder_acc_id','amrfinder_gene','amrfinder_score','resfinder_acc_id','resfinder_gene','resfinder_score','cbmar_acc_id','cbmar_gene','cbmar_score','megares_acc_id','megares_gene','megares_score'],axis=1)
    df.to_csv('results.txt')

    #  output cbmar protein results
    family_details_df=pandas.read_csv(FLAGS.FAMILY_DETAILS_PATH,encoding = "utf-8")
    drug_set=set()
    for f in glob.glob('%s/cns*fa.prot.fa.blast.txt' %FLAGS.output_folder):
        if os.path.getsize(f)>0:
            df = pandas.read_csv(f,sep='\t',header=None)
            df[2]=df[2].astype(int)
            df[3]=df[3].astype(int)
            df=df[(df[2]>=FLAGS.pident_cbmar) & (df[3]>=FLAGS.qcovs_cbmar)]
            for index, row in family_details_df.iterrows():
                if df[1].str.contains(row['Uniprot ID']).any():
                    drug_set.add(row['Hydrolytic profile'])
    with open('cbmar_protein_blasted_hydrolytic_profile.txt','w') as fo:
        for drug in drug_set:
            fo.write(str(drug)+'\n')


def main():
    os.makedirs(FLAGS.output_folder,exist_ok=True)
    bam_path = os.path.abspath(FLAGS.query_bam)
    os.chdir(FLAGS.output_folder)
    print("current directory: " + os.getcwd())

    if not os.path.exists("sample.sorted.bam") or not os.path.exists("sample.sorted.bam.bai") or not os.path.exists("header.sam") :
        os.system("samtools sort {bam_path} -o sample.sorted.bam;" 
                "samtools index sample.sorted.bam;"
                "samtools view -@ {threads} -H sample.sorted.bam > header.sam".format(bam_path=bam_path,threads=FLAGS.threads))

    os.makedirs("results",exist_ok=True)
    
    process_tax_id(bam_path)

    print("All results have been generated")
    os.chdir("results")
    print("current directory " + os.getcwd())

    try:
        dir_arr = next(os.walk('.'), ([],[],[]))[1]
    except subprocess.CalledProcessError:
        print("No accession ID is mapped with the sequence ID in bam.")
    else:
        print("Merging results")
        merge_results(dir_arr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MegaPath-Nano: AMR Detection') 
    parser.add_argument('--query_bam', required=True,help='Input bam')
    parser.add_argument('--output_folder', required=True,help='Output directory')
    parser.add_argument('--taxon', help='Taxon-specific options for AMRFinder, curated organisms: Campylobacter, Enterococcus_faecalis, Enterococcus_faecium, Escherichia, Klebsiella, Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius, Vibrio_cholerae')
    parser.add_argument('--threads', default=int(psutil.cpu_count(logical=True)/2), help='Num of threads')
    parser.add_argument('--pident_cbmar', default=90, help='The threshold of percentage of identical matches in blastp')
    parser.add_argument('--qcovs_cbmar', default=60, help='The threshold of percentage of query coverage in blastp')
    CWD=os.path.dirname(os.path.realpath(__file__))
    NANO_DIR=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--REFSEQ_PATH', default='%s/genomes/refseq/refseq.fna'%(NANO_DIR), help='The path of reference files. RefSeq by default')
    parser.add_argument('--NANO_DIR_PATH', default=NANO_DIR, help='The path of root directory of MegaPath-Nano')
    parser.add_argument('--CBMAR_PROT_DB_PATH', default='%s/bin/amr_db/cbmar/cbmar_prot.fsa'%(NANO_DIR), help='The path of betalactamase family details in protein, collected from http://proteininformatics.org/mkumar/lactamasedb/cllasification.html.')
    parser.add_argument('--FAMILY_DETAILS_PATH', default='%s/bin/amr_db/cbmar/family_details.csv'%(NANO_DIR), help='The path of betalactamase family details in protein, collected from http://proteininformatics.org/mkumar/lactamasedb/cllasification.html.')
    FLAGS = parser.parse_args()
    main()
