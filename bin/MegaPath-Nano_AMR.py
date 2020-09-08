#!/usr/bin/python
import psutil
import sys, getopt, subprocess, os, time, resource
import pysam
from concurrent.futures import ThreadPoolExecutor
import argparse 

FLAGS=None
def processTaxID(bam_path):
    """
    """
    try:
        acc_id_list = subprocess.check_output(
            "samtools view -@ {threads} {bam}|cut -f3|sort|uniq".format(threads=FLAGS.threads,bam=bam_path) ,encoding='utf-8',
            shell=True).strip().split("\n")
    except subprocess.CalledProcessError:
        print("No accession ID is mapped with the taxonomy ID.")
        exit()
    with pysam.FastxFile(FLAGS.REFSEQ_PATH) as fin:
        for entry in fin:
            if str(entry.name) in acc_id_list:
                with open( "ref_{acc_id}.fa".format(acc_id=entry.name), mode='a') as fout, open("ref_{acc_id}.bed".format(acc_id=entry.name), mode='a') as bed_out:
                    fout.write(str(entry)+'\n')
                    bed_out.write("{name}\t0\t{end}\n".format(name=entry.name,end=len(entry.sequence)))
    with ThreadPoolExecutor(int(FLAGS.threads)) as executor:
        for acc_id in acc_id_list:
            executor.submit(processAccessionNo,acc_id)
    pass

def processAccessionNo(acc_id):
    """
    """
    curr_time = time.asctime(time.localtime(time.time()))
    print("Processing accession ID - " + acc_id + " at " + curr_time)
    
    os.system("cp header.sam sample_{acc_id}.sam".format(acc_id=acc_id))
    os.system("samtools view sample.sorted.bam | awk '$3 ~ /{acc_id}/' >> sample_{acc_id}.sam".format(acc_id=acc_id))
    os.system("samtools view -@ {threads} -b sample_{acc_id}.sam > sample_{acc_id}.bam".format(threads=FLAGS.threads,acc_id=acc_id))
    os.system("bedtools bamtobed -i  sample_{acc_id}.bam |bedops -m - > sample_{acc_id}.merged.bed".format(acc_id=acc_id))
    os.system("bedops -d ref_{acc_id}.bed sample_{acc_id}.merged.bed > sample_{acc_id}.0cov.bed".format(acc_id=acc_id))
    os.system("bcftools mpileup -R sample_{acc_id}.merged.bed -Ou -f ref_{acc_id}.fa sample_{acc_id}.bam | bcftools call -Oz -mv --threads {threads} -o calls_{acc_id}.vcf.gz".format(acc_id=acc_id,threads=FLAGS.threads))
    os.system("tabix calls_{acc_id}.vcf.gz".format(acc_id=acc_id))
    os.system("cat ref_{acc_id}.fa | bcftools consensus -m sample_{acc_id}.0cov.bed calls_{acc_id}.vcf.gz > cns_{acc_id}.fa".format(acc_id=acc_id))
    os.makedirs("results/{acc_id}".format(acc_id=acc_id), exist_ok=True)

    reformat_command = "gawk -F'\\t' 'NR>1 {if($1 in g){if($2>s[$8][$1]){s[$8][$1]=$2}}else {g[$8][$1]=$1; s[$8][$1]=$2};}END {for (i in g) {for (j in g[i]){fg[i]=sprintf(\"%s;%s\", fg[i], g[i][j]);fs[i]=sprintf(\"%s;%s\", fs[i],s[i][j]);}}for (i in g){if(i==\"\"){print \"Unknown \\t \" fg[i] \"\\t\" fs[i]}else {print i \"\\t \" fg[i] \"\\t\" fs[i]}}}'"

    os.makedirs("results/{acc_id}/resfinder".format(acc_id=acc_id), exist_ok=True)
    r_subprocess = subprocess.Popen("python resfinder/resfinder.py -p amr_db/resfinder/ -i cns_{acc_id}.fa -o results/{acc_id}/resfinder -t 0.90 -l 0.60 ".format(acc_id=acc_id) + ' && ' + \
            reformat_command + " results/{acc_id}/resfinder/results_tab.txt > results/{acc_id}/resfinder/resf_temp.txt".format(acc_id=acc_id)+ ' && ' + \
            #  separate ", "
             "awk -F'\\t' -vOFS=, '{split($1,a,\", \"); for (i in a) print a[i]\"\\t\"$2\"\\t\"$3}'  results/{acc_id}/resfinder/resf_temp.txt >  results/{acc_id}/resfinder/resf_temp2.txt".format(acc_id=acc_id)+ ' && ' + \
            #  separate "and "
             "awk -F'\\t' -vOFS=, '{split($1,a,\"and \"); for (i in a) if(a[i]==\"\") {} else { print a[i]\"\\t\"$2\"\\t\"$3} ;}'  results/{acc_id}/resfinder/resf_temp2.txt >  results/{acc_id}/resfinder/resf.txt".format(acc_id=acc_id), shell=True)

    os.makedirs("results/{acc_id}/card".format(acc_id=acc_id), exist_ok=True)
    c_subprocess = subprocess.Popen("rgi main --input_sequence cns_{acc_id}.fa --output_file results/{acc_id}/card/results --input_type contig -n {threads}".format(acc_id=acc_id,threads=FLAGS.threads) + ' && ' + \
        "awk -F'\t' 'NR>1{print $15 \"|\" $9 \"|\" $10}' results/{acc_id}/card/results.txt > results/{acc_id}/card/card_temp.txt".format(acc_id=acc_id)  + ' && ' + \
        "awk -F'|' -vOFS=, '{split($1,a,\"; \"); for (i in a) if(a[i] !~ /antibiotic/){print a[i]\"|\"$2\"|\"$3} else {split(a[i],b,\" \"); print b[1]\"|\"$2\"|\"$3} ;}' results/{acc_id}/card/card_temp.txt > results/{acc_id}/card/card.txt".format(acc_id=acc_id), shell=True)


    os.makedirs("results/{acc_id}/amrfinder".format(acc_id=acc_id), exist_ok=True )
    add_flag=''
    if FLAGS.taxon != '':
        add_flag='-O '+FLAGS.taxon
    a_subprocess = subprocess.Popen("amrfinder {add_flag} -n cns_{acc_id}.fa > results/{acc_id}/amrfinder/results.txt".format(acc_id=acc_id,add_flag=add_flag)+ ' && ' + \
                "awk -F'\\t' 'NR>1{print $12 \"\\t\" $6 \"\\t\" $17}' results/{acc_id}/amrfinder/results.txt > results/{acc_id}/amrfinder/amrf_temp.txt".format(acc_id=acc_id)+ ' && ' + \
                "awk -F'\\t' -vOFS=, '{split($1,a,\"; \"); for (i in a) if(a[i] !~ /antibiotic/){print a[i]\"\\t\"$2\"\\t\"$3} else {split(a[i],b,\" \"); print b[1]\"\\t\"$2\"\\t\"$3} ;}' results/{acc_id}/amrfinder/amrf_temp.txt > results/{acc_id}/amrfinder/amrf.txt".format(acc_id=acc_id), shell=True)
    os.makedirs("results/{acc_id}/megares".format(acc_id=acc_id), exist_ok=True )
    m_subprocess = subprocess.Popen("python blast_amr.py -i cns_{acc_id}.fa -o results/{acc_id}/megares/ -d megares_database_v1.01_AMRDetection -p amr_db/megares_db -l 0.9 -t 0.6 ".format(acc_id=acc_id)+ ' && ' + \
            reformat_command+" results/{acc_id}/megares/results_tab.txt > results/{acc_id}/megares/megares.txt".format(acc_id=acc_id), shell=True)

    os.makedirs("results/{acc_id}/cbmar".format(acc_id=acc_id), exist_ok=True )
    cbmar_subprocess = subprocess.Popen("python blast_amr.py -i cns_{acc_id}.fa -o results/{acc_id}/cbmar/ -d cbmar -p amr_db/cbmar -l 0.9 -t 0.6 ".format(acc_id=acc_id)+ ' && ' + \
            reformat_command+" results/{acc_id}/cbmar/results_tab.txt > results/{acc_id}/cbmar/cbmar.txt".format(acc_id=acc_id), shell=True)

    r_subprocess.communicate()
    c_subprocess.communicate()
    a_subprocess.communicate()
    m_subprocess.communicate()
    cbmar_subprocess.communicate()


def mergeResults(dir_arr):
    #TODO threading
    """
    consolidating the results of all accession IDs
    :param dir_arr:
    :return:
    """
    c_accID = {}
    c_gene = {}
    c_score = {}

    r_accID = {}
    r_gene = {}
    r_score = {}

    a_accID = {}
    a_gene = {}
    a_score = {}

    m_accID = {}
    m_gene = {}
    m_score = {}

    cbmar_accID = {}
    cbmar_gene = {}
    cbmar_score = {}

    for accID in dir_arr:
        accID = accID[:-1]  #removing slash

        card = []
        if os.path.isfile("{accID}/card/card.txt".format(accID=accID)) and os.path.getsize("{accID}/card/card.txt".format(accID=accID)):
            card = subprocess.check_output("cat {accID}/card/card.txt".format(accID=accID),encoding='utf-8', shell=True).strip().split("\n")

        for rec in card:
            
            #col = rec.split("\t")
            col = rec.split("|")
            cdrug = col[0].strip()
            if cdrug not in c_accID:
                c_accID[cdrug] = accID
                c_gene[cdrug] = col[1]
                c_score[cdrug] = col[2]
            else:
                if accID not in c_accID[cdrug]:
                    c_accID[cdrug] = c_accID[cdrug] + "|" + accID
                    c_gene[cdrug] = c_gene[cdrug] + "|"
                    c_score[cdrug] = c_score[cdrug] + "|"
                if c_gene[cdrug][-1] != "|":
                    c_gene[cdrug] = c_gene[cdrug] + ";" + col[1]
                    c_score[cdrug] = c_score[cdrug] + ";" + col[2]
                elif c_gene[cdrug][-1] == "|":
                    c_gene[cdrug] = c_gene[cdrug] + col[1]
                    c_score[cdrug] = c_score[cdrug] + col[2]

        resf = []
        if os.path.isfile("{accID}/resfinder/resf.txt".format(accID=accID)) and os.path.getsize("{accID}/resfinder/resf.txt".format(accID=accID)):
            resf = subprocess.check_output("cat {accID}/resfinder/resf.txt".format(accID=accID),encoding='utf-8', shell=True).strip().split("\n")
        for rec in resf:
            
            col = rec.split("\t")
            rdrug_arr = col[0].split(" resistance")
            rdrug = rdrug_arr[0]
            rdrug = rdrug.lower().strip()
            if rdrug.endswith('s'):
                rdrug = rdrug[:-1]
            # merge drugs if multiple entries
            if rdrug not in r_accID:
                r_accID[rdrug] = accID
                r_gene[rdrug] = col[1]
                r_score[rdrug] = col[2]
            else:
                if accID not in r_accID[rdrug]:
                    r_accID[rdrug] = r_accID[rdrug] + "|" + accID
                    r_gene[rdrug] = r_gene[rdrug] + "|"
                    r_score[rdrug] = r_score[rdrug] + "|"

                if col[1].strip() != "":
                    if r_gene[rdrug][-1] != "|":
                        r_gene[rdrug] = r_gene[rdrug] + ";" + col[1]
                        r_score[rdrug] = r_score[rdrug] + ";" + col[2]
                    elif r_gene[rdrug][-1] == "|":
                        r_gene[rdrug] = r_gene[rdrug] + col[1]
                        r_score[rdrug] = r_score[rdrug] + col[2]
        mega=[]
        if os.path.isfile("{accID}/megares/megares.txt".format(accID=accID)) and os.path.getsize("{accID}/megares/megares.txt".format(accID=accID)):
            mega = subprocess.check_output("cat {accID}/megares/megares.txt".format(accID=accID),encoding='utf-8', shell=True).strip().split("\n")
        for rec in mega:
            
            col = rec.split("\t")
            mdrug = col[0]
            mdrug = mdrug.lower().strip()
            if mdrug.endswith('s'):
                mdrug = mdrug[:-1]
            # merge drugs if multiple entries
            if mdrug not in m_accID:
                m_accID[mdrug] = accID
                m_gene[mdrug] = col[1]
                m_score[mdrug] = col[2].strip()
            else:
                if accID not in m_accID[mdrug]:
                    m_accID[mdrug] = m_accID[mdrug] + "|" + accID
                    m_gene[mdrug] = m_gene[mdrug] + "|"
                    m_score[mdrug] = m_score[mdrug] + "|"

                if col[1].strip() != "":
                    if m_gene[mdrug][-1] != "|":
                        m_gene[mdrug] = m_gene[mdrug] + ";" + col[1]
                        m_score[mdrug] = m_score[mdrug] + ";" + col[2].strip()
                    elif m_gene[mdrug][-1] == "|":
                        m_gene[mdrug] = m_gene[mdrug] + col[1]
                        m_score[mdrug] = m_score[mdrug] + col[2].strip()
        
        cbmar=[]
        if os.path.isfile("{accID}/cbmar/cbmar.txt".format(accID=accID)) and os.path.getsize("{accID}/cbmar/cbmar.txt".format(accID=accID)):
            cbmar = subprocess.check_output("cat {accID}/cbmar/cbmar.txt".format(accID=accID),encoding='utf-8', shell=True).strip().split("\n")
        for rec in cbmar:
            
            col = rec.split("\t")
            cbmardrug = col[0]
            cbmardrug = cbmardrug.lower().strip()
            if cbmardrug.endswith('s'):
                cbmardrug = cbmardrug[:-1]
            # merge drugs if multiple entries
            if cbmardrug not in cbmar_accID:
                cbmar_accID[cbmardrug] = accID
                cbmar_gene[cbmardrug] = col[1]
                cbmar_score[cbmardrug] = col[2].strip()
            else:
                if accID not in cbmar_accID[cbmardrug]:
                    cbmar_accID[cbmardrug] = cbmar_accID[cbmardrug] + "|" + accID
                    cbmar_gene[cbmardrug] = cbmar_gene[cbmardrug] + "|"
                    cbmar_score[cbmardrug] = cbmar_score[cbmardrug] + "|"

                if col[1].strip() != "":
                    if cbmar_gene[cbmardrug][-1] != "|":
                        cbmar_gene[cbmardrug] = cbmar_gene[cbmardrug] + ";" + col[1]
                        cbmar_score[cbmardrug] = cbmar_score[cbmardrug] + ";" + col[2].strip()
                    elif cbmar_gene[cbmardrug][-1] == "|":
                        cbmar_gene[cbmardrug] = cbmar_gene[cbmardrug] + col[1]
                        cbmar_score[cbmardrug] = cbmar_score[cbmardrug] + col[2].strip()



        # merge the amrfinder results:
        amrf = []
        if os.path.isfile("{accID}/amrfinder/amrf.txt".format(accID=accID)) and os.path.getsize("{accID}/amrfinder/amrf.txt".format(accID=accID)):
            amrf = subprocess.check_output("cat {accID}/amrfinder/amrf.txt".format(accID=accID),encoding='utf-8', shell=True).strip().split("\n")
        for rec in amrf:
            
            col = rec.split("\t")
            adrug = col[0].strip().lower()
            if adrug not in a_accID:
                a_accID[adrug] = accID
                a_gene[adrug] = col[1]
                a_score[adrug] = col[2].strip()
            else:
                if accID not in a_accID[adrug]:
                    a_accID[adrug] = a_accID[adrug] + "|" + accID
                    a_gene[adrug] = a_gene[adrug] + "|"
                    a_score[adrug] = a_score[adrug] + "|"
                if a_gene[adrug][-1] != "|":
                    a_gene[adrug] = a_gene[adrug] + ";" + col[1]
                    a_score[adrug] = a_score[adrug] + ";" + col[2].strip()
                elif a_gene[adrug][-1] == "|":
                    a_gene[adrug] = a_gene[adrug] + col[1]
                    a_score[adrug] = a_score[adrug] + col[2].strip()



    # output results
    with open("results.txt", "w") as f:
        f.write(". \t resfinder \t . \t . \t CARD \t . \t . \t megares \t . \t . \t cbmar \t . \t . \t amrfinder \t . \t . \t \n")
        f.write(
            "Antibiotics Ineffective to Bacteria \t AccessionIDs \t Genes \t IDscore \t AccessionIDs \t Genes \t IDscore \t AccessionIDs \t Genes \t IDscore \t AccessionIDs \t Genes \t IDscore \t AccessionIDs \t Genes \t IDscore \n")

        for key, value in list(c_accID.items()):
            if key not in r_accID and key not in m_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\n" % (
                    key, c_accID[key], c_gene[key], c_score[key]))
            elif key not in r_accID and key not in m_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\n" % (
                    key,  c_accID[key], c_gene[key], c_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key in r_accID and key not in m_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\n" % (
                    key, r_accID[key], r_gene[key], r_score[key], c_accID[key], c_gene[key], c_score[key]))
            elif key in r_accID and key not in m_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\n" % (
                    key,  r_accID[key], r_gene[key], r_score[key],c_accID[key], c_gene[key], c_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key not in r_accID and key in m_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\n" % (
                    key, c_accID[key], c_gene[key], c_score[key], m_accID[key], m_gene[key], m_score[key]))
            elif key not in r_accID and key in m_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\n" % (
                    key,  c_accID[key], c_gene[key], c_score[key],m_accID[key], m_gene[key], m_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key in r_accID and key in m_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\n" % (
                    key, r_accID[key], r_gene[key], r_score[key], c_accID[key], c_gene[key], c_score[key], m_accID[key], m_gene[key], m_score[key]))
            elif key in r_accID and key in m_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\n" % (
                    key,  r_accID[key], r_gene[key], r_score[key],c_accID[key], c_gene[key], c_score[key],m_accID[key], m_gene[key], m_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key in r_accID and key in m_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\n" % (
                    key, r_accID[key], r_gene[key], r_score[key], c_accID[key], c_gene[key], c_score[key], m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key in r_accID and key in m_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    key,  r_accID[key], r_gene[key], r_score[key],c_accID[key], c_gene[key], c_score[key],m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key not in r_accID and key in m_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\n" % (
                    key, c_accID[key], c_gene[key], c_score[key], m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key not in r_accID and key in m_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    key,  c_accID[key], c_gene[key], c_score[key], m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key in r_accID and key not in m_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\n" % (
                    key, r_accID[key], r_gene[key], r_score[key], c_accID[key], c_gene[key], c_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key in r_accID and key not in m_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    key,  r_accID[key], r_gene[key], r_score[key],c_accID[key], c_gene[key], c_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))
            elif key not in r_accID and key not in m_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\n" % (
                    key,  c_accID[key], c_gene[key], c_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key not in r_accID and key not in m_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    key,  c_accID[key], c_gene[key], c_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))

        for key, value in list(r_accID.items()):
            if key not in c_accID and key not in m_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n" % (key, r_accID[key], r_gene[key], r_score[key]))
            elif key not in c_accID and key not in m_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\n" % (key, r_accID[key], r_gene[key], r_score[key],a_accID[key],a_gene[key],a_score[key]))

            elif key not in c_accID and key in m_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\n" % (key, r_accID[key], r_gene[key], r_score[key], m_accID[key], m_gene[key], m_score[key]))
            elif key not in c_accID and key in m_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\n" % (key, r_accID[key], r_gene[key], r_score[key], m_accID[key], m_gene[key], m_score[key],a_accID[key],a_gene[key],a_score[key]))


            elif key not in c_accID and key in m_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\n" % (key, r_accID[key], r_gene[key], r_score[key], m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key not in c_accID and key in m_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, r_accID[key], r_gene[key], r_score[key], m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))


            elif key not in c_accID and key not in m_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\n" % (key, r_accID[key], r_gene[key], r_score[key], cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key not in c_accID and key not in m_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, r_accID[key], r_gene[key], r_score[key], cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))

        for key, value in list(m_accID.items()):
            if key not in c_accID and key not in r_accID and key not in cbmar_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t.\t.\t.\n" % (key, m_accID[key], m_gene[key], m_score[key]))
            elif key not in c_accID and key not in r_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\t%s\t%s\t%s\n" % (key, m_accID[key], m_gene[key], m_score[key],a_accID[key],a_gene[key],a_score[key]))


            elif key not in c_accID and key not in r_accID and key in cbmar_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t.\n" % (key, m_accID[key], m_gene[key], m_score[key], cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key not in c_accID and key not in r_accID and key in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, m_accID[key], m_gene[key], m_score[key],cbmar_accID[key],cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))


        for key, value in list(cbmar_accID.items()):
            if key not in c_accID and key not in r_accID and key not in m_accID and key not in a_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t.\t.\t.\n" % (key, cbmar_accID[key],cbmar_gene[key],cbmar_score[key]))
            elif key not in c_accID and key not in r_accID and key not in cbmar_accID and key in a_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, cbmar_gene[key],cbmar_score[key],a_accID[key],a_gene[key],a_score[key]))
        for key, value in list(a_accID.items()):
            if key not in c_accID and key not in r_accID and key not in m_accID and key not in cbmar_accID:
                f.write("%s\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t%s\t%s\t%s\n" % (key, a_accID[key],a_gene[key],a_score[key]))

            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MegaPath-Nano: AMR Detection') 
    parser.add_argument('--query_bam', required=True,help='Input bam')
    parser.add_argument('--output_folder', required=True,help='Output directory')
    parser.add_argument('--taxon', help='Taxon-specific options for AMRFinder, curated organisms: Campylobacter, Enterococcus_faecalis, Enterococcus_faecium, Escherichia, Klebsiella, Salmonella, Staphylococcus_aureus, Staphylococcus_pseudintermedius, Vibrio_cholerae')
    parser.add_argument('--threads', default=psutil.cpu_count(logical=True), help='Num of threads')
    cwd=os.path.dirname(os.path.realpath(__file__))
    nano_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--REFSEQ_PATH', default=nano_dir+'/genomes/refseq/refseq.fna', help='The path of RefSeq')
    FLAGS = parser.parse_args()
    os.makedirs(FLAGS.output_folder,exist_ok=True)
    bam_path = os.path.abspath(FLAGS.query_bam)
    os.chdir(FLAGS.output_folder)
    print("current directory: " + os.getcwd())

    if not os.path.exists("sample.sorted.bam") or not os.path.exists("sample.sorted.bam.bai") or not os.path.exists("header.sam") :
        p = subprocess.Popen("samtools sort {bam_path} -o sample.sorted.bam".format(bam_path=bam_path), shell=True)
        p.communicate()
        p = subprocess.Popen("samtools index sample.sorted.bam", shell=True)

        # make the header for the sam file which will be used for each accession no
        s = subprocess.Popen("samtools view -@ {threads} -H sample.sorted.bam > header.sam".format(threads=FLAGS.threads), shell=True)

        # make sure the original bam is sorted and indexed and header template is ready
        p.communicate()
        print("Original bam file is sorted and indexed")
        s.communicate()
        print("Header template for sam files is ready")

    # make the output directory
    os.makedirs("results",exist_ok=True)
    
    processTaxID(bam_path)

    print("All results have been generated")
    os.chdir("results")
    print("current directory " + os.getcwd())

    try:
        dir_arr = next(os.walk('.'), ([],[],[]))[1]
    except subprocess.CalledProcessError:
        print("No accession ID is mapped with the sequence ID in bam.")
    else:
        print("Processing results of: {acc_list}".format(acc_list=", ".join(dir_arr)))
        mergeResults(dir_arr)
