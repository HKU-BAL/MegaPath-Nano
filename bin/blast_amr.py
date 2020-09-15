import resfinder.resfinder
from argparse import ArgumentParser
from cgecore.blaster import Blaster
import collections
import re
from tabulate import tabulate

class Blast_amr(resfinder.resfinder.ResFinder):

    def __init__(self, databases,db_path):
        self.databases=databases
        self.db_path=db_path
        self.phenos=dict()

    def write_results(self, out_path, result):
        result_str = self.results_to_str(
                                             results=result.results,
                                             databases=self.databases,
                                             query_align=result.gene_align_query,
                                             homo_align=result.gene_align_homo,
                                             sbjct_align=result.gene_align_sbjct)

        with open(out_path + "/results_tab.txt", "w") as fh:
            fh.write(result_str[0])
        with open(out_path + "/results_table.txt", "w") as fh:
            fh.write(result_str[1])
        with open(out_path + "/results.txt", "w") as fh:
            fh.write(result_str[2])
        with open(out_path + "/Resistance_gene_seq.fsa", "w") as fh:
            fh.write(result_str[3])
        with open(out_path + "/Hit_in_genome_seq.fsa", "w") as fh:
            fh.write(result_str[4])

    def results_to_str(self, results, databases, query_align=None,
                       homo_align=None, sbjct_align=None):


        # Write the header for the tab file
        tab_str_list=[]
        tab_str_list.append("Resistance gene\tIdentity\tAlignment Length/Gene Length\t"
                   "Coverage\tPosition in reference\tContig\t"
                   "Position in contig\tPhenotype\tAccession no.\n")

        table_str_list = []
        txt_str_list = []
        ref_str_list = []
        hit_str_list = []

        # Getting and writing out the results
        titles = dict()
        rows = dict()
        headers = dict()
        txt_file_seq_text = dict()
        split_print = collections.defaultdict(list)

        for db in results:
            if (db == "excluded"):
                continue

            # Clean up dbs with only excluded hits
            if isinstance(results[db], str):
                results[db] = "No hit found"
                print("No hits in: " + str(db))
                continue

            no_hits = True
            for hit in results[db]:
                if (hit not in results["excluded"]):
                    print("Hit not excluded: " + str(hit))
                    print("\tIn DB: " + str(db))
                    no_hits = False

            if (no_hits):
                print("No hits in: " + str(db) + " (hit(s) excluded)")
                results[db] = "No hit found"

            profile = str(db)
            if results[db] == "No hit found":
                table_str_list.append("%s\n%s\n\n" % (profile, results[db]))
            else:
                titles[db] = "%s" % (profile)
                headers[db] = ["Resistance gene", "Identity",
                               "Alignment Length/Gene Length", "Coverage",
                               "Position in reference", "Contig",
                               "Position in contig", "Phenotype",
                               "Accession no."]
                table_str_list.append("%s\n" % (profile))
                table_str_list.append("Resistance gene\tIdentity\t"
                              "Alignment Length/Gene Length\tCoverage\t"
                              "Position in reference\tContig\tPosition in contig\t"
                              "Phenotype\tAccession no.\n")

                rows[db] = list()
                txt_file_seq_text[db] = list()

                for hit in results[db]:
                    if (hit in results["excluded"]):
                        continue

                    res_header = results[db][hit]["sbjct_header"]
                    gene=res_header
                    ID = results[db][hit]["perc_ident"]
                    coverage = results[db][hit]["perc_coverage"]
                    sbjt_length = results[db][hit]["sbjct_length"]
                    HSP = results[db][hit]["HSP_length"]
                    positions_contig = "%s..%s" % (results[db][hit]["query_start"],
                                                   results[db][hit]["query_end"])
                    positions_ref = "%s..%s" % (results[db][hit]["sbjct_start"],
                                                results[db][hit]["sbjct_end"])
                    contig_name = results[db][hit]["contig_name"]
                    if 'cbmar' in databases:
                        pheno="Betalactamase"
                    elif 'megares' in databases:
                        pheno=""
                        drug_class=["Aminocoumarins", "Aminoglycosides", "Bacitracin", "betalactams", "Cationic_antimicrobial_peptides", "Elfamycins", "Fluoroquinolones", "Fosfomycin", "Fusidic_acid", "Glycopeptides", "Lipopeptides", "Metronidazole", "MLS", "Multi-drug_resistance", "Mycobacterium_tuberculosis-specific_Drug", "Phenicol", "Rifampin", "Sulfonamides", "Tetracyclines", "Thiostrepton", "Trimethoprim", "Tunicamycin"]
                        for drug in drug_class:
                            if drug in res_header:
                                pheno = drug
                    
                    acc = re.search("(NZ_){0,1}[a-zA-Z]{1,6}_{0,1}[0-9]{5,9}\.{0,1}[0-9]{0,2}", res_header).group()
                    if not acc:
                        acc = res_header

                    if "split_length" in results[db][hit]:
                        total_HSP = results[db][hit]["split_length"]
                        split_print[res_header].append([gene, ID, total_HSP,
                                                        sbjt_length, coverage,
                                                        positions_ref, contig_name,
                                                        positions_contig, pheno,
                                                        acc])
                        tab_str_list.append("%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                    % (gene, ID, HSP, sbjt_length, coverage,
                                       positions_ref, contig_name, positions_contig,
                                       pheno, acc)
                                    )
                    else:
                        # Write tabels
                        table_str_list.append("%s\t%.2f\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                      % (gene, ID, HSP, sbjt_length, coverage,
                                         positions_ref, contig_name,
                                         positions_contig, pheno, acc)
                                      )
                        tab_str_list.append("%s\t%.2f\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                    % (gene, ID, HSP, sbjt_length, coverage,
                                       positions_ref, contig_name, positions_contig,
                                       pheno, acc)
                                    )

                        # Saving the output to write the txt result table
                        hsp_length = "%s/%s" % (HSP, sbjt_length)
                        rows[db].append([gene, ID, hsp_length, coverage,
                                         positions_ref, contig_name, positions_contig,
                                         pheno, acc])

                    # Writing subjet/ref sequence
                    ref_seq = sbjct_align[db][hit]

                    ref_str_list.append(">%s_%s\n" % (gene, acc))
                    for i in range(0, len(ref_seq), 60):
                        ref_str_list.append("%s\n" % (ref_seq[i:i + 60]))

                    # Getting the header and text for the txt file output
                    sbjct_start = results[db][hit]["sbjct_start"]
                    sbjct_end = results[db][hit]["sbjct_end"]
                    text = ("%s, ID: %.2f %%, Alignment Length/Gene Length: %s/%s, "
                            "Coverage: %s, "
                            "Positions in reference: %s..%s, Contig name: %s, "
                            "Position: %s" % (gene, ID, HSP, sbjt_length, coverage,
                                              sbjct_start, sbjct_end, contig_name,
                                              positions_contig))
                    hit_str_list.append(">%s\n" % text)

                    # Writing query/hit sequence
                    hit_seq = query_align[db][hit]

                    for i in range(0, len(hit_seq), 60):
                        hit_str_list.append("%s\n" % (hit_seq[i:i + 60]))

                    # Saving the output to print the txt result file allignemts
                    txt_file_seq_text[db].append((text, ref_seq,
                                                      homo_align[db][hit], hit_seq))

                for res in split_print:
                    gene = split_print[res][0][0]
                    ID = split_print[res][0][1]
                    HSP = split_print[res][0][2]
                    sbjt_length = split_print[res][0][3]
                    coverage = split_print[res][0][4]
                    positions_ref = split_print[res][0][5]
                    contig_name = split_print[res][0][6]
                    positions_contig = split_print[res][0][7]
                    pheno = split_print[res][0][8]
                    acc = split_print[res][0][9]

                    total_coverage = 0

                    for i in range(1, len(split_print[res])):
                        ID = "%s, %.2f" % (ID, split_print[res][i][1])
                        total_coverage += split_print[res][i][4]
                        positions_ref = positions_ref + ", " + split_print[res][i][5]
                        contig_name = contig_name + ", " + split_print[res][i][6]
                        positions_contig = (positions_contig + ", "
                                            + split_print[res][i][7])

                    table_str_list.append("%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                  % (gene, ID, HSP, sbjt_length, coverage,
                                     positions_ref, contig_name, positions_contig,
                                     pheno, acc)
                                  )

                    hsp_length = "%s/%s" % (HSP, sbjt_length)

                    rows[db].append([gene, ID, hsp_length, coverage, positions_ref,
                                     contig_name, positions_contig, pheno, acc])

                table_str_list.append("\n")


        # Writing table txt for all hits
        for db in titles:
            # Txt file table
            table = Blast_amr.text_table(titles[db], headers[db], rows[db])
            txt_str_list.append(table)

        # Writing alignment txt for all hits
        for db in titles:
            # Txt file alignments
            txt_str_list.append("##################### %s #####################\n"
                        % (db))
            for text in txt_file_seq_text[db]:
                txt_str_list.append("%s\n\n" % (text[0]))
                for i in range(0, len(text[1]), 60):
                    txt_str_list.append("%s\n" % (text[1][i:i + 60]))
                    txt_str_list.append("%s\n" % (text[2][i:i + 60]))
                    txt_str_list.append("%s\n\n" % (text[3][i:i + 60]))
                txt_str_list.append("\n")

        tab_str=''.join(tab_str_list)
        table_str=''.join(table_str_list)
        txt_str=''.join(txt_str_list)
        ref_str=''.join(ref_str_list)
        hit_str=''.join(hit_str_list)
        return (tab_str, table_str, txt_str, ref_str, hit_str)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--inputfile",
                        dest="inputfile",
                        help="Input file",
                        default=None)
    parser.add_argument("-o", "--outputPath",
                        dest="out_path",
                        help="Path to blast output",
                        default='')
    parser.add_argument("-p", "--databasePath",
                        dest="db_path",
                        help="Path to the databases",
                        default='')
    parser.add_argument("-d", "--databases",
                        dest="databases",
                        help="Database name without suffix, separated by comma - if none are \
                              specified all are used. Databases are expected to have fsa suffix",
                        default=None)
    parser.add_argument("-l", "--min_cov",
                        dest="min_cov",
                        help="Minimum coverage",
                        default=0.60)
    parser.add_argument("-t", "--threshold",
                        dest="threshold",
                        help="Blast threshold for identity",
                        default=0.90)
    args = parser.parse_args()

    db_list=args.databases.split(",")
    
    blast_amr=Blast_amr(databases=db_list, db_path=args.db_path)
    blast_run = Blaster(inputfile=args.inputfile, databases=db_list, db_path=args.db_path,
                            out_path=args.out_path, min_cov=args.min_cov,
                            threshold=args.threshold)
    blast_amr.write_results(out_path=args.out_path,result=blast_run)
