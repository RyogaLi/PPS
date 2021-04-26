from __future__ import division
from conf import *
import pandas as pd

def count_gene(summary, plate_col_name, gene_col_name, snp_col=None, indel_col=None, alignment_rate=None):
    """
    for each plate in summary, count all the genes
    :param summary: summary file
    :param col: column that contains gene name
    :return:
    """

    df = pd.read_csv(summary,error_bad_lines=False)
    df.dropna(axis=1)

    plate_names = set(df[plate_col_name])
    # print len(plate_names)
    all_plates = {}
    for i in plate_names:
        # print i
        if (snp_col != None) and (alignment_rate != None): # select only the one with snp
            df[alignment_rate] = df[alignment_rate].apply(pd.to_numeric, errors='coerce')
            genes = list(df.loc[(df[plate_col_name] == i) & (df[snp_col] == "0") & (df[indel_col] == "0") & (
                df[alignment_rate] >= 0.9)][gene_col_name])
        elif alignment_rate != None:
            df[alignment_rate] = df[alignment_rate].apply(pd.to_numeric, errors='coerce')
            genes = list(df.loc[(df[plate_col_name] == i) & (df[alignment_rate] >= 0.9)][gene_col_name])
        else:
            genes = list(df.loc[df[plate_col_name] == i][gene_col_name])

        clean_genes = [x for x in genes if str(x) != 'nan']
        all_plates[i] = clean_genes
    return all_plates


def normalize_keys(dictionary, regex):
    """
    select string in regex and replace keys with a match
    :param dictionary:
    :param regex:
    :return:
    """
    for key in dictionary.keys():
        up = key.upper()
        try:
            result = re.match(regex, up)
            result = result.group(1)
            if "_0" in result:
                result = result.replace("_0","")
        except Exception as e:
            # print(key)
            pass
        dictionary[result] = dictionary.pop(key)
    return dictionary


def find_overlap(d1_ref, d2):
    """
    find overlap between values with the same key
    :param d1:
    :param d2:
    :return:
    """
    summary = {}
    for key in d1_ref.keys():
        ref_genes = set(map(int, d1_ref[key]))
        genes = set(map(int, d2[key]))
        print ref_genes
        print genes
        # number of genes we got
        recovered_genes = len(list(ref_genes&genes))
        print list(ref_genes&genes)
        # number of all the genes
        all_genes = len(ref_genes)

        summary[key] = [recovered_genes,all_genes]
    return summary


def reads_count(fastq_file_list):
    """
    count number of reads in each fastq file
    :param fastq: path to fastq file
    :return:
    """
    rc = {}

    for file in os.listdir(fastq_file_list):
        line = os.path.join(fastq, file)
        cmd = "cat "+line+" | "+ "wc"+" -l"
        temp = int(os.popen(cmd).read())/4
        plate_name = line.split("/")[-1].split("_")[0]
        if "Sup" in plate_name:
            plate_name = plate_name.split("-")[1]
        rc[plate_name] = temp
    return rc

def yeast_pps_main():
    ref_summary_hip = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/csv/hip_plates.csv"
    ref_summary_sup = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/csv/sup_plates.csv"

    hip_genes = count_gene(ref_summary_hip, "384-Pool Name", "ORF_NAME_NODASH")
    sup_genes = count_gene(ref_summary_sup, "Condensed plate", "orf_name")

    # merge 2 dict into one
    hip_genes.update(sup_genes)

    # get gene summary from pps summary
    pps_summary = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/csv/20170927_filtered_synonymous_summary.csv"
    # with_ambiguous = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/csv/filtered_synonymous_summary.csv"
    pps_genes_removed_snp = count_gene(pps_summary, "plate_name", "gene_name", "number_of_SNP", "number_of_INDEL",
                                       "alignment_rate")
    pps_genes_align = count_gene(pps_summary, "plate_name", "gene_name", alignment_rate="alignment_rate")
    pps_genes = count_gene(pps_summary, "plate_name", "gene_name")

    # print hip_genes.keys()
    # normalize the keys (plate names)
    regex = ".*((HIP|SUP)[-|_]?[0-9]*).*"
    hip_genes = normalize_keys(hip_genes, regex)
    pps_genes = normalize_keys(pps_genes, regex)
    pps_genes_removed_snp = normalize_keys(pps_genes_removed_snp, regex)
    pps_genes_align = normalize_keys(pps_genes_align, regex)

    # find overlap (with snp)
    overlap_dict = find_overlap(hip_genes, pps_genes)
    # calculate percentage recovered
    # plate_names = overlap_dict.keys()
    percent_recovered = []
    for key in overlap_dict:
        percent_recovered.append((key, overlap_dict[key][0] / overlap_dict[key][1]))
    percent_recovered = dict(percent_recovered)

    # find overlap (without snp)
    overlap_dict = find_overlap(hip_genes, pps_genes_removed_snp)
    # calculate percentage recovered
    # plate_names = overlap_dict.keys()
    percent_recovered_no_snp = []
    for key in overlap_dict:
        percent_recovered_no_snp.append((key, overlap_dict[key][0] / overlap_dict[key][1]))
    percent_recovered_no_snp = dict(percent_recovered_no_snp)

    # find overlap (without snp)
    overlap_dict = find_overlap(hip_genes, pps_genes_align)
    # calculate percentage recovered
    # plate_names = overlap_dict.keys()
    percent_recovered_align = []
    for key in overlap_dict:
        percent_recovered_align.append((key, overlap_dict[key][0] / overlap_dict[key][1]))
    percent_recovered_align = dict(percent_recovered_align)

    # read counts
    # fastq_files = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/fastq_list.txt"
    # all_reads = reads_count(fastq_files)
    # # print all_reads.keys()
    # #normalize keys
    # all_reads = normalize_keys(all_reads, regex)

    # percentage of genes recovered
    # x_keys =  sorted(all_reads, key=all_reads.get)
    # x = []
    # y = []
    # for key in x_keys:
    # 	x.append(percent_recovered_no_snp[key])
    # 	# print all_reads[key]
    # 	y.append(all_reads[key])
    # 	# labels.append(key.split("-")[-1])
    # # plot from (overlap and reads count)
    # plt.plot(x, y, ".")
    # # plt.xticks(x, x_keys, rotation=0, fontsize=6)
    # plt.ylabel("Total number of reads")
    # plt.xlabel("Percentage of genes found")
    # plt.savefig("read_depth.png")
    # plt.close()

    sorted_plates = sorted(percent_recovered, key=percent_recovered.get)
    sorted_values = [percent_recovered[i] for i in sorted_plates]

    sorted_plates_no_mut = sorted(percent_recovered_no_snp, key=percent_recovered_no_snp.get)
    sorted_values_no_mut = [percent_recovered_no_snp[i] for i in sorted_plates_no_mut]

    sorted_plates_align = sorted(percent_recovered_align, key=percent_recovered_align.get)
    sorted_values_align = [percent_recovered_align[i] for i in sorted_plates_align]

    print sorted_values_no_mut
    print sorted_values_align == sorted_values

    plt.figure(figsize=(15, 7.5))
    plt.plot(range(len(sorted_plates)), sorted_values, ".", ms=14, label="with no filter applied")
    plt.plot(range(len(sorted_plates)), sorted_values_no_mut, ".", ms=14, label="with all filter applied")
    plt.plot(range(len(sorted_plates)), sorted_values_align, ".", ms=14, label="alignment rate >= 90%")
    plt.legend(loc=4)
    plt.ylim((0, 1))
    plt.xlim((0, 23))
    plt.xlabel("plates", fontsize=14)
    plt.ylabel("percentage of genes recovered", fontsize=16)
    plt.xticks(range(len(sorted_plates)), sorted_plates, rotation=20)
    plt.savefig("compare_mutations.png")
    plt.close()

def huri_pps_main():
    huri_ref = "/Users/roujia/Documents/02_dev/01_KiloSEQ/PASS_3/20171011_bhORFeome_P3_entry_withseqs.csv"
    huri_summary = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_missing_orffiltered_synonymous_summary.csv"
    pps_fastq = "/Users/roujia/Documents/01_ngsdata/20171005_Pass3_PPS/combined/"

    ref_genes = count_gene(huri_ref, "384_plate_name", "orf_id")
    pps_genes = count_gene(huri_summary, "plate_name", "gene_name")


    # normalize keys
    regex = ".+((AD|DB|both)(_?0?[0-9]?)).?"
    ref_genes = normalize_keys(ref_genes, regex)
    pps_genes = normalize_keys(pps_genes, regex)

    overlap_dict = find_overlap(ref_genes, pps_genes)

    percent_recovered = []
    for key in overlap_dict:
        percent_recovered.append((key, overlap_dict[key][0] / overlap_dict[key][1]))
    percent_recovered = dict(percent_recovered)
    print percent_recovered

    # read count
    all_reads = reads_count(fastq)
    # print all_reads.keys()
    #normalize keys
    all_reads = normalize_keys(all_reads, regex)

    # percentage of genes recovered
    x_keys =  sorted(all_reads, key=all_reads.get)
    x = []
    y = []
    for key in x_keys:
        x.append(percent_recovered[key])
        # print all_reads[key]
        y.append(all_reads[key])
        # labels.append(key.split("-")[-1])
    # plot from (overlap and reads count)
    plt.plot(x, y, ".")
    # plt.xticks(x, x_keys, rotation=0, fontsize=6)
    plt.ylabel("Total number of reads")
    plt.xlabel("Percentage of genes found")
    plt.savefig("read_depth.png")
    plt.close()

def missing_clones_main():
    missing_clones_ref = ""
    PPS_summary = ""

    ref_genes = ""
    pps_genes = ""


    pass

if __name__ == '__main__':
    # yeast_pps_main()
    # huri_pps_main()
    missing_clones_main()
