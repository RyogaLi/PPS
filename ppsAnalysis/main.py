#!/usr/bin/env python3.7

"""
Main script for plasmid pool sequencing analysis

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

"""
import glob

import pandas as pd
import numpy as np
import os
import argparse
import seaborn as sns

import ppsAnalysis.alignment
import ppsAnalysis.cluster
import ppsAnalysis.yeast_variant_analysis
import logging.config
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted


def variants_main(arguments):

    orfs = check_args(arguments)
    # create output folder with user input name
    run_name = arguments.name
    output = os.path.join(arguments.output, run_name)
    if not os.path.isdir(output):
        os.mkdir(output)

    logging.config.fileConfig("./ppsAnalysis/logging.conf")
    main_logger = logging.getLogger("main")
    file_list = os.listdir(arguments.fastq)
    if arguments.align:
        align_log = logging.getLogger("align.log")
        # first align fastq files if user want to use alignment
        all_alignment_jobs = []
        for f in file_list:
            if not f.endswith(".fastq.gz"): continue
            align_log.info(f)
            # for all the fastq files in the dir
            # align the fastq files with given reference
            # note that fastq files and the corresponding reference file has the same id
            if arguments.mode == "human":
                # extract ID
                fastq_ID = f.split("-")[-2]
            #    ref = arguments.ref + "orf9-1_" + fastq_ID
            elif arguments.mode == "yeast":
                fastq_ID = f.split("_")[0]
                # fastq_ID = r1.split("-")[-2]
            #    ref = arguments.ref + "ORF_withpDONR"
            else:
                raise ValueError("Please provide valid mode: human or yeast")
            # mkae sub_output dir for this sample
            sub_output = os.path.join(os.path.abspath(output), fastq_ID)
            if not os.path.isdir(sub_output):
                os.mkdir(sub_output)
            # make sh file for submission in sub_output directory for alignment
            # this is developped for GALEN cluster
            sh_file = os.path.join(sub_output, f"{fastq_ID}.sh")
            f = os.path.join(arguments.fastq, f)
            alignment_obj = ppsAnalysis.alignment.Alignment(arguments.ref, arguments.mode, f, sub_output, sh_file, align_log)
            # the main function writes to the sh file ans submit the file to cluster
            # return job ID
            at = 9
            job_id = alignment_obj.main(at)
            all_alignment_jobs.append(job_id)
        # track all alignment jobs
        alignment_log = logging.getLogger("alignment.log")
        jobs_finished = ppsAnalysis.cluster.parse_jobs_galen(all_alignment_jobs, alignment_log)
        if jobs_finished:
            main_logger.info("Alignment jobs all finished")

    # for each sample, parse vcf files
    all_log = {"fastq_ID": [], "reads":[], "map_perc":[]}
    genes_found = []
    all_genes_summary = pd.DataFrame([], columns=["gene_ID", "gene_len", "total_rd", "avg_rd", "db", "count", "gene_name"])
    all_summary = os.path.join(output, "all_summary_subsetORF.csv")
    all_genes_summary.to_csv(all_summary, index=False)
    
    for f in file_list:
        if not f.endswith(".fastq.gz"): continue
        if arguments.mode == "human":
            # extract ID
            fastq_ID = f.split("-")[-2]
        elif arguments.mode == "yeast":
            fastq_ID = f.split("_")[0]
        else:
            raise ValueError("Wrong mode")
        sub_output = os.path.join(os.path.abspath(output), fastq_ID)

        # there should be only one log file in the dir
        log_file = glob.glob(f"{sub_output}/*.log")[0]

        if not os.path.isfile(log_file):
            main_logger.warning(f"log file does not exist: {log_file}")
            continue
        # get information from the log file to make a summary log file for all the samples
        with open(log_file, "r") as log_f:
            for line in log_f:
                if "reads;" in line:
                    n_reads = line.split(" ")[0]
                    all_log["reads"].append(n_reads)
                if "alignment rate" in line:
                    perc_aligned = line.split("%")[0]
                    all_log["map_perc"].append(perc_aligned)
        all_log["fastq_ID"] += [fastq_ID]*3
        # for each vcf file, get how many genes are fully aligned
        if arguments.mode == "human":
            # extract ID
            pass
        else:  # yeast
            # first get the genes that are fully covered in the fastq files
            orfs_df = orfs[orfs["plate"] == fastq_ID]
            raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_allORFs_raw.vcf")
            if os.path.isfile(raw_vcf_file):
                # analysis of ORFs aligned to all ref + backbone
                fully_covered, stats_list = analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
                fully_covered_file = os.path.join(sub_output, "fully_covered_allORFs.csv")
                fully_covered.to_csv(fully_covered_file, index=False)
                stats_list.append("allORFs")
                genes_found.append(stats_list)
            
            raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_allwithbackbone_raw.vcf")
            if os.path.isfile(raw_vcf_file):
                # analysis of ORFs aligned to all ref
                fully_covered, stats_list = analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
                fully_covered_file = os.path.join(sub_output, "fully_covered_allwithbackbone.csv")
                fully_covered.to_csv(fully_covered_file, index=False)
                stats_list = stats_list.append("allwithbackbone")
                genes_found.append(stats_list)

            raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_plateORFs_raw.vcf")
            if os.path.isfile(raw_vcf_file):
                # analysis of ORFs aligned to subgroup
                fully_covered, stats_list = analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
                fully_covered_file = os.path.join(sub_output, "fully_covered_plateORFs.csv")
                fully_covered.to_csv(fully_covered_file, index=False)
                fully_covered.to_csv(all_summary, index=False, header=False, mode="a")
                stats_list.append("plateORFs")
                genes_found.append(stats_list)

            raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_subsetORFs_raw.vcf")
            if not os.path.isfile(raw_vcf_file):
                # analysis of ORFs aligned to plate
                fully_covered, stats_list = analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
                fully_covered_file = os.path.join(sub_output, "fully_covered_subsetORFs.csv")
                fully_covered.to_csv(fully_covered_file, index=False)
                stats_list.append("subsetORFs")
                genes_found.append(stats_list)

    # process all log
    all_log = pd.DataFrame(all_log)
    all_log_file = os.path.join(output, "alignment_log.csv")
    all_log.to_csv(all_log_file, index=False)

    # process summary of number of genes found in each sample
    all_genes_stats = pd.DataFrame(genes_found, columns=["plate", "fully_aligned", "all_genes_found",
                                                         "all_targeted_on_plate", "all_targeted_full",
                                                         "n_genes_with_any_mut", "n_ref", "aligned_to"])
    genes_found_file = os.path.join(output, "genes_stats.csv")

    all_genes_stats["% on plate fully aligned"] = all_genes_stats["all_targeted_full"]/all_genes_stats["all_targeted_on_plate"]
    all_genes_stats.to_csv(genes_found_file, index=False)
    sns.set_theme(style="whitegrid", font_scale=1.5)
    plt.figure(figsize=(20,14))
    g = sns.barplot(data = all_genes_stats, x = "plate", y = "% on plate fully aligned", hue="aligned_to")
    plt.xticks(rotation=90, fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylabel("all the fully aligned unique ORFs", fontsize=15)
    plt.tight_layout()
    plt.savefig(os.path.join(output, "./perc_full_matched.png"))
    plt.close()

    # compare genes in all the targeted space (ORFs) vs all fully aligned 
    all_targeted_unique_db = orfs["ORF_NAME_NODASH"].dropna().unique()
    all_found = pd.read_csv(all_summary)
    all_found["gene_name"] = all_found["gene_name"].replace("-", "")
    all_found_genes = all_found["gene_name"].dropna().unique()
    venn2([set(all_targeted_unique_db), set(all_found_genes)], set_labels = ("all ORFs", "all_fully_aligned"))
    plt.savefig(os.path.join(output, "./allORFs_venn.png"))
    plt.close()

    # HIP subset
    all_HIP_targeted = orfs[orfs["db"] == "HIP"]["ORF_NAME_NODASH"].dropna().unique() 
    all_found_hip = all_found[all_found["db"] == "HIP"]
    all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
    all_found_genes = all_found_hip["gene_name"].dropna().unique()
    venn2([set(all_HIP_targeted), set(all_found_genes)], set_labels = ("all HIP ORFs", "all_fully_aligned"))
    plt.savefig(os.path.join(output, "./HIPORFs_venn.png"))
    plt.close()

    # SGD subset
    all_HIP_targeted = orfs[orfs["db"] == "SGD"]["ORF_NAME_NODASH"].dropna().unique() 
    all_found_hip = all_found[all_found["db"] == "SGD"]
    all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
    all_found_genes = all_found_hip["gene_name"].dropna().unique()
    venn2([set(all_HIP_targeted), set(all_found_genes)], set_labels = ("all SGD ORFs", "all_fully_aligned"))
    plt.savefig(os.path.join(output, "./SGDORFs_venn.png"))
    plt.close()

    # PROTGEN subset
    all_HIP_targeted = orfs[orfs["db"] == "PROTGEN"]["ORF_NAME_NODASH"].dropna().unique() 
    all_found_hip = all_found[all_found["db"] == "PROTGEN"]
    all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
    all_found_genes = all_found_hip["gene_name"].dropna().unique()
    venn2([set(all_HIP_targeted), set(all_found_genes)], set_labels = ("all PROTGEN ORFs", "all_fully_aligned"))
    plt.savefig(os.path.join(output, "./PROTORFs_venn.png"))
    plt.close()


def read_yeast_csv(HIP_target_ORFs, other_target_ORFs):
    """
    Join HIP data and other data into one df, remove unwanted columns
    :param HIP_target_ORFs: csv file contains which HIP ORF is in which sample
    :param other_target_ORFs: csv file contains which other ORF is in which sample
    :return: df with ORF name, db name and sample name
    """
    HIP_df = pd.read_csv(HIP_target_ORFs)
    other_target_ORFs = pd.read_csv(other_target_ORFs)

    HIP_df = HIP_df[["ORF_id", "ORF_NAME_NODASH", "len(seq)", "SYMBOL", "plate"]]
    HIP_df["db"] = "HIP"
    HIP_df = HIP_df.rename(columns={"ORF_id": "orf_name"})
    other_ORFs = other_target_ORFs[["orf_name", "ORF_NAME_NODASH", "src_collection", "plate"]]
    other_ORFs = other_ORFs.rename(columns={"src_collection": "db"})
    #other_ORFs['plate'] = 'scORFeome-' + other_ORFs['plate'].astype(str)
    combined = pd.concat([HIP_df, other_ORFs], axis=0, ignore_index=True)
    return combined


def analysisYeast(raw_vcf_file, fastq_ID, orfs_df):
    """
    Run yeast variants analysis, make files for each plate and save to corresponded dir
    also return dfs for combining
    """
    analysis = ppsAnalysis.yeast_variant_analysis.yeastAnalysis(raw_vcf_file, fastq_ID, orfs_df)
    full_cover_genes, gene_dict, ref_dict = analysis.get_full_cover()
        
    # all the genes with full coverage
    n_fully_aligned = len(full_cover_genes.keys())
    # all genes in ref fasta
    n_ref = len(ref_dict.keys())
    # all genes found in this fastq file 
    n_all_found = len(gene_dict.keys())

    # save all the genes that are fully covered to the output folder
    fully_covered = pd.DataFrame.from_dict(full_cover_genes, orient='index').reset_index()
    fully_covered.columns = ["gene_ID", "gene_len", "total_rd", "avg_rd"]
    # split gene ID col 
    #fully_covered["gene_ID"] = fully_covered["gene_ID"].str.replace(to_replace="-[A-G]", "A")
     
    #fully_covered = fully_covered.replace(to_replace ='-index[0-9]+', value = '', regex = True)
    fully_covered["db"] = fully_covered["gene_ID"].str.extract(r".*-([A-Z]+)-[1-9]")
    fully_covered["count"] = fully_covered["gene_ID"].str.extract(r".*-[A-Z]+-([1-9])")
    fully_covered["gene_name"] = fully_covered["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
    # filter vcf based on QUAL and DP
    mut_count_dict = analysis.filter_vcf()
    gene_mut_count = pd.DataFrame.from_dict(mut_count_dict, orient='index').reset_index()
    gene_mut_count.columns = ["gene_name", "n_variants"]
    n_mut_genes = len(mut_count_dict)

    # merge with target orfs
    merged_df = pd.merge(orfs_df, fully_covered, how="left", left_on="orf_name", right_on="gene_ID")
    # merged_file = os.path.join(sub_output, "merged_with_targets.csv")
    # merged_df.to_csv(merged_file, index=False)
    # merged_df.to_csv(all_summary, mode="a", index=False, header=False)
    n_targeted = orfs_df.shape[0]
    n_targeted_full = merged_df[~merged_df["gene_ID"].isnull()].shape[0]
    stats_list = [fastq_ID, n_fully_aligned, n_all_found, n_targeted, n_targeted_full,n_mut_genes, n_ref]

    # # merge with ref to get gene len
    # ref = pd.DataFrame.from_dict(ref_dict, orient='index').reset_index()
    # ref.columns = ["gene_name", "gene_len"]
    # merge_mut_count = pd.merge(gene_mut_count, ref, how="left", on="gene_name")
    # # split gene name col
    # merge_mut_count["gene_ID"] = merge_mut_count["gene_name"].str.replace("gene", "")
    # merge_mut_count = merge_mut_count.replace(to_replace ='-index[0-9]+', value = '', regex = True)
    # # merge this with targeted ORFs
    # target_gene_mut_count = pd.merge(orfs_df, merge_mut_count, how="left", left_on="orf_name", right_on="gene_ID")
    # # save to file
    # target_gene_mut_file = os.path.join(sub_output, "target_gene_mutcount.csv")
    # target_gene_mut_count.to_csv(target_gene_mut_file, index=False)

    return fully_covered, stats_list

def check_args(arguments):
    """
    Check user input arguments
    :param arguments: user input arguments, argparse struct
    :return: None
    """

    if not os.path.isdir(arguments.output):
        raise NotADirectoryError(f"{arguments.output} does not exist")

    if arguments.fastq and not os.path.isdir(arguments.fastq):
        raise NotADirectoryError(f"{arguments.fastq} does not exist")

    if arguments.align:
        # if alignment, must also provide path to fastq files and reference files
        if not arguments.fastq or not arguments.ref:
            raise ValueError("Please also provide reference dir and fastq dir")

    if arguments.mode == "yeast":
        # load yeast data
        # files contain ORF information in
        HIP_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/HIP_targeted_ORFs.csv"
        other_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/other_targeted_ORFs.csv"
        orfs = read_yeast_csv(HIP_target_ORFs, other_target_ORFs)
    elif arguments.mode == "human":
        orfs =""
    else:
        exit()

    return orfs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('--align', action="store_true", help='provide this argument if users want to start with '
                                                         'alignment, otherwise the program assumes alignment was '
                                                         'done and will analyze the vcf files.')
    parser.add_argument("-f", "--fastq", help="input fastq files")
    parser.add_argument("-m", "--mode", help="Human or Yeast PPS?")
    parser.add_argument("-r", "--ref", help="Path to referece files")
    parser.add_argument('-o', "--output", help='Output directory', required=True)
    parser.add_argument('-n', "--name", help='Name for this run', required=True)
    args = parser.parse_args()

    variants_main(args)
