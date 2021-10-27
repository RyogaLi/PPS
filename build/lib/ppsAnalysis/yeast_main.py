#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import glob

import pandas as pd
import numpy as np
import os
import argparse
import seaborn as sns

import ppsAnalysis.alignment
import ppsAnalysis.cluster
import ppsAnalysis.yeast_variant_analysis
import ppsAnalysis.human_variant_analysis
import ppsAnalysis.logthis
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted

"""
DEPRECIATED
"""

def variants_main(arguments):
    """
    Main function for pps analysis for yeast data
    :param arguments: user input arguments
    :return: None
    """
    # set log level
    loglevel = "DEBUG"
    orfs = check_args(arguments)
    # create output folder with user input name
    run_name = arguments.name
    output = os.path.join(arguments.output, run_name)
    if not os.path.isdir(output):
        os.mkdir(output)
    # make main log file
    main_log = os.path.join(output, "main.log")
    log_obj = ppsAnalysis.logthis.logit(log_f=main_log, log_level=loglevel)
    main_logger = log_obj.get_logger("main")
    # list all the files from input fastq files
    # the fastq IDs are used to create output folders for each plate
    file_list = os.listdir(arguments.fastq)
    # if user want to run alignment
    if arguments.align:
        align_log = log_obj.get_logger("align.log")
        # first align fastq files if user want to use alignment
        all_alignment_jobs = []
        for f in file_list:
            if not f.endswith(".fastq.gz"): continue
            align_log.info(f)
            # for all the fastq files in the dir
            # align the fastq files with given reference
            # note that fastq files and the corresponding reference file has the same id
            fastq_ID = f.split("_")[0]
            # make sub_output dir for this sample
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
        cluster_log = log_obj.get_logger("alignment.log")
        jobs_finished = ppsAnalysis.cluster.parse_jobs_galen(all_alignment_jobs, cluster_log)
        if jobs_finished:
            main_logger.info("Alignment jobs all finished")

    parse_vcf_files_yeast(output, file_list, orfs, main_logger)


def parse_vcf_files_yeast(output, file_list, orfs, logger):
    """
    Parse vcf files for each sample
    :param output:
    :param file_list: list of fastq files
    :param orfs: dataframe contains all the targeted ORFs
    :param logger: logging object
    :return:
    """
    # keep track of alignment log
    all_log = {"fastq_ID": [], "reads": [], "map_perc": []}
    genes_found = []
    # empty df to save all the variants
    # all_summary_file = os.path.join(output, "all_summary.csv")

    all_summary_plate, all_summary_all, all_summary_subset = [], [], []
    all_mut_df_plate, all_mut_df_all, all_mut_df_subset = [], [], []
    for f in file_list:
        if not f.endswith(".fastq.gz"): continue
        fastq_ID = f.split("_")[0]
        sub_output = os.path.join(os.path.abspath(output), fastq_ID)

        # get the latest log file from dir
        log_file = glob.glob(f"{sub_output}/*.log")[-1]
        if not os.path.isfile(log_file):
            logger.warning(f"log file does not exist: {fastq_ID}")
            continue

        # get information from the log file to make a summary log file for all the samples
        with open(log_file, "r") as log_f:
            n_sample = 0
            for line in log_f:
                if "reads;" in line:
                    n_sample +=1
                    n_reads = line.split(" ")[0]
                    all_log["reads"].append(n_reads)
                if "alignment rate" in line:
                    perc_aligned = line.split("%")[0]
                    all_log["map_perc"].append(perc_aligned)

        all_log["fastq_ID"] += [fastq_ID] * n_sample
        # for each vcf file, get how many genes are fully aligned
        orfs_df = orfs[orfs["plate"] == fastq_ID]
        raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_plateORFs_raw.vcf")
        if os.path.isfile(raw_vcf_file):
            # analysis of ORFs aligned to subgroup
            all_summary_df, stats_list, mut_df= analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
            # fully_covered_file = os.path.join(sub_output, "fully_covered_plateORFs.csv")
            # fully_covered.to_csv(fully_covered_file, index=False)
            all_group_summary_file = os.path.join(sub_output, "all_summary_plateORFs.csv")
            all_summary_df.to_csv(all_group_summary_file, index=False)            # all_found.to_csv(all_found_file, index=False)
            # fully_covered.to_csv(all_summary, index=False, header=False, mode="a")
            # all_found.to_csv(all_found_summary, index=False, header=False, mode="a")

            db = all_summary_df["db"].unique()
            stats_list.append("plateORFs")
            genes_found.append(stats_list)
            mut_df["plate"] = fastq_ID
            mut_df["db"] = db[0]
            all_mut_df_plate.append(mut_df)
            all_summary_df["db"] = db[0]
            all_summary_plate.append(all_summary_df)

        raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_allORFs_raw.vcf")
        if os.path.isfile(raw_vcf_file):
            # analysis of ORFs aligned to subgroup
            all_summary_df, stats_list, mut_df= analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
            # fully_covered_file = os.path.join(sub_output, "fully_covered_plateORFs.csv")
            # fully_covered.to_csv(fully_covered_file, index=False)

            all_group_summary_file = os.path.join(sub_output, "all_summary_allORFs.csv")
            all_mut_file = os.path.join(sub_output, "all_mutation_allORFs.csv")

            all_summary_df.to_csv(all_group_summary_file, index=False)            # all_found.to_csv(all_found_file, index=False)
            # fully_covered.to_csv(all_summary, index=False, header=False, mode="a")
            # all_found.to_csv(all_found_summary, index=False, header=False, mode="a")

            db = all_summary_df["db"].unique()
            stats_list.append("allORFs")
            genes_found.append(stats_list)
            mut_df["plate"] = fastq_ID
            mut_df["db"] = db[0]
            all_mut_df_all.append(mut_df)
            all_summary_df["db"] = db[0]
            all_summary_all.append(all_summary_df)
        
        raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_subsetORFs_raw.vcf")
        if os.path.isfile(raw_vcf_file):
            # analysis of ORFs aligned to subgroup
            all_summary_df, stats_list, mut_df= analysisYeast(raw_vcf_file, fastq_ID, orfs_df)
            # fully_covered_file = os.path.join(sub_output, "fully_covered_plateORFs.csv")
            # fully_covered.to_csv(fully_covered_file, index=False)

            all_group_summary_file = os.path.join(sub_output, "all_summary_subsetORFs.csv")
            all_mut_file = os.path.join(sub_output, "all_mutation_subsetORFs.csv")
            all_summary_df.to_csv(all_group_summary_file, index=False)            # all_found.to_csv(all_found_file, index=False)
            # fully_covered.to_csv(all_summary, index=False, header=False, mode="a")
            # all_found.to_csv(all_found_summary, index=False, header=False, mode="a")

            db = all_summary_df["db"].unique()
            stats_list.append("subsetORFs")
            genes_found.append(stats_list)
            mut_df["plate"] = fastq_ID
            mut_df["db"] = db[0]
            all_mut_df_subset.append(mut_df)
            all_summary_df["db"] = db[0]
            all_summary_subset.append(all_summary_df)

    # process all log
    all_log = pd.DataFrame(all_log)
    all_log_file = os.path.join(output, "alignment_log.csv")
    all_log.to_csv(all_log_file, index=False)

    # empty df to save all the variants
    all_summary_file_plate = os.path.join(output, "all_summary_plateORFs.csv")
    all_summary_file_all = os.path.join(output, "all_summary_allORFs.csv")
    all_summary_file_subset = os.path.join(output, "all_summary_subsetORFs.csv")

    all_mut_file_plate = os.path.join(output, "all_mutation_plateORFs.csv")
    all_mut_file_all = os.path.join(output, "all_mutation_allORFs.csv")
    all_mut_file_subset = os.path.join(output, "all_mutation_subsetORFs.csv")

    # process all summary and all mutation (plate)
    all_summary_df = pd.concat(all_summary_plate)
    all_summary_df.to_csv(all_summary_file_plate, index=False)
    # get all the mutations
    all_mut_df = pd.concat(all_mut_df_plate)
    # save to file
    #all_mut_file = os.path.join(output, "all_mutations.csv")
    all_mut_df.to_csv(all_mut_file_plate, index=False)

    # process all summary and all mutation (all)
    all_summary_df = pd.concat(all_summary_all)
    all_summary_df.to_csv(all_summary_file_all, index=False)
    # get all the mutations
    all_mut_df = pd.concat(all_mut_df_all)
    # save to file
    #all_mut_file = os.path.join(output, "all_mutations.csv")
    all_mut_df.to_csv(all_mut_file_all, index=False)

    # process all summary and all mutation (subset)
    all_summary_df = pd.concat(all_summary_subset)
    all_summary_df.to_csv(all_summary_file_subset, index=False)
    # get all the mutations
    all_mut_df = pd.concat(all_mut_df_subset)
    # save to file
    #all_mut_file = os.path.join(output, "all_mutations.csv")
    all_mut_df.to_csv(all_mut_file_subset, index=False)

    # process summary of number of genes found in each sample
    all_genes_stats = pd.DataFrame(genes_found, columns=["plate", "fully_aligned", "all_genes_found",
                                                         "all_targeted_on_plate", "all_targeted_full",
                                                         "n_fully_aligned_genes_with_any_mut", "n_ref", "aligned_to"])
    genes_found_file = os.path.join(output, "genes_stats.csv")

    all_genes_stats["% on plate fully aligned"] = all_genes_stats["all_targeted_full"] / all_genes_stats[
        "all_targeted_on_plate"]
    all_genes_stats.to_csv(genes_found_file, index=False)


def analysisYeast(raw_vcf_file, fastq_ID, orfs_df):
    """
    Run yeast variants analysis, make files for each plate and save to corresponded dir
    also return dfs for combining
    """
    analysis = ppsAnalysis.yeast_variant_analysis.yeastAnalysis(raw_vcf_file, fastq_ID, orfs_df)
    summary = analysis.get_full_cover()

    # all the genes with full coverage
    n_fully_aligned = summary[summary["fully_covered"] == "y"].shape[0]
    # all genes in ref fasta
    n_ref = summary.shape[0]
    # all genes found in this fastq file
    n_all_found = summary[summary["found"] == "y"].shape[0]

    # merge with target orfs
    merged_df = pd.merge(orfs_df, summary, how="left", left_on="orf_name", right_on="gene_ID")

    n_targeted = orfs_df.shape[0]
    n_targeted_full = merged_df[merged_df["fully_covered"] == "y"].shape[0]
    # split gene ID col
    # fully_covered["gene_ID"] = fully_covered["gene_ID"].str.replace(to_replace="-[A-G]", "A")
    # save all the genes that are found to output
    # save all the genes that are fully covered to the output folder
    # merge with target orfs
    # merged_df = pd.merge(orfs_df, all_found, how="left", left_on="orf_name", right_on="gene_ID")
    # merged_df = merged_df[~merged_df["gene_ID"].isnull()]
    merged_df["db"] = merged_df["gene_ID"].str.extract(r".*-([A-Z]+)-[1-9]")
    merged_df["count"] = merged_df["gene_ID"].str.extract(r".*-[A-Z]+-([1-9])")
    merged_df["gene_name"] = merged_df["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")

    merged_df = merged_df[["orf_name", "ORF_NAME_NODASH", "SYMBOL", "len(seq)", "plate", "db", "gene_name", "fully_covered", "found", "gene_len_mapped", "aligned_perc"]]
    # # fully_covered = fully_covered.replace(to_replace ='-index[0-9]+', value = '', regex = True)
    # fully_covered["db"] = fully_covered["gene_ID"].str.extract(r".*-([A-Z]+)-[1-9]")
    # fully_covered["count"] = fully_covered["gene_ID"].str.extract(r".*-[A-Z]+-([1-9])")
    # fully_covered["gene_name"] = fully_covered["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
    #
    # # merge with target orfs
    # merged_df_full = pd.merge(orfs_df, fully_covered.drop(['db'], axis=1), how="left", left_on="orf_name",
    #                      right_on="gene_ID")
    # merged_df_full = merged_df_full[~merged_df_full["gene_ID"].isnull()]
    # merged_df_full = merged_df_full[["orf_name", "ORF_NAME_NODASH", "SYMBOL", "len(seq)", "plate", "db", "gene_name"]]
    # merged_file = os.path.join(sub_output, "merged_with_targets.csv")
    # merged_df.to_csv(merged_file, index=False)
    # merged_df.to_csv(all_summary, mode="a", index=False, header=False)

    # filter vcf based on QUAL and DP
    mut_count_df = analysis.filter_vcf()
    if not mut_count_df.empty:
        # label mutations with syn/non-syn
        # load all sequences
        all_seq = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/all_sequence.csv"
        all_seq_df = pd.read_csv(all_seq)
        processed_mut = analysis.process_mut(all_seq_df, mut_count_df)
        # from fully aligned genes, select those with any mutations
        fully_aligned_with_mut = pd.merge(merged_df[["orf_name", "gene_name", "found", "fully_covered"]],
                                          processed_mut,
                                          how="left",
                                          left_on="orf_name",
                                          right_on="orf_name")
        mut_count_df = fully_aligned_with_mut[~fully_aligned_with_mut["ref"].isnull()]
        n_mut_genes_full = fully_aligned_with_mut[~fully_aligned_with_mut["ref"].isnull()]
        n_mut_genes_full = n_mut_genes_full["orf_name"].unique().shape[0]

    else:
        mut_count_df = pd.DataFrame({}, ["orf_name", "pos", "ref", "alt", "qual", "read_counts", "read_depth", "label",
                                         "type"])
        n_mut_genes_full = 0

    stats_list = [fastq_ID, n_fully_aligned, n_all_found, n_targeted, n_targeted_full, n_mut_genes_full, n_ref]

    return merged_df, stats_list, mut_count_df


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


    HIP_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/HIP_targeted_ORFs.csv"
    other_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/other_targeted_ORFs.csv"
    orfs = read_yeast_csv(HIP_target_ORFs, other_target_ORFs)

    return orfs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('--align', action="store_true", help='provide this argument if users want to start with '
                                                         'alignment, otherwise the program assumes alignment was '
                                                         'done and will analyze the vcf files.')
    parser.add_argument("-f", "--fastq", help="input fastq files", default="/home/rothlab/rli/01_ngsdata/yeast_PPS_seqs/merged_r1r2/")
    parser.add_argument("-r", "--ref", help="Path to referece files", default="/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/yeast_ref_all/")
    parser.add_argument('-o', "--output", help='Output directory', default="/home/rothlab/rli/02_dev/06_pps_pipeline/output/")
    parser.add_argument('-n', "--name", help='Name for this run', default="yeast_alignedwithmoreinfo")
    args = parser.parse_args()

    variants_main(args)
