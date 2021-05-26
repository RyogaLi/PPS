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
import ppsAnalysis.human_variant_analysis
import ppsAnalysis.logthis
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted


def variants_main(arguments):
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
    file_list = os.listdir(arguments.fastq)
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
            fastq_ID = f.split(".")[0]
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
        cluster_log = log_obj.get_logger("alignment.log")
        jobs_finished = ppsAnalysis.cluster.parse_jobs_galen(all_alignment_jobs, cluster_log)
        if jobs_finished:
            main_logger.info("Alignment jobs all finished")

    parse_vcf_files_human(output, file_list, arguments, orfs, main_logger)

def parse_vcf_files_human(output, file_list, arguments, orfs, logger):
    # for each sample, parse vcf files
    all_log = {"fastq_ID": [], "reads": [], "map_perc": []}
    genes_found = []
    all_genes_summary = pd.DataFrame([],
                                     columns=["orf_id", 'entrez_gene_id', 'Pool group #', 'entrez_gene_symbol', 'Mapped reads', 'Verified', '# mut', 'orf_name', 'gene_ID', 'gene_len'])
    all_found_summary = os.path.join(output, "all_found_summary.csv")
    all_genes_summary.to_csv(all_found_summary, index=False)

    all_full_summary = os.path.join(output, "all_full_summary.csv")
    all_genes_summary.to_csv(all_full_summary, index=False)

    all_mut_df = []
    for f in file_list:
        if not f.endswith(".fastq.gz"): continue
        fastq_ID = f.split(".")[0]
        sub_output = os.path.join(os.path.abspath(output), fastq_ID)

        # there should be only one log file in the dir
        log_file = glob.glob(f"{sub_output}/*.log")[0]
        if not os.path.isfile(log_file):
            logger.warning(f"log file does not exist: {log_file}")
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
        all_log["fastq_ID"] += [fastq_ID] * 2

        # for each vcf file, get how many genes are fully aligned
        # get only the subset that are in the group 
        group_ID = fastq_ID.split("_")[-1][-1]
        orfs_df = orfs[orfs["Pool group #"] == int(group_ID)]
        raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_group_spec_orfs_raw.vcf")
        if os.path.isfile(raw_vcf_file):
            # analysis of ORFs aligned to group specific reference
            all_found, fully_covered, stats_list, mut_df = analysisHuman(raw_vcf_file, fastq_ID, orfs_df)
            fully_covered_file = os.path.join(sub_output, "fully_covered_groupSpecORFs.csv")
            fully_covered.to_csv(fully_covered_file, index=False)
            all_found_file = os.path.join(sub_output, "all_found_groupSpecORFs.csv")
            all_found.to_csv(all_found_file, index=False)

            fully_covered.to_csv(all_full_summary, index=False, header=False, mode="a")
            all_found.to_csv(all_found_summary, index=False, header=False, mode="a")

            stats_list.append("groupSpecORFs")
            genes_found.append(stats_list)
            mut_df["sample"] = fastq_ID
            all_mut_df.append(mut_df)

    # process all log
    all_log = pd.DataFrame(all_log)
    all_log_file = os.path.join(output, "alignment_log.csv")
    all_log.to_csv(all_log_file, index=False)

    # get all the mutations
    all_mut_df = pd.concat(all_mut_df)
    # save to file
    all_mut_file = os.path.join(output, "all_mutations.csv")
    all_mut_df.to_csv(all_mut_file, index=False)
    # process summary of number of genes found in each sample
    all_genes_stats = pd.DataFrame(genes_found, columns=["plate", "fully_aligned", "all_genes_found",
                                                         "all_targeted_on_plate", "all_targeted_full",
                                                         "n_fully_aligned_genes_with_any_mut", "n_ref", "aligned_to"])
    genes_found_file = os.path.join(output, "genes_stats.csv")

    all_genes_stats["% in group fully aligned"] = all_genes_stats["all_targeted_full"] / all_genes_stats[
        "all_targeted_on_plate"]
    all_genes_stats.to_csv(genes_found_file, index=False)
        #exit()


def read_human_csv(human91_ORFs):
    humanallORF = pd.read_csv(human91_ORFs)
    humanallORF = humanallORF[['orf_id', 'entrez_gene_id', 'Pool group #', 'entrez_gene_symbol', 'Mapped reads', 'Verified',
                                '# mut']]
    humanallORF = humanallORF.fillna(-1)
    humanallORF["entrez_gene_id"] = humanallORF["entrez_gene_id"].astype(int)
    humanallORF['orf_name'] = humanallORF['orf_id'].astype(str) + "_" + humanallORF['entrez_gene_id'].astype(str) + "_G0" + humanallORF['Pool group #'].astype(str) + "_" + humanallORF['entrez_gene_symbol'].astype(str)  

    return humanallORF


def analysisHuman(raw_vcf_file, fastq_ID, orfs_df):
    """

    """
    analysis = ppsAnalysis.human_variant_analysis.humanAnalysis(raw_vcf_file, fastq_ID, orfs_df)
    full_cover_genes, gene_dict, ref_dict = analysis.get_full_cover()

    # all the genes with full coverage
    n_fully_aligned = len(full_cover_genes.keys())
    # all genes in ref fasta
    n_ref = len(ref_dict.keys())
    # all genes found in this fastq file
    n_all_found = len(gene_dict.keys())

    # save all the genes that are fully covered to the output folder
    fully_covered = pd.DataFrame.from_dict(full_cover_genes, orient='index').reset_index()
    fully_covered.columns = ["gene_ID", "gene_len"]

    # save all the genes that are found to output
    # save all the genes that are fully covered to the output folder
    all_found = pd.DataFrame.from_dict(gene_dict, orient='index').reset_index()
    all_found.columns = ["gene_ID", "gene_len"]

    # merge with target orfs
    merged_df = pd.merge(orfs_df, all_found, how="left", left_on="orf_name", right_on="gene_ID")
    merged_df = merged_df[~merged_df["gene_ID"].isnull()]

    # merge with target orfs
    merged_df_full = pd.merge(orfs_df, fully_covered, how="left", left_on="orf_name", right_on="gene_ID")
    merged_df_full = merged_df_full[~merged_df_full["gene_ID"].isnull()]
    n_targeted = orfs_df.shape[0]
    n_targeted_full = merged_df_full.shape[0]

    mut_df = analysis.filter_vcf()
    mut_df = pd.DataFrame(mut_df)
    mut_df.columns = ["gene_ID", "pos", "ref", "alt", "qual", "read_counts", "read_depth", "label"]

    # merge mut_df with fully covered
    merge_mut = pd.merge(mut_df, fully_covered, how="left", on="gene_ID")
    merge_mut_fully_covered = merge_mut[~merge_mut["gene_len"].isnull()]
    processed_mut = analysis.process_mut(mut_df)
    # from fully aligned genes, select those with any mutations
    fully_aligned_with_mut = pd.merge(fully_covered, processed_mut, how="left", left_on="gene_ID",
                                      right_on="gene_ID")
    mut_count_df = fully_aligned_with_mut[~fully_aligned_with_mut["ref"].isnull()]
    n_mut_genes_full = fully_aligned_with_mut[~fully_aligned_with_mut["ref"].isnull()]
    n_mut_genes_full = n_mut_genes_full["gene_ID"].unique().shape[0]
    # count how many ORFs have variants
    n_orf_with_v = len(merge_mut_fully_covered["gene_ID"].unique())

    # from fully aligned genes, select those with any mutations
    stats_list = [fastq_ID, n_fully_aligned, n_all_found, n_targeted, n_targeted_full, n_orf_with_v, n_ref]
    return merged_df, merged_df_full, stats_list, mut_count_df


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

    human_91ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/human_91/20161117_ORFeome91_seqs.csv"
    orfs = read_human_csv(human_91ORFs)

    return orfs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('--align', action="store_true", help='provide this argument if users want to start with '
                                                         'alignment, otherwise the program assumes alignment was '
                                                         'done and will analyze the vcf files.')
    parser.add_argument("-f", "--fastq", help="input fastq files", default="/home/rothlab/rli/01_ngsdata/PPS_data/orfPool/merged_pool9-1/")
    parser.add_argument("-r", "--ref", help="Path to referece files", default="/home/rothlab/rli/human_test0/")
    parser.add_argument('-o', "--output", help='Output directory', default="/home/rothlab/rli/02_dev/06_pps_pipeline/output/")
    parser.add_argument('-n', "--name", help='Name for this run', default="/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/human_91/")
    args = parser.parse_args()

    variants_main(args)
