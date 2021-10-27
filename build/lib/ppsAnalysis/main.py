#!/usr/bin/env python3.7
import os
import glob
import pandas as pd
import argparse
import seaborn as sns

import ppsAnalysis.alignment
import ppsAnalysis.cluster
import ppsAnalysis.yeast_variant_analysis
import ppsAnalysis.human_variant_analysis
import ppsAnalysis.logthis
# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca


def variants_main(arguments):
    """
    Main function for pps analysis for yeast data
    :param arguments: user input arguments
    :return: None
    """
    # set log level
    loglevel = arguments.log
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
            at = 9
            job_id = alignment_obj.main(at)
            all_alignment_jobs.append(job_id)
        # track all alignment jobs
        cluster_log = log_obj.get_logger("alignment.log")
        jobs_finished = ppsAnalysis.cluster.parse_jobs_galen(all_alignment_jobs, cluster_log)
        if jobs_finished:
            main_logger.info("Alignment jobs all finished")
    # start barcode counting
    if arguments.mode == "human":
        parse_vcf_files(output, file_list, arguments, orfs, main_logger)
    elif arguments.mode == "yeast":
        parse_vcf_files(output, file_list, arguments, orfs, main_logger)


def parse_vcf_files(output, file_list, arguments, orfs, logger):
    """
    Process each vcf file based on mode, call human/yeast analysis script separately
    :param output: output directory
    :param file_list: file list (fastq files)
    :param arguments: user input arguments
    :param logger:
    :param mode:
    :return:
    """
    # for each sample, parse vcf files
    all_log = {"fastq_ID": [], "reads": [], "map_perc": []}
    genes_found = []
    all_summary_file = os.path.join(output, "all_summary.csv")
    all_summary = []
    all_mut_df = []
    # go through files in the file list
    for f in file_list:
        if not f.endswith(".fastq.gz"): continue
        fastq_ID = f.split("_")[0]
        logger.info(f"Processing {fastq_ID}")
        sub_output = os.path.join(os.path.abspath(output), fastq_ID)

        # there should be only one log file in the dir
        # this log file contains sam alignment rate
        try:
            log_file = glob.glob(f"{sub_output}/*.log")[0]
        except:
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
        # depends on which data set we are processing
        if arguments.mode == "human":
            group_ID = fastq_ID.split("_")[-1][-1]
            # for human sequencing data, we process one group at a time
            orfs_df = orfs[orfs["Pool group #"] == int(group_ID)]
            raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_group_spec_orfs_raw.vcf")
            # if the file is not found, raise error
            if not os.path.isfile(raw_vcf_file):
                raise FileNotFoundError(f"{raw_vcf_file}")
            # analysis of ORFs aligned to group specific reference
            all_summary_df, stats_list, mut_df = analysisHuman(raw_vcf_file, fastq_ID, orfs_df, sub_output, arguments.refName, logger)

            all_group_summary_file = os.path.join(sub_output, "all_summary_groupSpecORFs.csv")
            all_summary_df.to_csv(all_group_summary_file, index=False)

            stats_list.append("groupSpecORFs")
            genes_found.append(stats_list)
            mut_df["sample"] = fastq_ID
            all_summary_df["sample"] = fastq_ID
            all_mut_df.append(mut_df)
            all_summary.append(all_summary_df)

        elif arguments.mode == "yeast":
            # for each vcf file, get how many genes are fully aligned
            # for yeast, we aligned all the orfs to ORFs on each plate
            orfs_df = orfs[orfs["plate"] == fastq_ID]
            raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_plateORFs_raw.vcf")

            # analysis of ORFs aligned to subgroup
            all_summary_df, stats_list, mut_df= analysisYeast(raw_vcf_file, fastq_ID, orfs_df, logger)
            all_group_summary_file = os.path.join(sub_output, "all_summary_plateORFs.csv")
            all_summary_df.to_csv(all_group_summary_file, index=False)

            db = all_summary_df["db"].unique()
            stats_list.append("plateORFs")
            genes_found.append(stats_list)
            mut_df["plate"] = fastq_ID
            mut_df["db"] = db[0]
            all_mut_df.append(mut_df)
            all_summary_df["db"] = db[0]
            all_summary.append(all_summary_df)

    # process all summary
    all_summary_df = pd.concat(all_summary)
    all_summary_df = all_summary_df.reset_index(drop=True)
    all_summary_df.to_csv(all_summary_file, index=False)
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

    all_genes_stats["% on plate fully aligned"] = all_genes_stats["all_targeted_full"] / all_genes_stats[
        "all_targeted_on_plate"]
    all_genes_stats.to_csv(genes_found_file, index=False)


def analysisHuman(raw_vcf_file, fastq_ID, orfs_df, suboutput, ref, logger):
    """

    :param raw_vcf_file: input vcf file
    :param fastq_ID: ID for the fastq file
    :param orfs_df: dataframe contains all input orfs (group specific)
    :param suboutput: sub output folder
    :param ref:
    :return:
    """
    analysis = ppsAnalysis.human_variant_analysis.humanAnalysis(raw_vcf_file, fastq_ID, orfs_df, ref)
    logger.info("Getting fully covered stats")
    summary = analysis.get_full_cover()
    # all the genes with full coverage
    n_fully_aligned = summary[summary["fully_covered"] == "y"].shape[0]
    # all genes in ref fasta
    n_ref = summary.shape[0]
    logger.info(f"Total targeted ORFs in this group: {n_ref}")
    # all genes found in this fastq file
    n_all_found = summary[summary["found"] == "y"].shape[0]
    logger.info(f"Total ORFs found in this group: {n_all_found}")
    # merge with target orfs
    merged_df = pd.merge(orfs_df, summary, how="left", left_on="orf_name", right_on="gene_ID")
    # total number of targeted ORFs in this group
    n_targeted = orfs_df.shape[0]
    # total number of targeted ORFs in this group that are fully covered
    n_targeted_full = merged_df[merged_df["fully_covered"] == "y"].shape[0]
    # filter vcf file to remove variants that didn't pass filter
    # return df with columns: ["gene_ID", "pos", "ref", "alt", "qual", "read_counts", "read_depth", "label"]
    mut_df = analysis.filter_vcf()
    # merge mut_df with fully covered
    merge_mut = pd.merge(mut_df, merged_df, how="left", on="gene_ID")
    merge_mut_fully_covered = merge_mut[merge_mut["fully_covered"] == "y"]
    mut_file = os.path.join(suboutput, "all_mut.csv")
    # only process mut if the mut file for this subgroup doesn't exist, to save time
    # because process_mut gets data from gnomAD

    if not os.path.isfile(mut_file) or os.stat(mut_file).st_size == 0:
        logger.info("Process mutations")
        processed_mut = analysis._process_mut(mut_df)
        processed_mut.to_csv(mut_file)
    else:
        processed_mut = pd.read_csv(mut_file)
        processed_mut = processed_mut.drop(processed_mut.columns[0], axis=1)

    # from fully aligned genes, select those with any mutations
    fully_aligned_with_mut = pd.merge(merged_df[["gene_ID", "entrez_gene_symbol", "found", "fully_covered", "gene_len", "gene_len_mapped", "aligned_perc"]],
                                      processed_mut,
                                      how="left",
                                      on="gene_ID")
    mut_count_df = fully_aligned_with_mut[~fully_aligned_with_mut["ref"].isnull()]
    # n_mut_genes_full = fully_aligned_with_mut[~fully_aligned_with_mut["ref"].isnull()]["gene_ID"].unique().shape[0]
    # count how many ORFs have variants
    n_orf_with_v = len(merge_mut_fully_covered["gene_ID"].unique())

    # from fully aligned genes, select those with any mutations
    stats_list = [fastq_ID, n_fully_aligned, n_all_found, n_targeted, n_targeted_full, n_orf_with_v, n_ref]
    return merged_df, stats_list, mut_count_df


def analysisYeast(raw_vcf_file, fastq_ID, orfs_df, logger):
    """
    Run yeast variants analysis, make files for each plate and save to corresponded dir
    also return dfs for combining
    """
    analysis = ppsAnalysis.yeast_variant_analysis.yeastAnalysis(raw_vcf_file, fastq_ID, orfs_df)
    logger.info("Getting fully covered ORFs..")
    summary = analysis.get_full_cover()

    # all the genes with full coverage
    n_fully_aligned = summary[summary["fully_covered"] == "y"].shape[0]
    # all genes in ref fasta
    n_ref = summary.shape[0]
    logger.info(f"Total targeted ORFs on this plate: {n_ref}")
    
    # all genes found in this fastq file
    n_all_found = summary[summary["found"] == "y"].shape[0]
    logger.info(f"Total ORFs found on this plate: {n_all_found}")
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
        #all_seq = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/all_sequence.csv"
        #all_seq_df = pd.read_csv(all_seq)
        #print(all_seq_df)
        processed_mut = analysis.process_mut(orfs_df, mut_count_df)
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


def check_args(arguments):
    """
    Check user input arguments
    :param arguments: user input arguments, argparse struct
    :return: None
    """

    if arguments.mode == "human":
        # if we are running analysis for human
        # it is required to also specify which reference to align to
        if arguments.refName == "":
            raise ValueError("Please also specify --refName grch37, grch38 or cds_seq")
        # orfs = read_human91(human_91ORFs)
        # orfs = pd.read_csv(arguments.summaryFile)
    elif arguments.mode == "yeast":
        pass
    else:
        raise ValueError("Wrong mode, please select human or yeast (case sensitive)")

    if not os.path.isfile(arguments.summaryFile):
        raise ValueError("Please provide summary file")

    orfs = pd.read_csv(arguments.summaryFile)

    if not os.path.isdir(arguments.output):
        raise NotADirectoryError(f"{arguments.output} does not exist")

    if arguments.fastq and not os.path.isdir(arguments.fastq):
        raise NotADirectoryError(f"{arguments.fastq} does not exist")

    if arguments.align:
        # if alignment, must also provide path to fastq files and reference files
        if not arguments.fastq or not arguments.ref:
            raise ValueError("Please also provide reference dir and fastq dir")

    return orfs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('--align', action="store_true", help='provide this argument if users want to start with '
                                                         'alignment, otherwise the program assumes alignment was '
                                                         'done and will analyze the vcf files.')
    parser.add_argument("-f", "--fastq", help="input fastq files", required=True)
    parser.add_argument("-n", "--name", help="Run name")
    parser.add_argument("-m", "--mode", help="Yeast or Human", required=True)
    parser.add_argument('-o', "--output", help='Output directory', required=True)
    parser.add_argument('-r', "--ref", help='Path to reference', required=True)
    parser.add_argument("--refName", help="grch37, grch38, human91")
    parser.add_argument("--summaryFile", help="Summary file contains ORF information")
    #parser.add_argument("--orfseq", help="File contains ORF sequences")
    parser.add_argument("--log", help="set log level", default="info")
    args = parser.parse_args()

    variants_main(args)
