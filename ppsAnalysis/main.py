#!/usr/bin/env python3.7

"""
Main script for plasmid pool sequencing analysis

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

"""
import glob

import pandas as pd
import os
import argparse

import ppsAnalysis.alignment
import ppsAnalysis.cluster
import ppsAnalysis.yeast_variant_analysis
import logging.config


def write_full_cover(plate_name, all_genes, full_cover_genes, snp, indel, ref_dict, output_file):
    """

    :param plate_name: str
    :param all_genes: dictionary
    :param full_cover_genes:
    :param snp:
    :param indel:
    :param ref_dict:
    :param output_file:
    :return:
    """
    with open(output_file, "a") as output_file:
        output_file.write("plate_name,gene_name,gene_length,aligned_length,alignment_rate,total_read_count,average_read_depth,number_of_SNP,number_of_INDEL\n")
        for gene in all_genes.keys():
            if gene in full_cover_genes:
                all_genes = full_cover_genes
            if gene in snp.keys():

                line = [plate_name,
                        gene,
                        str(ref_dict[gene]),
                        str(all_genes[gene][0]),
                        str(all_genes[gene][0]/ref_dict[gene]),
                        str(all_genes[gene][1]),
                        str(all_genes[gene][2]),
                        str(len(snp[gene])),
                        "0"]
            else:
                line = [plate_name,
                        gene,
                        str(ref_dict[gene]),
                        str(all_genes[gene][0]),
                        str(all_genes[gene][0]/ref_dict[gene]),
                        str(all_genes[gene][1]),
                        str(all_genes[gene][2]),
                        "0",
                        "0"]
            if gene in indel.keys():
                line[-1] = str(indel[gene])
            output_file.write(",".join(line)+"\n")


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
                ref = arguments.ref + "orf9-1_" + fastq_ID
            elif arguments.mode == "yeast":
                fastq_ID = f.split("_")[0]
                # fastq_ID = r1.split("-")[-2]
                ref = arguments.ref + "ORF_withpDONR"
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
            alignment_obj = ppsAnalysis.alignment.Alignment(ref, f, sub_output, sh_file, align_log)
            # the main function writes to the sh file ans submit the file to cluster
            # return job ID
            at = 6
            job_id = alignment_obj._main(at)
            all_alignment_jobs.append(job_id)
        # track all alignment jobs
        alignment_log = logging.getLogger("alignment.log")
        jobs_finished = ppsAnalysis.cluster.parse_jobs_galen(all_alignment_jobs, alignment_log)
        if jobs_finished:
            main_logger.info("Alignment jobs all finished")
    # for each sample, parse vcf files
    all_log = []
    genes_found = []
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
        raw_vcf_file = os.path.join(sub_output, f"{fastq_ID}_L001_raw.vcf")
        # there should be only one log file in the dir
        log_file = glob.glob(f"{sub_output}/*.log")[0]
        if not os.path.isfile(raw_vcf_file):
            main_logger.warning(f"VCF file does not exist: {raw_vcf_file}")
            continue
        if not os.path.isfile(log_file):
            main_logger.warning(f"log file does not exist: {log_file}")
            continue
        # get information from the log file to make a summary log file for all the samples
        with open(log_file, "r") as log_f:
            for line in log_f:
                if "reads;" in line:
                    n_reads = line.split(" ")[0]
                if "alignment rate" in line:
                    perc_aligned = line.split("%")[0]
            all_log.append([fastq_ID, n_reads, perc_aligned])

        # for each vcf file, get how many genes are fully aligned
        if arguments.mode == "human":
            # extract ID
            pass
        else:  # yeast
            # first get the genes that are fully covered in the fastq files
            orfs_df = orfs[orfs["plate"] == fastq_ID]
          
            analysisYeast = ppsAnalysis.yeast_variant_analysis.yeastAnalysis(raw_vcf_file, fastq_ID, orfs_df)
            full_cover_genes, gene_dict, ref_dict = analysisYeast.get_full_cover()

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
            fully_covered["gene_ID"] = fully_covered["gene_ID"].str.replace("gene", "")
            fully_covered = fully_covered.replace(to_replace ='-index[0-9]+', value = '', regex = True)
            print(fully_covered["gene_ID"]) 
            try: 
                fully_covered[["gene_ID", "db"]] = fully_covered["gene_ID"].str.split("_", expand=True)
            except:
                fully_covered["db"] = "smorf"
            fully_covered_file = os.path.join(sub_output, "fully_covered.csv")
            fully_covered.to_csv(fully_covered_file, index=False)

            # select from fully_covered based on which db the orf is from
            if "HIP" in fastq_ID:
                filtered_db = fully_covered[fully_covered["db"] == "HIP"]
            else:
                filtered_db = fully_covered[fully_covered["db"] != "HIP"]

            # merge with target orfs
            merged_df = pd.merge(orfs_df, filtered_db, how="left", left_on="orf_name", right_on="gene_ID")
            merged_file = os.path.join(sub_output, "merged_with_targets.csv")
            merged_df.to_csv(merged_file, index=False)
            genes_found.append([fastq_ID, n_fully_aligned, n_all_found, n_ref])

            # merge fully covered gene to the targeted genes in this sample 

    
    # process all log
    all_log_df = pd.DataFrame(all_log, columns=["sample", "total reads", "alignment rate"])
    all_log_file = os.path.join(output, "alignment_log.csv")
    all_log_df.to_csv(all_log_file, index=False)

    # process summary of number of genes found in each sample
    all_genes_stats = pd.DataFrame(genes_found, columns=["sample", "fully_aligned", "all_genes_found", "n_ref"])
    genes_found_file = os.path.join(output, "genes_stats.csv")
    all_genes_stats.to_csv(genes_found_file, index=False)


def read_yeast_csv(HIP_target_ORFs, other_target_ORFs):
    """
    Join HIP data and other data into one df, remove unwanted columns
    :param HIP_target_ORFs: csv file contains which HIP ORF is in which sample
    :param other_target_ORFs: csv file contains which other ORF is in which sample
    :return: df with ORF name, db name and sample name
    """
    HIP_df = pd.read_csv(HIP_target_ORFs)
    other_target_ORFs = pd.read_csv(other_target_ORFs)

    HIP_df = HIP_df[["ORF_NAME_NODASH", "len(seq)", "SYMBOL", "plate"]]
    HIP_df["db"] = "HIP"
    HIP_df = HIP_df.rename(columns={"ORF_NAME_NODASH": "orf_name"})
    other_ORFs = other_target_ORFs[["orf_name", "src_collection", "plate"]]
    other_ORFs = other_ORFs.rename(columns={"src_collection": "db"})
    other_ORFs['plate'] = 'scORFeome-' + other_ORFs['plate'].astype(str)
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
