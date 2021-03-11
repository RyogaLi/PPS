#!/usr/bin/env python3.7

"""
Main script for plasmid pool sequencing analysis

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

"""

import pandas as pd
import os
import argparse

import ppsAnalysis.alignment
import ppsAnalysis.cluster
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
    # create output folder with user input name
    run_name = arguments.n
    output = os.path.join(arguments.output, run_name)
    if not os.path.isdir(output):
        os.mkdir(output)

    logging.config.fileConfig("./ppsAnalysis/logging.conf")
    main_logger = logging.getLogger("main")

    if arguments.align:
        align_log = logging.getLogger("align.log")
        # first align fastq files if user want to use alignment
        file_list = os.listdir(arguments.fastq)
        all_alignment_jobs = []
        for f in file_list:
            if not f.endswith(".fastq.gz"): continue
            f = os.path.join(arguments.fastq, f)
            # for all the fastq files in the dir
            # align the fastq files with given reference
            # note that fastq files and the corresponding reference file has the same id
            if arguments.mode == "human":
                # extract ID
                fastq_ID = f.split("-")[-2]
                ref = arguments.ref + "orf9-1_" + fastq_ID
            elif arguments.mode == "yeast":
                fastq_ID = f.split(".")[0]
                # fastq_ID = r1.split("-")[-2]
                ref = arguments.ref + "ORF_combined_ref"
            else:
                raise ValueError("Please provide valid mode: human or yeast")
            # mkae sub_output dir for this sample
            sub_output = os.path.join(output, fastq_ID)
            if not os.path.isdir((sub_output)):
                os.mkdir(sub_output)
            # make sh file for submission in sub_output directory for alignment
            # this is developped for GALEN cluster
            sh_file = os.path.join(sub_output, f"{fastq_ID}.sh")
            alignment_obj = ppsAnalysis.alignment.Alignment(ref, f, sub_output, sh_file, align_log)
            # the main function writes to the sh file ans submit the file to cluster
            # return job ID
            at = 12
            job_id = alignment_obj._main(at)
            all_alignment_jobs.append(job_id)
        # track all alignment jobs
        alignment_log = logging.getLogger("alignment.log")
        jobs_finished = ppsAnalysis.cluster.parse_jobs_galen(all_alignment_jobs, alignment_log)
        if jobs_finished:
            main_logger.info("Alignment jobs all finished")

    # dir_list = os.listdir(output)
    # print(output)
    # for dir in dir_list:
    #     if not os.path.isdir(output + "/" + dir):
    #         continue
    #     os.chdir(os.path.join(output, dir))
    #     for file in os.listdir("."):
    #         # step 4
    #         # analysis
    #         # after variant calling
    #         # filter vcf file and get fully aligned genes with their avg read depth
    #         if file.endswith(".log") and os.stat(file).st_size == 0:
    #             os.remove(file)
    #         if file.endswith(".raw.vcf"):
    #             print(file)
    #             # full_cover is a dictionary contains key=gene ids for genes that are fully covered
    #             # value = [fully covered gene length, average read depth for this gene]
    #             logger.info("Analyzing %s ...", file)
    #             full_cover, total_gene_count, all_gene_dict, ref_dict = get_full_cover(file)
    #             logger.info("Got fully aligned genes in %s", file)
    #             snp, indel, read_depth = filter_vcf(file, all_gene_dict)
    #             logger.info("Filtered %s", file.replace(".raw.vcf", "_filtered.vcf"))
    #
    #             if clinvar:
    #                 VA = SnpAnalysis(orf_id_convert, clinvar_db)
    #                 all_clinvar_snp= VA._main(file.replace(".raw.vcf", "_filtered.vcf"))
    #                 all_clinvar_snp.drop_duplicates()
    #                 # write this information to file
    #                 all_clinvar_snp.to_csv(path_or_buf="./clinvar_variants.txt", sep="\t")
    #                 print(all_clinvar_snp)
    #             else:
    #                 out_file = output + "summary.csv"
    #                 write_full_cover(file, all_gene_dict, full_cover, snp, indel, ref_dict, out_file)
    #
    #             # take top 5 genes based on read depth
    #             # plot_top_n(snp, indel, read_depth, 3)
    #             # logger.info("Plot generated")
    #             logger.info("Analaysis done for %s \n", file)

        # break


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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('--align', action="store_true", help='provide this argument if users want to start with '
                                                         'alignment, otherwise the program assumes alignment was '
                                                         'done and will analyze the vcf files.')
    parser.add_argument("-f", "--fastq", help="input fastq files")
    parser.add_argument("-m", "--mode", help="Human or Yeast PPS?")
    parser.add_argument('-o', "--output", help='Output directory', required=True)
    parser.add_argument('-n', "--name", help='Name for this run', required=True)
    args = parser.parse_args()

    variants_main(args)