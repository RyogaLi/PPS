#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('--align', action="store_true", help='provide this argument if users want to start with '
                                                         'alignment, otherwise the program assumes alignment was '
                                                         'done and will analyze the vcf files.')
    parser.add_argument("-f", "--fastq", help="input fastq files", default="/home/rothlab/rli/01_ngsdata/PPS_data/orfPool/merged_pool9-1/")
    parser.add_argument("-n", "--name", help="Run name", default="pps")
    parser.add_argument('-o', "--output", help='Output directory', default="/home/rothlab/rli/02_dev/06_pps_pipeline/output/", required=True)
    parser.add_argument('-r', "--ref", help='Path to reference', default="/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/human_91/", required=True)
    parser.add_argument("--refName", help="grch37, grch38, cds_seq. Required if mode == human")
    parser.add_argument("--summaryFile", help = "Yeast or Human summary file")
    parser.add_argument("--orfseq", help="File contains ORF sequences")
    args = parser.parse_args()

    variants_main(args)
