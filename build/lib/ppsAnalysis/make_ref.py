#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')
import os
import glob
import argparse
import pandas as pd

"""
Make reference fasta files for yeast or human, 
this step needs to be done before submitting jobs for 
alignment 
"""

def make_yeast_fasta(output):
    """

    :param output: output directory for fasta files
    :return:
    """
    HIP_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/HIP_targeted_ORFs.csv"
    other_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/other_targeted_ORFs.csv"

    hip_df = pd.read_csv(HIP_target_ORFs)
    other_df = pd.read_csv(other_target_ORFs)

    hip_df = hip_df[["ORF_id", "plate", "SEQ"]]
    # make fasta file for all HIP ORFs
    fasta_hip = os.path.join(output, "hip_all.fasta")
    with open(fasta_hip, "w") as output_hip:
        for index, row in hip_df.iterrows():
            seq_id = f"{row['ORF_id']}"
            output_hip.write(f">{seq_id}\n")
            output_hip.write(f"{row['SEQ']}\n")

    # make fasta file for plate specific ORFs
    all_plates = hip_df["plate"].unique().tolist()
    for p in all_plates:
        plate_hip = hip_df[hip_df["plate"] == p]
        plate_fasta = os.path.join(output, f"{p}.fasta")
        with open(plate_fasta, "w") as platefile:
            for index, row in plate_hip.iterrows():
                seq_id = f"{row['ORF_id']}"
                platefile.write(f">{seq_id}\n")
                platefile.write(f"{row['SEQ']}\n")

    all_sequence = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/all_sequence.txt"
    all_seq = pd.read_csv(all_sequence, sep="\t")
    PROTGEN = all_seq[all_seq["source"] == "PROTGEN"]
    SGD = all_seq[all_seq["source"] == "SGD"]

    # make fasta file for all sequences
    fasta_all = os.path.join(output, "all_seq.fasta")
    with open(fasta_all, "w") as all_output:
        for index, row in all_seq.iterrows():
            all_output.write(f">{row['ORF_id']}\n")
            all_output.write(f"{row['cds_seq']}\n")

    # make fasta file for all other protgen
    fasta_prot = os.path.join(output, "PROTGEN_all.fasta")
    with open(fasta_prot, "w") as output_prot:
        for index, row in PROTGEN.iterrows():
            output_prot.write(f">{row['ORF_id']}\n")
            output_prot.write(f"{row['cds_seq']}\n")

    # make fasta file for all other SGD
    fasta_prot = os.path.join(output, "SGD_all.fasta")
    with open(fasta_prot, "w") as output_prot:
        for index, row in SGD.iterrows():
            output_prot.write(f">{row['ORF_id']}\n")
            output_prot.write(f"{row['cds_seq']}\n")

    # sort SGD orfs and PROTGEN orfs into plate specific groups 
    # using other orfs df and all_sequence
    other = all_seq[(all_seq["source"] == "PROTGEN") | (all_seq["source"] == "SGD")]
    merged_seq_other = pd.merge(other_df, other, how="left", on="orf_name")
    other_plates = merged_seq_other["plate"].unique()
    for p in other_plates:
        plate_other = merged_seq_other[merged_seq_other["plate"] == p]
        plate_fasta = os.path.join(output, f"{p}.fasta")
        with open(plate_fasta, "w") as platefile:
            for index, row in plate_other.iterrows():
                seq_id = f"{row['ORF_id']}"
                platefile.write(f">{seq_id}\n")
                platefile.write(f"{row['cds_seq']}\n")

    # build bowtie2 index for later use 
    all_fasta = glob.glob(f"{output}/*.fasta")
    for f in all_fasta:
        f_id = os.path.basename(f).split(".")[0]
        cmd = f"bowtie2-build {f} {output}/{f_id}"
        os.system(cmd)
        
def main(mode, output):
    """
    Make fasta files in output dir
    :param mode: yeast or human
    :param output: output directory for fasta files
    :return: None
    """
    if mode == "human":
        pass
    else:
        make_yeast_fasta(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('-m', help='human or yeast')
    parser.add_argument("-o", help="output dir for fasta files")

    args = parser.parse_args()
    main(args.m, args.o)

