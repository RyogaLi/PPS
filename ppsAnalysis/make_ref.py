#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')
import os
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

    hip_df = hip_df[["ORF_NAME", "plate", "SEQ"]]
    # make fasta file for all HIP ORFs
    fasta_hip = os.path.join(output, "hip_all.fasta")
    for index, row in hip_df.iterrows():
        with open(fasta_hip, "w") as output_hip:
            seq_id = f"{row['ORF_NAME']}-HIP"
            output_hip.write(f">{seq_id}\n")
            output_hip.write(f"{row['SEQ']}\n")

    # make fasta file for plate specific ORFs
    all_plates = hip_df["plate"].unique().tolist()
    for p in all_plates:
        plate_hip = hip_df[hip_df["plate"] == p]
        plate_fasta = os.path.join(output, f"{p}.fasta")
        for index, row in plate_hip:
            with open(plate_fasta, "w") as platefile:
                seq_id = f"{row['ORF_NAME']}-HIP"
                platefile.write(f">{seq_id}\n")
                platefile.write(f"{row['SEQ']}\n")

    all_sequence = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/all_sequence.txt"
    all_seq = pd.read_csv(all_sequence, sep="\t")
    PROTGEN = all_seq[all_seq["source"] == "PROTGEN"]
    SGD = all_seq[all_seq["source"] == "SGD"]

    # make fasta file for all other SGD
    for index, row in PROTGEN.iterrows():

    # make fasta file for all other PROTGEN



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


    main()

