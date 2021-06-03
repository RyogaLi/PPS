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


def make_human_fasta(output):
    """
    Make fasta files for human 9.1
    One fasta for all the sequences
    Also group specific sequences
    :param ref_91: csv file contains human 9.1 reference sequences
    :return: None
    """
    ref_91 = "/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/human_91/20161117_ORFeome91_seqs.csv"
    ref_df_91 = pd.read_csv(ref_91)
    ref_df_91 = ref_df_91.fillna(-1)
    # make all ref fasta
    all_ref_output = os.path.join(output, "all_ref_human.fasta")
    with open(all_ref_output, "w") as all_ref:
        for index, row in ref_df_91.iterrows():
            id_line = f">{row['orf_id']}_{int(row['entrez_gene_id'])}_G0{row['Pool group #']}_{row['entrez_gene_symbol']}\n"
            seq = row["cds_seq"]+"\n"
            all_ref.write(id_line)
            all_ref.write(seq)

    # make group sepecific fasta
    # get all groups
    groups = ref_df_91["Pool group #"].unique().tolist()
    for g in groups:
        # make fasta
        group_fasta = os.path.join(output, f"group_ref_G0{g}.fasta")
        # select subset of orfs belongs to this group
        subset = ref_df_91[ref_df_91["Pool group #"] == g]
        with open(group_fasta, "w") as g_fasta:
            for index, row in subset.iterrows():
                id_line = f">{row['orf_id']}_{int(row['entrez_gene_id'])}_G0{row['Pool group #']}_{row['entrez_gene_symbol']}\n"
                seq = row["cds_seq"] + "\n"
                g_fasta.write(id_line)
                g_fasta.write(seq)

    # build bowtie2 index for later use 
    all_fasta = glob.glob(f"{output}/*.fasta")
    for f in all_fasta:
        f_id = os.path.basename(f).split(".")[0]
        cmd = f"bowtie2-build {f} {output}/{f_id}"
        os.system(cmd)


def make_human_fasta_ensembl(output):
    """
    Make reference fasta file with ensembl ref seq
    :return:
    """
    ref_91 = "/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/human_91/20161117_ORFeome91_seqs.csv"
    ref_ensembl = "/home/rothlab/rli/02_dev/06_pps_pipeline/publicdb/merged_ensembl_sequence.csv"
    ref_df_91 = pd.read_csv(ref_91)
    ref_df_ensembl = pd.read_csv(ref_ensembl)
    ref_df_91 = ref_df_91.fillna(-1)
    print(ref_df_91.shape)
    # merge this two df together
    # check if there are NAs in entrez gene ID and entrez gene symbol
    print(ref_df_91[ref_df_91[["entrez_gene_id", "entrez_gene_symbol"]].duplicated()])
    ref_df_ensembl = ref_df_ensembl.drop_duplicates(subset=["entrez_gene_id", "symbol"])
    print(ref_df_ensembl.shape)
    merged_df = pd.merge(ref_df_91, ref_df_ensembl, left_on=["entrez_gene_id", "entrez_gene_symbol"], right_on=["entrez_gene_id", "symbol"], how="left")
    
    # make grch37 and 38 output dir if not exist
    grch37_output = os.path.join(output, "grch37")
    if not os.path.isdir(grch37_output):
        os.makedir(grch37_output)

    grch38_output = os.path.join(output, "grch38")
    if not os.path.isdir(grch38_output):
        os.makedir(grch38_output)

    # make group sepecific fasta
    # get all groups
    groups = merged_df["Pool group #"].unique().tolist()
    for g in groups:
        # make fasta for grch37
        # for missing values in cds_seq37, fill with original cds_seq
        merged_df["cds_seq37_filled"] = merged_df["cds_seq37"].fillna(merged_df["cds_seq"])
        group_fasta = os.path.join(grch37_output, f"group_ref_G0{g}.fasta")
        # select subset of orfs belongs to this group
        subset = merged_df[merged_df["Pool group #"] == g]
        with open(group_fasta, "w") as g_fasta:
            for index, row in subset.iterrows():
                id_line = f">{row['orf_id']}_{int(row['entrez_gene_id'])}_G0{row['Pool group #']}_{row['entrez_gene_symbol']}\n"
                # remove stop codon
                seq = row["cds_seq37_filled"][:-3] + "\n"
                g_fasta.write(id_line)
                g_fasta.write(seq)

        # make fasta for grch38
        # for missing values in cds_seq38, fill with original cds_seq
        merged_df["cds_seq38_filled"] = merged_df["cds_seq38"].fillna(merged_df["cds_seq"])
        group_fasta = os.path.join(grch38_output, f"group_ref_G0{g}.fasta")
        # select subset of orfs belongs to this group
        subset = merged_df[merged_df["Pool group #"] == g]
        with open(group_fasta, "w") as g_fasta:
            for index, row in subset.iterrows():
                id_line = f">{row['orf_id']}_{int(row['entrez_gene_id'])}_G0{row['Pool group #']}_{row['entrez_gene_symbol']}\n"
                # remove stop codon
                seq = row["cds_seq38_filled"][:-3] + "\n"
                g_fasta.write(id_line)
                g_fasta.write(seq)

    # build bowtie2 index for later use 
    all_fasta = glob.glob(f"{grch37_output}/*.fasta")
    for f in all_fasta:
        f_id = os.path.basename(f).split(".")[0]
        cmd = f"bowtie2-build {f} {grch37_output}/{f_id}"
        os.system(cmd)

    all_fasta = glob.glob(f"{grch38_output}/*.fasta")
    for f in all_fasta:
        f_id = os.path.basename(f).split(".")[0]
        cmd = f"bowtie2-build {f} {grch38_output}/{f_id}"
        os.system(cmd)


def compare_human_ref(output):
    """

    :param output:
    :return:
    """
    """
        Make reference fasta file with ensembl ref seq
        :return:
        """
    ref_91 = "/home/rothlab/rli/02_dev/06_pps_pipeline/fasta/human_91/20161117_ORFeome91_seqs.csv"
    ref_ensembl = "/home/rothlab/rli/02_dev/06_pps_pipeline/publicdb/merged_ensembl_sequence.csv"
    ref_df_91 = pd.read_csv(ref_91)
    ref_df_ensembl = pd.read_csv(ref_ensembl)
    ref_df_91 = ref_df_91.fillna(-1)
    print(ref_df_91.shape)
    # merge this two df together
    # check if there are NAs in entrez gene ID and entrez gene symbol
    print(ref_df_91[ref_df_91[["entrez_gene_id", "entrez_gene_symbol"]].duplicated()])
    ref_df_ensembl = ref_df_ensembl.drop_duplicates(subset=["entrez_gene_id", "symbol"])
    print(ref_df_ensembl.shape)
    merged_df = pd.merge(ref_df_91, ref_df_ensembl, left_on=["entrez_gene_id", "entrez_gene_symbol"],
                         right_on=["entrez_gene_id", "symbol"], how="left")


def main(mode, output):
    """
    Make fasta files in output dir
    :param mode: yeast or human
    :param output: output directory for fasta files
    :return: None
    """
    if mode == "human":
        make_human_fasta_ensembl(output)
    else:
        make_yeast_fasta(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis')
    parser.add_argument('-m', help='human or yeast')
    parser.add_argument("-o", help="output dir for fasta files")

    args = parser.parse_args()
    main(args.m, args.o)

