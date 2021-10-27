#!/usr/bin/env python3.7
import pandas as pd
# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

## process files for yeast and human to make the final input file for PPS analysis
# the final input file can be found here:


def read_yeast_orf(HIP_target_ORFs, other_target_ORFs):
    """
    Join HIP data and other data into one df, remove unwanted columns. Save the merged df to file
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
    
    output_file = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/yeast_summary.csv"
    combined.to_csv(output_file)
    return combined


def read_human_orf(human_ref_with_seq, human91_ref):
    """
    human ref with enst/ensg ID
    """
    ref_df_91 = pd.read_csv(human91_ref)
    ref_df_ensembl = pd.read_csv(human_ref_with_seq)
    ref_df_91 = ref_df_91.fillna(-1)

    # merge this two df together
    # check if there are NAs in entrez gene ID and entrez gene symbol
    merged_df = pd.merge(ref_df_91, ref_df_ensembl, left_on=["entrez_gene_id", "entrez_gene_symbol"],
                         right_on=["entrez_gene_id", "symbol"], how="left")
    merged_df["grch37_filled"] = merged_df["cds_seq37"].fillna(merged_df["cds_seq"])
    merged_df["grch38_filled"] = merged_df["cds_seq38"].fillna(merged_df["cds_seq"])

    merged_df["entrez_gene_id"] = merged_df["entrez_gene_id"].astype(int)
    merged_df['orf_name'] = merged_df['orf_id'].astype(str) + "_" + merged_df['entrez_gene_id'].astype(str) + "_G0" + merged_df['Pool group #'].astype(str) + "_" + merged_df['entrez_gene_symbol'].astype(str)

    # humanallORF = pd.read_csv(human_ref)
    # humanallORF = humanallORF[["ORFID", "ensembl_transcript_id", "ensembl_protein_id", "ensembl_gene_id", "uniprot_AC_iso", "symbol", "entrez_gene_id", "CDS"]]
    #
    # humanallORF["entrez_gene_id"] = humanallORF["entrez_gene_id"].astype(int)
    # humanallORF['orf_name'] = humanallORF['entrez_gene_id'].astype(str) + "_" + humanallORF['entrez_gene_symbol'].astype(str)
    output_file = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/human_summary.csv"
    merged_df.to_csv(output_file)
    return merged_df

if __name__ == '__main__':
    # process yeast file
    HIP_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/HIP_targeted_ORFs.csv"
    other_target_ORFs = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/other_targeted_ORFs.csv"
    read_yeast_orf(HIP_target_ORFs, other_target_ORFs)

    # process human file
    ref_91 = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/20161117_ORFeome91_seqs.csv"
    ref_ensembl = "/home/rothlab/rli/02_dev/06_pps_pipeline/publicdb/merged_ensembl_sequence.csv"
    read_human_orf(ref_ensembl, ref_91)
