#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')

import os
import re
import pandas as pd


"""
Analysis of yeast PPS variants (vcf files) generated with samtools
"""

class yeastAnalysis(object):

    def __init__(self, input_vcf, basename, orfs_df):
        """
        Take vcf file from input dir
        :param input_: Input dir contains vcf files, sam files and bam files
        """
        self._rawvcf = input_vcf
        self._vvcf = input_vcf.replace("_raw", "_variants")
        self._basename = basename
        # contains all the targeted orfs for this sample
        self._orfs = orfs_df

    def get_full_cover(self):
        """
        Get a dictionary of gene names which are fully covered(aligned) in vcf file
        :return: dictionary with keys = gene names; value = gene length
        """
        with open(self._rawvcf, "r") as raw:
            gene_dict = {}
            ref_dict = {}
            temp_tracker = {}
            for line in raw:
                # in vcf header, grep gene names and gene len
                id_line = re.search("<ID=(.+?),length=(.+?)>", line)
                # assign gene length to geneID
                if id_line:
                    ref_dict[id_line.group(1)] = int(id_line.group(2))
                # not a header line
                if "#" not in line:
                    line = line.split()
                    if gene_dict.get(line[0], -1) == -1:
                        # grep read depth information from INFO section
                        # rd = re.search("DP=([0-9]+)", line[7])
                        # rd = rd.group(1)
                        gene_dict[line[0]] = 1
                        temp_tracker[line[0]] = [int(line[1])]
                    else:
                        # grep read depth information from INFO section
                        # rd = re.search("DP=([0-9]+)", line[7])
                        # rd = rd.group(1)
                        if int(line[1]) not in temp_tracker[line[0]]:
                            gene_dict[line[0]] += 1
                            temp_tracker[line[0]].append(int(line[1]))
                        # gene_dict[line[0]][1] += int(rd)
            removed_genes = gene_dict.copy()
            for key in gene_dict.keys():
                if gene_dict[key] < int(ref_dict[key]):
                    del removed_genes[key]
                # else:
                #     avg_rd = remove_genes[key][1] / remove_genes[key][0]
                #     remove_genes[key].append(avg_rd)
        # save all the genes that are fully covered to the output folder
        fully_covered = pd.DataFrame.from_dict(removed_genes, orient='index').reset_index()
        fully_covered.columns = ["gene_ID", "gene_len"]
        fully_covered["fully_covered"] = "y"

        # save all the genes that are found to output
        # save all the genes that are fully covered to the output folder
        all_found = pd.DataFrame.from_dict(gene_dict, orient='index').reset_index()
        all_found.columns = ["gene_ID", "gene_len_mapped"]
        all_found["found"] = "y"

        # save all the genes that are found to output
        # save all the genes that are fully covered to the output folder
        all_ref = pd.DataFrame.from_dict(ref_dict, orient='index').reset_index()
        all_ref.columns = ["gene_ID", "gene_len"]

        # join three dfs into one
        # with column names = ["gene_ID", "gene_len", "fully covered", "found"]
        # merge all found to all_ref
        summary = pd.merge(all_ref, all_found[["gene_ID", "gene_len_mapped", "found"]], how="outer", on="gene_ID")
        summary = pd.merge(summary, fully_covered[["gene_ID", "fully_covered"]], how="left", on="gene_ID")
        summary["aligned_perc"] = summary["gene_len_mapped"]/summary["gene_len"]
        return summary

    def filter_vcf(self):
        """
        for a given vcf (with variants only), filter the variants based on QUAL and DP
        write passed filter variants to a new vcf file 
        """
        
        filtered_vcf = self._vvcf.replace("_variants", "_filtered")
        mut_count = []
        with open(self._vvcf, "r") as rawvcf:
            with open(filtered_vcf, "w") as filteredvcf:
                for line in rawvcf:
                    if "#" in line: continue
                    l = line.split()
                    # get DP for this position 
                    info_title = l[-2].split(":")
                    info = l[-1].split(":")
                    info_dict = dict(zip(info_title, info))
                    if int(info_dict["DP"]) < 10:  # informative read depth
                        continue
                    # get variant call with quality > 20
                    try:
                        qual = float(l[5])
                    except:
                        continue
                    if qual < 20: continue
                    # if variant have two alt, split them and use the one with the most read counts
                    alt_bases = l[4].split(",")
                    alt_bases = [l[3]] + alt_bases
                    AD = info_dict["AD"].split(",")
                    alt_depth = dict(zip(alt_bases, AD))
                    df = pd.DataFrame(alt_depth.items())
                    df.columns = ["alt_base", "read_count"]
                    df["perc"] = df["read_count"].astype(float) / float(info_dict["DP"])
                    # select alt bases greater than 80%
                    df = df[df["perc"] > 0.8]
                    if df.empty:
                        continue
                    if l[3] in df["alt_base"].tolist():
                        continue
                    mut_base = df["alt_base"].tolist()[0]
                    mut_counts = df["read_count"].tolist()[0]
                    if len(l[3]) > 1:
                        label = "indel"
                    elif len(mut_base) > 1:
                        label = "indel"
                    else:
                        label = "SNP"
                    # track how many variants for each gene (with more than 10 reads mapped to it)
                    mut_count.append([l[0], l[1], l[3], mut_base, l[5], mut_counts, info_dict["DP"], label])
                    filteredvcf.write(line)
        mut_df = pd.DataFrame(mut_count)
        mut_cols = ["orf_name", "pos", "ref", "alt", "qual", "read_counts", "read_depth", "label"]

        if mut_df.empty:
            mut_df = pd.DataFrame({}, columns=mut_cols)
        else:
            mut_df.columns = mut_cols
        return mut_df
    
    def process_mut(self, all_df, mut_df):
        """
        Based on the coding sequence and protein sequence, determine if the variant is a syn/non-syn variant
        :param all_df: data frame contains all the coding sequence and protein sequences
        :param mut_df: data frame contains all the mutations
        :return: mut_df with syn label
        """
        # select subset of orfs with mut on this plate
        merge_mut = pd.merge(mut_df, all_df, how="left", on="orf_name")
        # for each pos, assign codon
        codon = [(int(i)//3)+1 if (int(i)%3 != 0) else int(i)/3 for i in merge_mut["pos"].tolist()]
        merge_mut["codon"] = codon
        # first group by ORF name
        # for each group, assign codon
        # SNP
        snp = merge_mut[merge_mut["label"] == "SNP"]
        # add a column for type
        snp["type"] = None
        grouped_snp = snp.groupby(["orf_name", "codon"]).apply(self._assign_syn)
        # get indel table
        indel = merge_mut[merge_mut["label"] == "indel"]
        indel["type"] = "indel"

        # join two table
        joined = pd.concat([grouped_snp, indel])
        return joined

    def _assign_syn(self, group):
        """
        Assign syn/non syn to each variant
        """
        table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }

        codon_seq = group["cds_seq"].values[0][
                        int(group["codon"].values[0] - 1) * 3:int(group["codon"].values[0]) * 3]
        if "N" in codon_seq:
            group["type"] = None
            return group
        pro = table[codon_seq]
        mut_codon = list(codon_seq)
        if group.shape[0] == 1:
            mut_pos = int(group["pos"]) % 3 - 1
            mut_codon[mut_pos] = group["alt"].values[0]
            mut_pro = table["".join(mut_codon)]
            if pro == mut_pro:
                group["type"] = "syn"
            else:
                    # # check if this maps the grch37 reference sequence
                    # grch37_refseq = group["grch37_filled"].values[0][:-3]
                    # codon_37 = grch37_refseq[int(group["codon"].values[0] - 1) * 3:int(group["codon"].values[0]) * 3]
                    # if mut_codon[mut_pos] == codon_37[mut_pos]:
                    #     track_syn.append("mapped_diffref")
                    # else:
                    #print(group.gene_ID.values[0])
                group["type"] = "non_syn"
        else: # two or three variants in the same codon
            group["mut_pos"] = group["pos"].astype(int) % 3 - 1
            for i, r in group.iterrows():
                mut_codon[r["mut_pos"]] = r["alt"]
            mut_pro = table["".join(mut_codon)]
            if pro == mut_pro:
                group["type"] = "syn"
            else:
                    # get grch37 codon
                    # check if this maps the grch37 reference sequence
                    # grch37_refseq = group["grch37_filled"].values[0][:-3]
                    # codon_37 = grch37_refseq[int(group["codon"].values[0] - 1) * 3:int(group["codon"].values[0]) * 3]
                    #
                    # if "".join(mut_codon) == codon_37:
                    #     track_syn += ["mapped_diffref"] * group["mut_pos"].shape[0]
                    # else:
                group["type"] = "non_syn"

        return group
