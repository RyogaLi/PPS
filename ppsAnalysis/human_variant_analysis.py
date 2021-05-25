#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')
import os
import operator
import re
import pandas as pd


class humanAnalysis(object):

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
                        #rd = re.search("DP=([0-9]+)", line[7])
                        #rd = rd.group(1)
                        gene_dict[line[0]] = 1
                    else:
                        # grep read depth information from INFO section
                        #rd = re.search("DP=([0-9]+)", line[7])
                        #rd = rd.group(1)
                        gene_dict[line[0]] += 1
                        #gene_dict[line[0]][1] += int(rd)
            remove_genes = gene_dict.copy()
            for key in gene_dict.keys():
                if gene_dict[key] < int(ref_dict[key]):
                    del remove_genes[key]
                #else:
                #    avg_rd = remove_genes[key][1] / remove_genes[key][0]
                #    remove_genes[key].append(avg_rd)
        return remove_genes, gene_dict, ref_dict

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
                    # don't use the raw DP
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
                    df["perc"] = df["read_count"].astype(float)/float(info_dict["DP"])
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
                    mut_count.append([l[0], l[1], l[3], mut_base,  l[5], mut_counts, info_dict["DP"], label])
                    filteredvcf.write(line)
        return mut_count

    def _process_mut(self, mut_df):
        """
        Based on the coding sequence and protein sequence, determine if the variant is a syn/non-syn variant,
        also call gnomAD API to see if the variant is common
        :param mut_df: data frame contains all the mutations
        :return: mut_df with syn label
        """
        # select subset of orfs with mut on this plate
        merge_mut = pd.merge(mut_df, self._orfs, how="left", left_on="gene_ID", right_on="ORF_id")
        # for each pos, assign codon
        codon = [(int(i) // 3) + 1 if (int(i) % 3 != 0) else int(i) / 3 for i in merge_mut["pos"].tolist()]
        merge_mut["codon"] = codon
        # first group by ORF name
        # for each group, assign codon
        # SNP
        snp = merge_mut[merge_mut["label"] == "SNP"]
        grouped = snp.groupby(["gene_ID", "codon"])
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
        # go through each codon, change the base
        track_syn = []
        for name, group in grouped:
            codon_seq = group["cds_seq"].values[0][
                        int(group["codon"].values[0] - 1) * 3:int(group["codon"].values[0]) * 3]
            if "N" in codon_seq:
                track_syn.append("NA")
                continue
            pro = table[codon_seq]
            mut_codon = list(codon_seq)
            if group.shape[0] == 1:
                mut_pos = int(group["pos"]) % 3 - 1
                mut_codon[mut_pos] = group["alt"].values[0]
                mut_pro = table["".join(mut_codon)]
                if pro == mut_pro:
                    track_syn.append("syn")
                else:
                    track_syn.append("non_syn")
            else:
                group["mut_pos"] = group["pos"].astype(int) % 3 - 1
                for i, r in group.iterrows():
                    mut_codon[r["mut_pos"]] = r["alt"]
                mut_pro = table["".join(mut_codon)]
                if pro == mut_pro:
                    track_syn += ["syn"] * group["mut_pos"].shape[0]
                else:
                    track_syn += ["non_syn"] * group["mut_pos"].shape[0]

        snp["type"] = track_syn
        # get indel table
        indel = merge_mut[merge_mut["label"] == "indel"]
        indel["type"] = "indel"

        # join two table
        joined = pd.concat([snp, indel])

    def _get_gnomAD(self):
        """
        Use gnomad API to get variants from gnomAD
        :return:
        """
        pass

    def _get_clinvar(self):
        """
        Process clinvar db and get variants from clinvar
        :return:
        """
        pass
