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
                    else:
                        # grep read depth information from INFO section
                        # rd = re.search("DP=([0-9]+)", line[7])
                        # rd = rd.group(1)
                        gene_dict[line[0]] += 1
                        # gene_dict[line[0]][1] += int(rd)
            removed_genes = gene_dict.copy()
            for key in gene_dict.keys():
                if gene_dict[key] < int(ref_dict[key]):
                    del removed_genes[key]
                # else:
                #     avg_rd = remove_genes[key][1] / remove_genes[key][0]
                #     remove_genes[key].append(avg_rd)
        return removed_genes, gene_dict, ref_dict

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
                    if len(l[3]) > 1:
                        label = "indel"
                    elif len(mut_base) > 1:
                        label = "indel"
                    else:
                        label = "SNP"
                    # track how many variants for each gene (with more than 10 reads mapped to it)
                    mut_count.append([l[0], l[1], l[3], mut_base, l[5], label])
                    filteredvcf.write(line)
        return mut_count
    
    def _process_target_genes(self):
        pass

    def filter_full_vcf(self):
        """
        Filter vcf file (full, including ambiguous and non-variant sites)
        :return: mut_count_dict
        """
        filtered_vcf = self._rawvcf.replace("_raw", "_raw_filtered")
        mut_count_dict = {}
        with open(self._rawvcf, "r") as rawvcf:
            with open(filtered_vcf, "w") as filteredvcf:
                for line in rawvcf:
                    if "#" in line: continue
                    l = line.split()
                    # if alt is ref
                    if l[4] == "<*>":
                        continue
                    # get variant call with quality > 20
                    try:
                        qual = float(l[5])
                    except:
                        continue
                    if qual < 20: continue
                    # get DP for this position
                    # match_str = "DP=([0-9]+)"
                    # m = re.search(match_str, l[7])
                    # if m:
                    #     dp = int(m.group(1))
                    # if dp < 10: continue
                    info_title = l[-2].split(":")
                    info = l[-1].split(":")
                    info_dict = dict(zip(info_title, info))
                    print(info_dict)
                    if int(info_dict["DP"]) < 10: # informative read depth
                        continue
                    print("passed dp")
                    alt_bases = l[4].split(",")
                    AD = info_dict["AD"].split(",")
                    alt_depth = dict(zip(alt_bases, AD))
                    pass_alt_filter = False
                    print(all_depth)
                    for bp in alt_depth.keys():
                        freq = float(alt_depth[bp]) / float(info_dict["DP"])
                        if bp == "<*>" and freq > 0.5:
                            print(all_depth)
                            # more than half the reads support wt
                            continue
                        pass_alt_filter = True
                    if pass_alt_filter:
                        filtered_vcf.write(line)
                        # track how many variants for each gene (with more than 10 reads mapped to it)
                        if mut_count_dict.get(l[0], -1) == -1:
                            mut_count_dict[l[0]] = 1
                        else:
                            mut_count_dict[l[0]] += 1
        return  mut_count_dict
