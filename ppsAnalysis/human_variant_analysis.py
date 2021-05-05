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
                    print(df)

                    if len(l[3]) > 1:
                        label = "indel"
                    elif len(l[4].split(",")) == 1 and len(l[4]) > 1:
                        label = "indel"
                    else:
                        label = "SNP"
                    # track how many variants for each gene (with more than 10 reads mapped to it)
                    mut_count.append([l[0], l[1], l[3], l[4], l[5], label])
                    filteredvcf.write(line)
        return mut_count