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

    def __init__(self, input_vcf, basename):
        """
        Take vcf file from input dir
        :param input_dir: Input dir contains vcf files, sam files and bam files
        """
        self._rawvcf = input_vcf
        self._basename = basename

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
                        rd = re.search("DP=([0-9]+)", line[7])
                        rd = rd.group(1)
                        gene_dict[line[0]] = [1, int(rd)]
                    else:
                        # grep read depth information from INFO section
                        rd = re.search("DP=([0-9]+)", line[7])
                        rd = rd.group(1)
                        gene_dict[line[0]][0] += 1
                        gene_dict[line[0]][1] += int(rd)
            remove_genes = gene_dict.copy()
            for key in gene_dict.keys():
                if gene_dict[key][0] < int(ref_dict[key]):
                    del remove_genes[key]
                else:
                    avg_rd = remove_genes[key][1] / remove_genes[key][0]
                    remove_genes[key].append(avg_rd)
        return remove_genes, gene_dict, ref_dict

    def _process_target_genes(self):
        pass
