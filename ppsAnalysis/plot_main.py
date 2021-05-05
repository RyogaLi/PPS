#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted, venn3, venn3_circles

class PlotObj(object):

    def __init__(self, inputdir):
        """
        Initialize plot object
        :param inputdir: input directory contains all the output files from the server
        """
        self._dir = inputdir

    def make_venn(self, orfs):
        """

        :return:
        """
        all_summary = os.path.join(self._dir, "all_summary_subsetORF.csv")
        # compare genes in all the targeted space (ORFs) vs all fully aligned
        all_targeted_unique_db = orfs["ORF_NAME_NODASH"].dropna().unique()
        all_found = pd.read_csv(all_summary)
        all_found["gene_name"] = all_found["gene_name"].replace("-", "")
        all_found_genes = all_found["gene_name"].dropna().unique()
        venn2([set(all_targeted_unique_db), set(all_found_genes)], set_labels=("all ORFs", "all_fully_aligned"))
        plt.savefig(os.path.join(self._dir, "./allORFs_venn.png"))
        plt.close()

        # HIP subset
        all_HIP_targeted = orfs[orfs["db"] == "HIP"]["ORF_NAME_NODASH"].dropna().unique()
        all_found_hip = all_found[all_found["db"] == "HIP"]
        all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
        all_found_genes_hip = all_found_hip["gene_name"].dropna().unique()
        venn2([set(all_HIP_targeted), set(all_found_genes_hip)], set_labels=("all HIP ORFs", "all_fully_aligned"))
        plt.savefig(os.path.join(self._dir, "./HIPORFs_venn.png"))
        plt.close()

        # SGD subset
        all_SGD_targeted = orfs[orfs["db"] == "SGD"]["ORF_NAME_NODASH"].dropna().unique()
        all_found_hip = all_found[all_found["db"] == "SGD"]
        all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
        all_found_genes_sgd = all_found_hip["gene_name"].dropna().unique()
        venn2([set(all_SGD_targeted), set(all_found_genes_sgd)], set_labels=("all SGD ORFs", "all_fully_aligned"))
        plt.savefig(os.path.join(self._dir, "./SGDORFs_venn.png"))
        plt.close()

        # PROTGEN subset
        all_PROT_targeted = orfs[orfs["db"] == "PROTGEN"]["ORF_NAME_NODASH"].dropna().unique()
        all_found_hip = all_found[all_found["db"] == "PROTGEN"]
        all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
        all_found_genes_prot = all_found_hip["gene_name"].dropna().unique()
        venn2([set(all_PROT_targeted), set(all_found_genes_prot)], set_labels=("all PROTGEN ORFs", "all_fully_aligned"))
        plt.savefig(os.path.join(self._dir, "./PROTORFs_venn.png"))
        plt.close()

        all_mut_summary = os.path.join(self._dir, "all_mutations.csv")
        all_mut = pd.read_csv(all_mut_summary)
        all_mut["gene_name"] = all_mut["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
        all_mut_hip = all_mut[all_mut["db"] == "HIP"]
        all_mut_hip["gene_name"] = all_mut_hip["gene_name"].replace("-", "")
        all_mut_hip_genes = all_mut_hip["gene_name"].dropna().unique()

        venn3([set(all_HIP_targeted), set(all_found_genes_hip), set(all_mut_hip_genes)], set_labels=("all HIP ORFs",
                                                                                                 "fully aligned",
                                                                                                 "fully aligned with "
                                                                                                 "mut"))
        venn3_circles([set(all_HIP_targeted), set(all_found_genes_hip), set(all_mut_hip_genes)], linestyle='dashed', linewidth=1, color="grey")
        plt.savefig(os.path.join(self._dir, "./HIPORFs_venn3.png"))
        plt.close()


        all_mut_sgd = all_mut[all_mut["db"] == "SGD"]
        all_mut_sgd["gene_name"] = all_mut_sgd["gene_name"].replace("-", "")
        all_mut_sgd_genes = all_mut_sgd["gene_name"].dropna().unique()

        venn3([set(all_SGD_targeted), set(all_found_genes_sgd), set(all_mut_sgd_genes)], set_labels=("all SGD ORFs",
                                                                                                 "fully aligned",
                                                                                                 "fully aligned with "
                                                                                                 "mut"))
        plt.savefig(os.path.join(self._dir, "./SGDORFs_venn3.png"))
        plt.close()

        all_mut_prot = all_mut[all_mut["db"] == "PROTGEN"]
        all_mut_prot["gene_name"] = all_mut_prot["gene_name"].replace("-", "")
        all_mut_prot_genes = all_mut_prot["gene_name"].dropna().unique()
        venn3([set(all_PROT_targeted), set(all_found_genes_prot), set(all_mut_prot_genes)], set_labels=("all PROT "
                                                                                                        "ORFs",
                                                                                                 "fully aligned",
                                                                                                 "fully aligned with "
                                                                                                 "mut"))
        plt.savefig(os.path.join(self._dir, "./PROTORFs_venn3.png"))
        plt.close()

    def make_fully_covered_bar_plot(self):
        """
        Make bar plot for all the fully aligned ORFs
        :return:
        """
        genes_found_file = os.path.join(self._dir, "genes_stats.csv")
        all_genes_stats = pd.read_csv(genes_found_file)
        sns.set_theme(style="whitegrid", font_scale=1.5)
        plt.figure(figsize=(20, 14))
        g = sns.barplot(data=all_genes_stats, x="plate", y="% on plate fully aligned", hue="aligned_to")
        plt.xticks(rotation=90, fontsize=15)
        plt.yticks(fontsize=15)
        plt.ylabel("all the fully aligned unique ORFs", fontsize=15)
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "./perc_full_matched.png"))
        plt.close()

    def make_fully_covered_withmut_bar_plot(self):
        """

        :return:
        """
        genes_found_file = os.path.join(self._dir, "genes_stats.csv")
        all_genes_stats = pd.read_csv(genes_found_file)
        all_genes_stats["xlabel"] = all_genes_stats["plate"].str.replace("scORFeome-", "")
        perc_variant = all_genes_stats["n_fully_aligned_genes_with_any_mut"] / all_genes_stats["all_targeted_on_plate"]
        sns.set_theme(style="whitegrid", font_scale=1.5)
        plt.figure(figsize=(20, 14))
        # set width of bars
        barWidth = 0.45
        # Set position of bar on X axis
        r1 = np.arange(len(all_genes_stats["% on plate fully aligned"]))
        r2 = [x + barWidth for x in r1]

        # Make the plot
        plt.bar(r1, all_genes_stats["% on plate fully aligned"], color='#67A7E0', width=barWidth, edgecolor='white',
                label='% fully aligned')
        plt.bar(r2, perc_variant, color='#FFC96F', width=barWidth, edgecolor='white', label='% fully aligned with '
                                                                                            'variants')
        plt.xticks(r1, labels=all_genes_stats["xlabel"].tolist(), rotation=30, fontsize=27)
        plt.yticks(fontsize=27)
        plt.ylabel("Fraction of targeted ORFs on each plate", fontsize=30)
        # Create legend & Show graphic
        plt.legend(fontsize=29)
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "full_with_variants_barplot.png"))
        plt.close()


def read_yeast_csv(HIP_target_ORFs, other_target_ORFs):
    """
    Join HIP data and other data into one df, remove unwanted columns
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
    return combined

def plot_main(inputdir):
    plot_obj = PlotObj(inputdir)
    HIP_target_ORFs = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/yeast_reference/HIP_targeted_ORFs.csv"
    other_target_ORFs = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/yeast_reference/other_targeted_ORFs.csv"
    orfs = read_yeast_csv(HIP_target_ORFs, other_target_ORFs)
    plot_obj.make_fully_covered_withmut_bar_plot()
    # plot_obj.make_venn_variants(orfs)
    plot_obj.make_venn(orfs)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis (plots)')
    parser.add_argument('-i', help='input dir')

    args = parser.parse_args()

    plot_main(args.i)
