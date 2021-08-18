#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca
import ast
import sys
sys.path.append('..')
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted, venn3, venn3_circles, _venn3
pd.options.mode.chained_assignment = None  # default='warn'

class PlotObjYeast(object):

    def __init__(self, inputdir, summary_file, mutation_file):
        """
        Initialize plot object
        :param inputdir: input directory contains all the output files from the server
        """
        self._dir = inputdir
         # file contains all fully covered genes
        # self._all_summary = pd.read_csv(os.path.join(self._dir, "all_summary.csv"))
        # # file contains all found genes
        # self._all_mut = pd.read_csv(os.path.join(self._dir, "all_mutations.csv"))
        #
        self._all_summary = summary_file
        self._all_mut = mutation_file

    def make_venn(self, orfs):
        """

        :return:
        """
        # file contains all fully covered genes
        all_summary = os.path.join(self._dir, "all_summary_plateORF.csv")
        # file contains all found genes
        all_found_summary = os.path.join(self._dir, "all_found_summary_plateORF.csv")

        # compare genes in all the targeted space (ORFs) vs all fully aligned
        all_targeted_unique_db = orfs["ORF_NAME_NODASH"].dropna().unique()

        all_found = pd.read_csv(all_found_summary)
        all_fully_covered = pd.read_csv(all_summary)

        all_found["gene_name"] = all_found["gene_name"].replace("-", "")
        all_found_genes = all_found["gene_name"].dropna().unique()

        all_fully_covered["gene_name"] = all_fully_covered["gene_name"].replace("-", "")
        all_fully_covered_genes = all_fully_covered["gene_name"].dropna().unique()
        venn3([set(all_targeted_unique_db), set(all_found_genes), set(all_fully_covered_genes)], set_labels=("all ORFs",
                                                                                                   "all found",
                                                                                                   "fully covered"))
        venn3_circles([set(all_targeted_unique_db), set(all_found_genes), set(all_fully_covered_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.savefig(os.path.join(self._dir, "./all_venn3_allfound.png"))
        plt.close()

        all_mut_summary = os.path.join(self._dir, "all_mutations.csv")
        all_mut = pd.read_csv(all_mut_summary)
        all_mut = all_mut[(all_mut["type"] != "syn") & (all_mut["type"] != "NA")]
        all_mut["gene_name"] = all_mut["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
        all_mut_genes = all_mut["gene_name"].dropna().unique()
        venn3([set(all_targeted_unique_db), set(all_found_genes), set(all_mut_genes)], set_labels=("all ORFs",
                                                                                                   "fully covered",
                                                                                                   "fully covered; "
                                                                                                   "with non-syn variants"))
        venn3_circles([set(all_targeted_unique_db), set(all_found_genes), set(all_mut_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.savefig(os.path.join(self._dir, "./all_venn3_nonsyn.png"))
        plt.close()

        # HIP subset
        all_HIP_targeted = orfs[orfs["db"] == "HIP"]["ORF_NAME_NODASH"].dropna().unique()
        all_found_hip = all_found[all_found["db"] == "HIP"]
        all_found_hip["gene_name"] = all_found_hip["gene_name"].replace("-", "")
        all_found_genes_hip = all_found_hip["gene_name"].dropna().unique()
        all_full_hip = all_fully_covered[all_fully_covered["db"] == "HIP"]
        all_full_hip["gene_name"] = all_full_hip["gene_name"].replace("-", "")
        all_fully_covered_genes_hip = all_full_hip["gene_name"].dropna().unique()
        plt.figure(figsize=(8, 5))
        vd = venn3([set(all_HIP_targeted), set(all_found_genes_hip), set(all_fully_covered_genes_hip)],
                   set_labels=(f"all HIP ORFs: {len(all_HIP_targeted)}", f"all found: {len(all_found_genes_hip)}", f"fully covered: {len(all_fully_covered_genes_hip)}"))
        venn3_circles([set(all_HIP_targeted), set(all_found_genes_hip), set(all_fully_covered_genes_hip)], linestyle='dashed',
                      linewidth=1, color="black")
        vd.get_label_by_id("100").set_x(-0.6)
        vd.get_label_by_id("100").set_y(0.3)
        vd.get_label_by_id("111").set_x(0.2)
        vd.get_label_by_id("11").set_y(-0.1)
        plt.savefig(os.path.join(self._dir, "./all_found_full_hip.png"))
        plt.close()

        # SGD subset
        all_SGD_targeted = orfs[orfs["db"] == "SGD"]["ORF_NAME_NODASH"].dropna().unique()
        all_found_sgd = all_found[all_found["db"] == "SGD"]
        all_found_sgd["gene_name"] = all_found_sgd["gene_name"].replace("-", "")
        all_found_genes_sgd = all_found_sgd["gene_name"].dropna().unique()

        all_full_sgd = all_fully_covered[all_fully_covered["db"] == "SGD"]
        all_full_sgd["gene_name"] = all_full_sgd["gene_name"].replace("-", "")
        all_fully_covered_genes_sgd = all_full_sgd["gene_name"].dropna().unique()
        plt.figure(figsize=(9, 5))
        vd = venn3([set(all_SGD_targeted), set(all_found_genes_sgd), set(all_fully_covered_genes_sgd)],
              set_labels=(f"all supp-SGD ORFs: {len(all_SGD_targeted)}", f"all found: {len(all_found_genes_sgd)}",
                          f"fully covered: {len(all_fully_covered_genes_sgd)}"))
        venn3_circles([set(all_SGD_targeted), set(all_found_genes_sgd), set(all_fully_covered_genes_sgd)], linestyle='dashed',
                      linewidth=1, color="black")
        vd.get_label_by_id("100").set_x(-0.6)
        vd.get_label_by_id("100").set_y(0.3)
        vd.get_label_by_id("111").set_x(0.2)
        vd.get_label_by_id("11").set_y(-0.1)
        plt.savefig(os.path.join(self._dir, "./all_found_full_sgd.png"))
        plt.close()

        # PROTGEN subset
        all_PROT_targeted = orfs[orfs["db"] == "PROTGEN"]["ORF_NAME_NODASH"].dropna().unique()
        all_found_prot = all_found[all_found["db"] == "PROTGEN"]
        all_found_prot["gene_name"] = all_found_prot["gene_name"].replace("-", "")
        all_found_genes_prot = all_found_prot["gene_name"].dropna().unique()

        all_full_prot = all_fully_covered[all_fully_covered["db"] == "PROTGEN"]
        all_full_prot["gene_name"] = all_full_prot["gene_name"].replace("-", "")
        all_fully_covered_genes_prot = all_full_prot["gene_name"].dropna().unique()
        plt.figure(figsize=(9, 5))
        vd = venn3([set(all_PROT_targeted), set(all_found_genes_prot), set(all_fully_covered_genes_prot)],
              set_labels=(f"all supp-PROT ORFs: {len(all_PROT_targeted)}", f"all found: {len(all_found_genes_prot)}",
                          f"fully covered: {len(all_fully_covered_genes_prot)}"))
        venn3_circles([set(all_PROT_targeted), set(all_found_genes_prot), set(all_fully_covered_genes_prot)], linestyle='dashed',
                      linewidth=1, color="black")
        vd.get_label_by_id("100").set_x(-0.6)
        vd.get_label_by_id("100").set_y(0.3)
        vd.get_label_by_id("111").set_x(0.2)
        vd.get_label_by_id("11").set_y(-0.1)
        plt.savefig(os.path.join(self._dir, "./all_found_full_prot.png"))
        plt.close()


        all_mut_summary = os.path.join(self._dir, "all_mutations.csv")
        all_mut = pd.read_csv(all_mut_summary)
        all_mut = all_mut[(all_mut["type"] != "syn") & (all_mut["type"] != "NA")]
        all_mut["gene_name"] = all_mut["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
        all_mut_genes = all_mut["gene_name"].dropna().unique()
        plt.figure(figsize=(8, 5))
        venn3([set(all_targeted_unique_db), set(all_fully_covered_genes), set(all_mut_genes)],
              set_labels=(f"all ORFs:{len(all_targeted_unique_db)}", f"fully covered: {len(all_fully_covered_genes)}",
                          f"fully covered; \n with non-syn variants: {len(all_mut_genes)}"))
        venn3_circles([set(all_targeted_unique_db), set(all_fully_covered_genes), set(all_mut_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.savefig(os.path.join(self._dir, "./all_venn3_nonsyn.png"))
        plt.close()


        all_mut_hip = all_mut[all_mut["db"] == "HIP"]
        all_mut_hip["gene_name"] = all_mut_hip["gene_name"].replace("-", "")
        all_mut_hip_genes = all_mut_hip["gene_name"].dropna().unique()
        plt.figure(figsize=(8, 5))
        venn3([set(all_HIP_targeted), set(all_fully_covered_genes_hip), set(all_mut_hip_genes)],
              set_labels=(f"all HIP ORFs:{len(all_HIP_targeted)}", f"fully covered: {len(all_fully_covered_genes_hip)}",
                          f"fully covered; \n with non-syn variants: {len(all_mut_hip_genes)}"))
        venn3_circles([set(all_HIP_targeted), set(all_fully_covered_genes_hip), set(all_mut_hip_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.savefig(os.path.join(self._dir, "./HIPORFs_venn3_nonsyn.png"))
        plt.close()

        all_mut_sgd = all_mut[all_mut["db"] == "SGD"]
        all_mut_sgd["gene_name"] = all_mut_sgd["gene_name"].replace("-", "")
        all_mut_sgd_genes = all_mut_sgd["gene_name"].dropna().unique()
        plt.figure(figsize=(9, 5))
        venn3([set(all_SGD_targeted), set(all_fully_covered_genes_sgd), set(all_mut_sgd_genes)],
              set_labels=(f"all supp-SGD ORFs: {len(all_SGD_targeted)}", f"fully covered:"
                                                                    f" {len(all_fully_covered_genes_sgd)}",
                          f"fully covered; \n with non-syn variants: {len(all_mut_sgd_genes)}"))
        venn3_circles([set(all_SGD_targeted), set(all_fully_covered_genes_sgd), set(all_mut_sgd_genes)], linestyle='dashed',
                      linewidth=1, color="black")

        plt.savefig(os.path.join(self._dir, "./SGDORFs_venn3_nonsyn.png"))
        plt.close()

        all_mut_prot = all_mut[all_mut["db"] == "PROTGEN"]
        all_mut_prot["gene_name"] = all_mut_prot["gene_name"].replace("-", "")
        all_mut_prot_genes = all_mut_prot["gene_name"].dropna().unique()
        plt.figure(figsize=(9, 5))
        venn3([set(all_PROT_targeted), set(all_fully_covered_genes_prot), set(all_mut_prot_genes)],
              set_labels=(f"all supp-PROT ORFs: {len(all_PROT_targeted)}", f"fully covered:"
                                                                        f" {len(all_fully_covered_genes_prot)}",
                          f"fully covered; \n with non-syn variants: {len(all_mut_prot_genes)}")
              )
        venn3_circles([set(all_PROT_targeted), set(all_fully_covered_genes_prot), set(all_mut_prot_genes)], linestyle='dashed',
                      linewidth=1, color="black")

        plt.savefig(os.path.join(self._dir, "./PROTORFs_venn3_nonsyn.png"))
        plt.close()

        ################################################################################
        # all mut
        all_mut_summary = os.path.join(self._dir, "all_mutations.csv")
        all_mut = pd.read_csv(all_mut_summary)
        all_mut["gene_name"] = all_mut["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
        all_mut_genes = all_mut["gene_name"].dropna().unique()
        venn3([set(all_targeted_unique_db), set(all_found_genes), set(all_mut_genes)], set_labels=("all ORFs",
                                                                                                   "fully covered",
                                                                                                   "fully covered; "
                                                                                                     "with any "
                                                                                                   "variants"))
        venn3_circles([set(all_targeted_unique_db), set(all_found_genes), set(all_mut_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.savefig(os.path.join(self._dir, "./all_venn3.png"))
        plt.close()

        all_mut_hip = all_mut[all_mut["db"] == "HIP"]
        all_mut_hip["gene_name"] = all_mut_hip["gene_name"].replace("-", "")
        all_mut_hip_genes = all_mut_hip["gene_name"].dropna().unique()

        venn3([set(all_HIP_targeted), set(all_found_genes_hip), set(all_mut_hip_genes)], set_labels=("all HIP ORFs",
                                                                                                     "fully covered",
                                                                                                     "fully covered; "
                                                                                                     "with "
                                                                                                     "any variants"))
        venn3_circles([set(all_HIP_targeted), set(all_found_genes_hip), set(all_mut_hip_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.savefig(os.path.join(self._dir, "./HIPORFs_venn3.png"))
        plt.close()

        all_mut_sgd = all_mut[all_mut["db"] == "SGD"]
        all_mut_sgd["gene_name"] = all_mut_sgd["gene_name"].replace("-", "")
        all_mut_sgd_genes = all_mut_sgd["gene_name"].dropna().unique()

        venn3([set(all_SGD_targeted), set(all_found_genes_sgd), set(all_mut_sgd_genes)], set_labels=("all Supp-SGD "
                                                                                                     "ORFs",
                                                                                                     "fully covered",
                                                                                                     "fully covered; "
                                                                                                     "with any "
                                                                                                     "variants"))
        venn3_circles([set(all_SGD_targeted), set(all_found_genes_sgd), set(all_mut_sgd_genes)], linestyle='dashed',
                      linewidth=1, color="black")

        plt.savefig(os.path.join(self._dir, "./SGDORFs_venn3.png"))
        plt.close()

        all_mut_prot = all_mut[all_mut["db"] == "PROTGEN"]
        all_mut_prot["gene_name"] = all_mut_prot["gene_name"].replace("-", "")
        all_mut_prot_genes = all_mut_prot["gene_name"].dropna().unique()
        venn3([set(all_PROT_targeted), set(all_found_genes_prot), set(all_mut_prot_genes)], set_labels=("all Supp-PROT "
                                                                                                        "ORFs",
                                                                                                        "fully covered",
                                                                                                        "fully covered; "
                                                                                                     "with any "
                                                                                                        "variants"))
        venn3_circles([set(all_PROT_targeted), set(all_found_genes_prot), set(all_mut_prot_genes)], linestyle='dashed',
                      linewidth=1, color="black")

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
        g = sns.barplot(data=all_genes_stats, x="plate", y="% on plate fully covered", hue="aligned_to")
        plt.xticks(rotation=90, fontsize=15)
        plt.yticks(fontsize=15)
        plt.ylabel("all the fully covered unique ORFs", fontsize=15)
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
                label='% of ORFs fully coverd')
        plt.bar(r2, perc_variant, color='#FFC96F', width=barWidth, edgecolor='white', label='% fully covered; with '
                                                                                            'variants')
        plt.xticks(r1, labels=all_genes_stats["xlabel"].tolist(), rotation=30, fontsize=27)
        plt.yticks(fontsize=27)
        plt.ylabel("Fraction of targeted ORFs on each plate that are fully covered", fontsize=30)
        # Create legend & Show graphic
        plt.legend(fontsize=29)
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "full_with_variants_barplot.png"))
        plt.close()

    def plot_n_variants(self):
        # all_mut_summary = os.path.join(self._dir, "all_mutations.csv")
        all_mut = pd.read_csv(self._all_mut)
        all_mut = all_mut[(all_mut["type"] != "syn") & (all_mut["type"] != "NA")]
        all_mut["gene_name"] = all_mut["gene_ID"].str.extract(r"(.*)-[A-Z]+-[1-9]")
        v_counts = all_mut[["gene_name", "db"]].value_counts().to_frame().reset_index()
        print(v_counts)
        v_counts.columns = ["gene_name", "db", "v_count"]
        v_counts["category"] = "1"
        v_counts.loc[v_counts["v_count"] == 2, "category"] = "2"
        v_counts.loc[v_counts["v_count"] > 2, "category"] = "3+"

        v_plot = v_counts[["db", "category"]].value_counts().to_frame().reset_index()
        v_plot.columns = ["subset", "category", "n_gene"]
        v_plot["subset"] = v_plot["subset"].replace({"SGD": "Supp-SGD", "PROTGEN": "Supp-PROT"})
        print(v_plot)
        sns.set_theme(style="whitegrid", font_scale=1.5)
        plt.figure(figsize=(8, 6))
        g = sns.barplot(data=v_plot, x="category", y="n_gene", hue="subset", palette="Blues")
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.ylabel("Number of unique ORFs", fontsize=20)
        plt.xlabel("Number of non-syn variants", fontsize=20)
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "./n_gene_with_variants.png"))
        plt.close()

    def make_perc_coverred_dist(self):
        """
        from all summary, plot perc covered by sequencing data on bar plot
        For All, HIP, Supp-SGD and Supp-PROT
        :return:
        """
        print(self._all_summary.columns)
        print(self._all_mut.columns)
        not_fully_covered = self._all_summary.loc[(self._all_summary["aligned_perc"] != 1) & (self._all_summary["found"] == "y")]
        HIP = not_fully_covered[not_fully_covered["db"] == "HIP"]
        SGD = not_fully_covered[not_fully_covered["db"] == "SGD"]
        PROT = not_fully_covered[not_fully_covered["db"] == "PROTGEN"]
        print(self._all_summary.db.value_counts())

        fig, ax = plt.subplots(figsize=(10, 12))
        sns.displot(HIP.aligned_perc*100, bins=40, ax=ax, color="#084c61", edgecolor="#084c61")
        plt.title("HIP ORFs")
        plt.xlabel("Percent of ORF len aligned")
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "nfully_HIP_perc_dist.png"))
        plt.close()

        fig, ax = plt.subplots(figsize=(10, 12))
        sns.displot(SGD.aligned_perc*100, bins=40, ax=ax, color="#084c61", edgecolor="#084c61")
        plt.title("Supp-SGD ORFs")
        plt.xlabel("Percent of ORF len aligned")
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "nfully_SGD_perc_dist.png"))

        fig, ax = plt.subplots(figsize=(10, 12))
        sns.displot(PROT.aligned_perc*100, bins=40, ax=ax, color="#084c61", edgecolor="#084c61")
        plt.title("Supp-PROT ORFs")
        plt.xlabel("Percent of ORF len aligned")
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "nfully_PROT_perc_dist.png"))

        plt.close()


class PlotObjHuman(object):

    def __init__(self, inputdir):
        """
        Initialize plot object
        :param inputdir: input directory contains all the output files from the server
        """
        self._dir = inputdir
         # file contains all fully covered genes
        self._all_summary = pd.read_csv(os.path.join(self._dir, "all_summary.csv"))
        # file contains all found genes
        self._all_mut = pd.read_csv(os.path.join(self._dir, "all_mutations.csv"))

    def make_fully_covered_withmut_bar_plot(self):
        """

        :return:
        """
        genes_found_file = os.path.join(self._dir, "genes_stats.csv")
        all_genes_stats = pd.read_csv(genes_found_file)
        all_genes_stats["xlabel"] = all_genes_stats["plate"].str.replace("9-1-", "")
        perc_variant = all_genes_stats["n_fully_aligned_genes_with_any_mut"] / all_genes_stats["all_targeted_on_plate"]
        sns.set_theme(style="whitegrid", font_scale=1.5)
        plt.figure(figsize=(20, 14))
        # set width of bars
        barWidth = 0.45
        # Set position of bar on X axis
        r1 = np.arange(len(all_genes_stats["% in group fully aligned"]))
        r2 = [x + barWidth for x in r1]

        # Make the plot
        plt.bar(r1, all_genes_stats["% in group fully aligned"], color='#67A7E0', width=barWidth, edgecolor='white',
                label='% of ORFs fully covered')
        plt.bar(r2, perc_variant, color='#FFC96F', width=barWidth, edgecolor='white', label='% of ORFs fully covered; '
                                                                                            'with '
                                                                                            'any variants')
        plt.xticks(r1, labels=all_genes_stats["xlabel"].tolist(), rotation=30, fontsize=27)
        plt.yticks(fontsize=27)
        plt.ylabel("Fraction of targeted ORFs on each plate that are fully covered", fontsize=30)
        # Create legend & Show graphic
        plt.legend(fontsize=29)
        plt.ylim((0,1))
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "full_with_variants_barplot.png"))
        plt.close()

    def make_venn(self):
        """

        :return:
        """
        human_91 = "/Users/roujia/Documents/02_dev/04_HuRI/human ORF/20161117_ORFeome91_seqs.csv"
        human_ORFs = read_human_csv(human_91)
        all_summary_file = os.path.join(self._dir, "all_summary.csv")
        all_summary = pd.read_csv(all_summary_file)
        # all_found_summary = os.path.join(self._dir, "all_found_summary.csv")
        # compare genes in all the targeted space (ORFs) vs all fully aligned
        all_targeted_unique_db = human_ORFs[human_ORFs["entrez_gene_symbol"] != '-1']["entrez_gene_symbol"].dropna().unique()

        all_found_genes = all_summary[(all_summary["found"] == 'y') & (all_summary["entrez_gene_symbol"] != '-1')]["entrez_gene_symbol"].dropna().unique()
        all_fully_aligned_genes = all_summary[(all_summary["fully_covered"] == 'y') & (all_summary["entrez_gene_symbol"] != '-1')]["entrez_gene_symbol"].dropna(

        ).unique()
        # print(all_fully_aligned)
        # print(len(all_found_genes))
        # print(len(all_fully_aligned_genes))
        # print(len(all_targeted_unique_db))
        # print(set(all_found_genes) ^ set(all_targeted_unique_db))
        print(len(set.intersection(set(all_found_genes), set(all_targeted_unique_db))))
        plt.figure(figsize=(8, 5))
        vd = venn3([set(all_targeted_unique_db), set(all_found_genes), set(all_fully_aligned_genes)],
                   set_labels=(f"all ORFs: {len(all_targeted_unique_db)}", f"all found: {len(all_found_genes)}",
                               f"fully covered: {len(all_fully_aligned_genes)}"))
        venn3_circles([set(all_targeted_unique_db), set(all_found_genes), set(all_fully_aligned_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        vd.get_label_by_id("100").set_x(-0.6)
        vd.get_label_by_id("100").set_y(0.3)
        vd.get_label_by_id("11").set_y(-0.2)
        vd.get_label_by_id("111").set_x(0.2)
        vd.get_label_by_id("11").set_y(-0.1)
        plt.savefig(os.path.join(self._dir, "./allORFs_venn.png"))
        plt.close()

        # all_found_genes = all_fully_aligned["orf_name"].dropna().unique()
        venn2([set(all_targeted_unique_db), set(all_found_genes)], set_labels=("all ORFs", "all_fully_aligned"))
        plt.savefig(os.path.join(self._dir, "./allORFs_fullaligned_venn.png"))
        plt.close()

        all_mut_summary = os.path.join(self._dir, "all_mutations.csv")


        # all_mut = pd.read_csv(all_mut_summary)
        # all_mut = all_mut[(all_mut["type"] != "syn") & (all_mut["type"] != "NA") & (all_mut["type"] != "mapped_diffref")]
        # also filter out gnomad common variants
        # all_mut.exome = all_mut.exome.fillna("{'af': 0.000000001}").astype(str)
        # all_mut.genome = all_mut.genome.fillna("{'af': 0.000000001}").astype(str)
        # all_mut["exome_af"] = pd.DataFrame(all_mut.exome.values.tolist())
        # all_mut["exome_af"] = all_mut["exome"].str.extract(r'(\d+.\d+e?-?\d+)')
        # all_mut["genome_af"] =  all_mut["genome"].str.extract(r'(\d+.\d+e?-?\d+)')
        # print(all_mut[all_mut["exome"].notnull()][["exome", "exome_af"]])
        # # print(all_mut[all_mut["af"].notnull()]["af"])
        # all_mut["filled_af"] = all_mut["genome_af"].fillna(all_mut["exome_af"])
        #
        # all_mut["filled_af"] = all_mut["filled_af"].astype(float)
        # all_mut["filled_af"] = all_mut["filled_af"].fillna(0.0000001)
        # print(all_mut.shape)
        # # filter common variants
        # all_mut = all_mut[all_mut["filled_af"] == 0.0000001]
        # print(all_mut.shape)
        # all_mut_genes = all_mut["gene_ID"].dropna().unique()
        # plt.figure(figsize=(10, 5))
        # vd = venn3([set(all_targeted_unique_db), set(all_found_genes), set(all_mut_genes)],
        #            set_labels=(f"all ORFs: {len(all_targeted_unique_db)}", f"fully covered: {len(all_found_genes)}",
        #                        f"fully covered; \nwith non-syn mut; \nfiltered variants based on latest ENSEMBL "
        #                        f"refseq; \nremoved gnomAD common variants:\n"
        #                        f" {len(all_mut_genes)}"))
        # venn3_circles([set(all_targeted_unique_db), set(all_found_genes), set(all_mut_genes)], linestyle='dashed',
        #               linewidth=1, color="black")
        # plt.tight_layout()
        # plt.savefig(os.path.join(self._dir, "./allORFs_venn3.png"))
        # plt.close()

        all_mut = pd.read_csv(all_mut_summary)
        all_mut = all_mut[(all_mut["type"] != "syn") & (all_mut["type"] != "NA") & (all_mut["fully_covered"] == "y")]
        # also filter out gnomad common variants
        # all_mut.exome = all_mut.exome.fillna("{'af': 0.000000001}").astype(str)
        # all_mut.genome = all_mut.genome.fillna("{'af': 0.000000001}").astype(str)
        # all_mut["exome_af"] = pd.DataFrame(all_mut.exome.values.tolist())
        all_mut["exome_af"] = all_mut["exome"].str.extract(r'(\d+.\d+e?-?\d+)')
        # all_mut["exome_af"] = all_mut["exome"]
        all_mut["genome_af"] =  all_mut["genome"].str.extract(r'(\d+.\d+e?-?\d+)')
        # print(all_mut[all_mut["exome"].notnull()][["exome", "exome_af"]])
        # print(all_mut[all_mut["af"].notnull()]["af"])
        all_mut["filled_af"] = all_mut["genome_af"].fillna(all_mut["exome_af"])

        all_mut["filled_af"] = all_mut["filled_af"].astype(float)
        all_mut["filled_af"] = all_mut["filled_af"].fillna(0.0000001)
        # print(all_mut[all_mut["filled_af"].notnull()])
        # filter common variants
        all_mut = all_mut[all_mut["filled_af"] < 0.001]

        all_mut_genes = all_mut[all_mut["entrez_gene_symbol"] != '-1']["entrez_gene_symbol"].dropna().unique()
        plt.figure(figsize=(10, 5))

        vd = venn3([set(all_targeted_unique_db), set(all_fully_aligned_genes), set(all_mut_genes)],
                   set_labels=(f"all ORFs: {len(all_targeted_unique_db)}", f"fully covered: {len(all_fully_aligned_genes)}",
                               f"fully covered; \nwith non-syn mut; \nremoved gnomAD common variants:\n"
                               f" {len(all_mut_genes)}"))
        venn3_circles([set(all_targeted_unique_db), set(all_fully_aligned_genes), set(all_mut_genes)], linestyle='dashed',
                      linewidth=1, color="black")
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "./allORFs_venn3_onlygnomad.png"))
        plt.close()

        # genes with non_syn_ref
        all_mut_nonsyn_ref = all_mut[(all_mut["type"] == "non_syn_ref") & (all_mut["entrez_gene_symbol"] != '-1')][
            "entrez_gene_symbol"].dropna().unique()
        vd = venn3([set(all_targeted_unique_db), set(all_mut_nonsyn_ref), set(all_mut_genes)],
                   set_labels=(
                   f"all ORFs: {len(all_targeted_unique_db)}", f"{len(all_mut_nonsyn_ref)}",
                   f"fully covered; \nwith non-syn mut; \nremoved gnomAD common variants:\n"
                   f" {len(all_mut_genes)}"))
        venn3_circles([set(all_targeted_unique_db), set(all_mut_nonsyn_ref), set(all_mut_genes)],
                      linestyle='dashed',
                      linewidth=1, color="black")
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "./allORFs_venn3_variants.png"))
        plt.close()

    def make_perc_coverred_dist(self):
        """
        from all summary, plot perc covered by sequencing data on bar plot
        For All, HIP, Supp-SGD and Supp-PROT
        :return:
        """
        print(self._all_summary.columns)
        print(self._all_mut.columns)
        not_fully_covered = self._all_summary.loc[(self._all_summary["aligned_perc"] < 1) & (self._all_summary["found"] == "y")]
        print(self._all_summary["aligned_perc"])
        print(not_fully_covered.shape)
        fig, ax = plt.subplots(figsize=(10, 12))
        sns.displot(not_fully_covered.aligned_perc*100, bins=40, ax=ax, color="#084c61", edgecolor="#084c61")
        plt.title("Human 9.1 ORFs")
        plt.xlabel("Percent of ORF len aligned")
        plt.tight_layout()
        plt.savefig(os.path.join(self._dir, "nfully_human91_perc_dist.png"))


def read_human_csv(human91_ORFs):
    humanallORF = pd.read_csv(human91_ORFs)
    humanallORF = humanallORF[['orf_id', 'entrez_gene_id', 'Pool group #', 'entrez_gene_symbol', 'Mapped reads', 'Verified',
                                '# mut']]
    humanallORF = humanallORF.fillna(-1)
    humanallORF["entrez_gene_id"] = humanallORF["entrez_gene_id"].astype(int)
    humanallORF['orf_name'] = humanallORF['orf_id'].astype(str) + "_" + humanallORF['entrez_gene_id'].astype(str) + "_G0" + humanallORF['Pool group #'].astype(str) + "_" + humanallORF['entrez_gene_symbol'].astype(str)

    return humanallORF


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


def plot_main(inputdir, summary, mut):
    if args.m == "yeast":

        plot_obj = PlotObjYeast(inputdir, summary, mut)
        # HIP_target_ORFs = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/yeast_reference/HIP_targeted_ORFs.csv"
        # other_target_ORFs = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/yeast_reference/other_targeted_ORFs.csv"
        # orfs = read_yeast_csv(HIP_target_ORFs, other_target_ORFs)
        # plot_obj.make_fully_covered_withmut_bar_plot()
        # # # plot_obj.make_venn_variants(orfs)
        # plot_obj.make_venn(orfs)
        # plot_obj.plot_n_variants()
        plot_obj.make_perc_coverred_dist()
    else:
        plot_obj = PlotObjHuman(inputdir)
        # plot_obj.make_fully_covered_withmut_bar_plot()
        # plot_obj.make_venn()
        plot_obj.make_perc_coverred_dist()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plasmid pool sequencing analysis (plots)')
    parser.add_argument('-i', help='input dir')
    parser.add_argument('-s', help='summary')
    parser.add_argument('-mut', help='mut_summary')
    parser.add_argument('-m', help='mode')

    args = parser.parse_args()

    plot_main(args.i, args.s, args.mut)
