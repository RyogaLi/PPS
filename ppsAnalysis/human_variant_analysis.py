#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys
sys.path.append('..')
import os
import operator
import re
import pandas as pd
import requests
import time

class humanAnalysis(object):

    def __init__(self, input_vcf, basename, orfs_df, ref):
        """
        Take vcf file from input dir
        :param input_: Input dir contains vcf files, sam files and bam files
        """
        self._rawvcf = input_vcf
        self._vvcf = input_vcf.replace("_raw", "_variants")
        self._basename = basename
        # contains all the targeted orfs for this sample
        self._orfs = orfs_df
        if ref == "grch37":
            self._seq_col = "grch37_filled"
        elif ref == "grch38":
            self._seq_col = "grch38_filled"
        else:
            self._seq_col = "cds_seq"
            
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
                    mut_count.append([l[0], l[1], l[3], mut_base, l[5], mut_counts, info_dict["DP"], label])
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
        merge_mut = pd.merge(mut_df, self._orfs, how="left", left_on="gene_ID", right_on="orf_name")
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
            codon_seq = group[self._seq_col].values[0][
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
                    # # check if this maps the grch37 reference sequence
                    # grch37_refseq = group["grch37_filled"].values[0][:-3]
                    # codon_37 = grch37_refseq[int(group["codon"].values[0] - 1) * 3:int(group["codon"].values[0]) * 3]
                    # if mut_codon[mut_pos] == codon_37[mut_pos]:
                    #     track_syn.append("mapped_diffref")
                    # else:
                    print(codon_seq, mut_codon)
                    print(pro, mut_pro)

                    track_syn.append("non_syn")
            else: # two or three variants in the same codon 
                group["mut_pos"] = group["pos"].astype(int) % 3 - 1
                for i, r in group.iterrows():
                    mut_codon[r["mut_pos"]] = r["alt"]
                mut_pro = table["".join(mut_codon)]
                if pro == mut_pro:
                    track_syn += ["syn"] * group["mut_pos"].shape[0]
                else:
                    # get grch37 codon 
                    # check if this maps the grch37 reference sequence
                    # grch37_refseq = group["grch37_filled"].values[0][:-3]
                    # codon_37 = grch37_refseq[int(group["codon"].values[0] - 1) * 3:int(group["codon"].values[0]) * 3]
                    #
                    # if "".join(mut_codon) == codon_37:
                    #     track_syn += ["mapped_diffref"] * group["mut_pos"].shape[0]
                    # else:
                    track_syn += ["non_syn"] * group["mut_pos"].shape[0]

        snp["type"] = track_syn
        # get indel table
        indel = merge_mut[merge_mut["label"] == "indel"]
        indel["type"] = "indel"

        # join two table
        joined = pd.concat([snp, indel])

        # get gnomad variants for genes in the table
        merge_gnomad = []
        gene_list = joined.gene_ID.unique().tolist()
        for gene in gene_list:
            gnomAD_variants = self._get_gnomAD(gene)
            pps_variants = joined[joined["gene_ID"] == gene]
            merge_df = pd.merge(pps_variants, gnomAD_variants[["ref", "alt", "cds_pos", "exome", "genome"]], how="left", left_on=["ref", "alt", "pos"], right_on=["ref", "alt", "cds_pos"], suffixes=["_pps", "_gnomad"])
            merge_gnomad.append(merge_df)
            # label variants with matching gnomAD ref and if they are common
        joined = pd.concat(merge_gnomad)
        
        return joined


    def _get_gnomAD(self, gene_ID):
        """
        Use gnomad API to get variants from gnomAD
        :return:
        """
        gene_name = gene_ID.split("_")[-1]
        print(gene_name)
        # use transcript id instead
        q = """
        {
            gene(reference_genome: GRCh37, gene_symbol: "%s"){
            variants(dataset: gnomad_r2_1) {
            consequence
            pos
            variantId
            hgvs
            ref
            alt
            hgvsc
            hgvsp
            genome {
            af
            }
            exome {
            af
            }
            }
            }

        }""" % gene_name

        # send request
        r = requests.post("https://gnomad.broadinstitute.org/api", json={'query': q})
        while r.status_code != 200:
            print(r.status_code)
            time.sleep(60)
            r = requests.post("https://gnomad.broadinstitute.org/api", json={'query': q})

        variants = r.json()
        if variants["data"]["gene"] is None:
            return pd.DataFrame(columns=['consequence', 'pos', 'variantId', 'hgvs', 'ref', 'alt', 'hgvsc', 'hgvsp', 'genome', 'exome', 'cds_pos'])

        variants_dict = variants["data"]["gene"]["variants"]
        # convert response to dataframe
        df = pd.DataFrame.from_dict(variants_dict)
        if df.empty:
            return pd.DataFrame(columns=['consequence', 'pos', 'variantId', 'hgvs', 'ref', 'alt', 'hgvsc', 'hgvsp', 'genome', 'exome', 'cds_pos'])

        coding_variants = df[df.hgvsp.notnull()]
        # extract cds position using regex
        coding_variants["cds_pos"] = coding_variants['hgvsc'].str.extract('(\d+)', expand=True)
        return coding_variants


    def _get_clinvar(self):
        """
        Process clinvar db and get variants from clinvar
        :return:
        """
        pass
