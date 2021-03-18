#

import pandas as pd
import os

class SnpAnalysis(object):

    def __init__(self, orf_id, clinvar):
        # file contains orf_id
        self._orf_id = orf_id
        # clinvar database
        self._clinvar = clinvar

    def _orf_names(self):
        orfs_df = pd.read_table(self._orf_id)
        self._orf_ID_name = orfs_df.set_index('orf_id')['ensemble_ID'].to_dict()

    def _load_clinvar(self):
        clin_raw = pd.read_table(self._clinvar)
        # print(clin_raw)
        # select only snp
        snp = clin_raw.loc[clin_raw['Type'] == "single nucleotide variant"]
        # useful cols: ["ID", "Name", "GeneSymbol", "ClinicalSignificance", "Start", "Stop", "ReferenceAllele", "AlternateAllele"]
        snp_simp = snp[["#AlleleID", "Name", "GeneSymbol", "ClinicalSignificance", "ReferenceAllele", "AlternateAllele"]]
        # print snp_simp
        # pathogenic / likely pathogenic
        snp_simp["Name_tmp"] = snp_simp["Name"].str.split(":").str.get(1)
        snp_simp["Name_nuc"], snp_simp["Name_aa"] = snp_simp["Name_tmp"].str.split(" ",1).str
        snp_simp = snp_simp.drop(["Name_tmp"], axis=1)
        snp_simp["Name_nuc"] = snp_simp["Name_nuc"].str.replace("c.","")
        snp_simp["Name_aa"] = snp_simp["Name_aa"].str.replace("(", "")
        snp_simp["Name_aa"] = snp_simp["Name_aa"].str.replace(")", "")

        rege = r'[0-9]+\+[0-9]+'
        rpl = lambda m:str(int(m.group(0).split("+")[0])+int(m.group(0).split("+")[1]))
        snp_simp["Name_nuc"] = snp_simp["Name_nuc"].str.replace(rege, rpl)
        # rpl = lambda m:m.group(0).split("+")[0]
        rege = r'[0-9]+\-[0-9]+'
        rpl = lambda m: str(int(m.group(0).split("-")[0]) - int(m.group(0).split("-")[1]))
        # print rpl
        snp_simp["Name_nuc"] = snp_simp["Name_nuc"].str.replace(rege, rpl)


        self._clinvar_patho = snp_simp.loc[snp_simp["ClinicalSignificance"].str.contains("Pathogenic|pathogenic") == True]
        self._clinvar_benign = snp_simp.loc[snp_simp["ClinicalSignificance"].str.contains("Benign|benign") == True]

        # print self._clinvar_patho

    def _load_vcf(self, vcf):
        # print(vcf)
        # read vcf files
        vcf = pd.read_csv(vcf, header=None, error_bad_lines=False)
        # print(vcf)
        # assign col names
        vcf_simp = vcf[[0, 1, 3, 4, 7]]
        colnames = ["orf_id", "position", "reference", "alter", "info"]
        vcf_simp.columns = colnames
        # for snps, only select those with >2 reads supports
        # DP=2+
        vcf_snp_only = vcf_simp.loc[(vcf_simp["alter"].str.contains("\*") == False) &  (vcf_simp["info"].str.contains("INDEL") == False) & (vcf_simp["info"].str.contains("DP=1|DP=2") == False)]
        # modify the vcf format so it has the same cols as clinvar
        # gene_name posREF>ALT
        colnames = ["gene_name", "orf_id", "positionREF>ALT", "info"]
        transfromed = []
        # transfromed.append(colnames)
        for index, row in vcf_snp_only.iterrows():
            entry = []
            # find gene name
            try:
                gene_name = self._orf_ID_name[row["orf_id"]]
            except Exception:
                continue
            try:
                alt = row["alter"].split(",")
                for a in alt:
                    pos_RA = str(row["position"])+row["reference"]+">"+a
                    entry = [gene_name, row["orf_id"], pos_RA, row["info"]]
                    # print entry
                    transfromed.append(entry)
            except Exception:
                pos_RA = str(row["position"]) + row["reference"] + ">" + row["alt"]
                entry = [gene_name, row["orf_id"], pos_RA, row["info"]]
                transfromed.append(entry)
                # print entry

        vcf_snp_only = pd.DataFrame(transfromed, columns=colnames)

        return vcf_snp_only

    def _find_patho(self, vcf_snp_only):
        # get all targeted gene_names
        # change col name
        self._clinvar_patho = self._clinvar_patho.rename(index=str, columns={"GeneSymbol": "gene_name", "Name_nuc": "positionREF>ALT"})
        patho_snp = pd.merge(self._clinvar_patho, vcf_snp_only, how='inner', on=['gene_name', 'positionREF>ALT'])
        return patho_snp

    def _find_benign(self, vcf_snp_only):
        self._clinvar_benign = self._clinvar_benign.rename(index=str, columns={"GeneSymbol": "gene_name","Name_nuc": "positionREF>ALT"})
        benign_snp = pd.merge(self._clinvar_benign, vcf_snp_only, how='inner', on=['gene_name', 'positionREF>ALT'])
        return benign_snp

    def _main(self, vcf):
        # VA = SnpAnalysis(orf_id, clinvar)
        self._orf_names()
        self._load_clinvar()
        # print("File:",vcf)
        vcf_snp = self._load_vcf(vcf)
        # print vcf_snp.shape
        patho_snp = self._find_patho(vcf_snp)
        # print patho_snp.shape
        benign_snp = self._find_benign(vcf_snp)
        # print benign_snp.shape
        all_snp_clinvar = patho_snp.append(benign_snp)
        #
        return all_snp_clinvar

class VCFParser(object):

    def __init__(self, vcf_file):
        # vcf file from pileups
        self._vcf_file = vcf_file

    def _load(self):
        vcf = pd.read_table(self._vcf_file, header=None)
        # assign col names
        vcf_simp = vcf[[0, 1, 3, 4, 7]]
        colnames = ["orf_id", "position", "reference", "alter", "info"]
        vcf_simp.columns = colnames



if __name__ == '__main__':

    orf_id = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/new_orf_id_convert.txt"
    clinvar = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/180213_clinvar_variant_summary.txt"
    vcf = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/test_filtered.vcf"
    VA = SnpAnalysis(orf_id, clinvar)
    VA._main(vcf)
