#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys

sys.path.append('..')
import requests
import os
import pandas as pd
from ppsAnalysis import logthis

class getCDS(object):

    def __init__(self, input_file, data_path, loglevel="DEBUG"):
        """
        :param input_file: contains all targeted genes and ensembl ID
        :param data_path: path contains ccds data downloaded from ncbi db
        :param loglevel
        """
        self._input_df = pd.read_csv(input_file)
        self._data_path = data_path
        input_dir = os.path.dirname(input_file)
        main_log = os.path.join(input_dir, "getCDS.log")
        log_obj = logthis.logit(log_f=main_log, log_level=loglevel)
        self._logger = log_obj.get_logger("main")

    def _get_ensembl_dna_38(self):
        """
        go through the input df, get cds sequence from ensembl 
        """
        server = "https://rest.ensembl.org"
        missing = []
        output_file = os.path.join(self._data_path, "ensembl_seq_grch38.csv")
        if os.path.isfile(output_file):
            return output_file
        with open(output_file, "w") as output_f:
            output_f.write("enst_id,enst_version,cds_seq\n")
            for enst in self._input_df["ensembl_transcript_id"].tolist():
                enst_id = enst.split(".")[0]
                ext = f"/sequence/id/{enst_id}?type=cds"
                r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

                if not r.ok:
                    missing.append(enst)
                    continue 
                decoded = r.json()
                output_f.write(f"{decoded['id']},{decoded['version']},{decoded['seq']}\n")
        missing_38 = os.path.join(self._data_path, "grch38_missing.csv")
        missing_df = pd.DataFrame({"grch38_missing":missing})
        missing_df.to_csv(missing_38)
        return output_file

    def _get_ensembl_dna_37(self):
        """
        go through the input df, get cds sequence from ensembl
        """
        server = "https://grch37.rest.ensembl.org"
        missing = []
        output_file = os.path.join(self._data_path, "ensembl_seq_grch37.csv")
        if os.path.isfile(output_file):
            return output_file
        with open(output_file, "w") as output_f:
            output_f.write("enst_id,enst_version,cds_seq\n")
            for enst in self._input_df["ensembl_transcript_id"].tolist():

                enst_id = enst.split(".")[0]
                ext = f"/sequence/id/{enst_id}?type=cds"

                r = requests.get(server + ext, headers={"Content-Type": "application/json"})

                if not r.ok:
                    print(enst)
                    missing.append(enst)
                    continue
                decoded = r.json()
                output_f.write(f"{decoded['id']},{decoded['version']},{decoded['seq']}\n")
        missing_37 = os.path.join(self._data_path, "grch37_missing.csv")
        missing_df = pd.DataFrame({"grch37_missing":missing})
        missing_df.to_csv(missing_37)
        return output_file

    def merge_seq(self, grch37_file, grch38_file):
        """
        Merge sequence found in grch37 and grch38 into our cds seq file
        :param grch37_file: file contains grch37 genes
        :param grch38_file: file contains grch38 genes
        :return: None
        """
        self._input_df[["enst", "version"]] = self._input_df["ensembl_transcript_id"].str.split(".",expand=True)

        # load grch37
        grch37_df = pd.read_csv(grch37_file)
        grch37_df.columns = ["enst37", "enst37_version", "cds_seq37"]
        grch37_df["ensembl_transcript_grch37"] = grch37_df['enst37'].astype(str) + "." + grch37_df['enst37_version'].astype(str)

        # load grch38
        grch38_df = pd.read_csv(grch38_file)
        grch38_df.columns = ["enst38", "enst38_version", "cds_seq38"]
        grch38_df["ensembl_transcript_grcht38"] = grch38_df['enst38'].astype(str) + "." + grch38_df['enst38_version'].astype(str)

        # merge cds to input file
        merge = pd.merge(self._input_df, grch37_df, how="left", left_on = "enst", right_on= "enst37")
        merge = pd.merge(merge, grch38_df, how="left", left_on="enst", right_on="enst38")

        # safe merge to file
        output_file = os.path.join(self._data_path, "merged_ensembl_sequence.csv")
        merge.to_csv(output_file, index=False)

if __name__ == "__main__":
    humanref_withseq = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/20180524_DK_ReferenceORFeome_human_withensemblID.csv"
    data_path = "/home/rothlab/rli/02_dev/06_pps_pipeline/publicdb"

    get_cds = getCDS(humanref_withseq, data_path)
    grch38_seqs = get_cds._get_ensembl_dna_38()
    grch37_seqs = get_cds._get_ensembl_dna_37()

    get_cds.merge_seq(grch37_seqs, grch38_seqs)
