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

    def _get_ensembl_dna(self):
        """
        go through the input df, get cds sequence from ensembl 
        """
        server = "https://rest.ensembl.org"
        missing = []
        output_file = os.path.join(self._data_path, "ensembl_seq.csv")
        with open(output_file, "w") as output_f:
            output_f.write("enst_id,enst_version,cds_seq\n")
            for enst in self._input_df["ensembl_transcript_id"].tolist():

                enst_id = enst.split(".")[0]
                ext = f"/sequence/id/{enst_id}?type=cdna"
 
                r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
                if not r.ok:
                    print(enst)
                    missing.append(enst)
                    continue 
                decoded = r.json()
                output_f.write(f"{decoded['id']},{decoded['version']},{decoded['seq']}\n")
        print(missing)

if __name__ == "__main__":
    humanref_withseq = "/home/rothlab/rli/02_dev/06_pps_pipeline/target_orfs/20180524_DK_ReferenceORFeome_human_withensemblID.csv"
    data_path = "/home/rothlab/rli/02_dev/06_pps_pipeline/publicdb"

    get_cds = getCDS(humanref_withseq, data_path)
    get_cds._get_ensembl_dna()
