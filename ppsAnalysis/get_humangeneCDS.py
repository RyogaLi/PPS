#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys

sys.path.append('..')

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
        input_dir = os.path.dirname(self._input)
        main_log = os.path.join(input_dir, "getCDS.log")
        log_obj = logthis.logit(log_f=main_log, log_level=loglevel)
        self._logger = log_obj.get_logger("main")

    def _get_exons(self):
        """
        Get exon information from ensembl
        """
        # read input df

        self._ccds2seq = os.path.join(self._data_path, "CCDS2Sequence.current.txt")

        if not os.path.exists(self._ccds2seq):
            download = "https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS2Sequence.current.txt"
            cmd = f"wget -O {self._ccds2seq} {download}"
            os.system(cmd)

        ccds_convert = pd.read_csv(self._ccds2seq, sep="\t")
        print(ccds_convert)
        print(self._input_df)
        print(ccds_convert["nucleotide_ID"])
        exit()
        ccds_id = ccds_convert[ccds_convert["nucleotide_ID"] == self._enst_can]

        user_enst = ""
        if ccds_id.empty:
            self._logger.info(f"cannot find ccds ID for {self._enst_can}")
            # try one version bak because CCDS db are not uptodate
            enst_V = self._enst_can.split(".")[1]
            self._enst_can = f"{self._enst_can.split('.')[0]}.{int(enst_V)-1}"
            self._logger.info(f"Trying an older version .. {self._enst_can}")
            ccds_id = ccds_convert[ccds_convert["nucleotide_ID"] == self._enst_can]

        while ccds_id.empty:
            self._logger.info(f"cannot find ccds ID for {self._enst_can}")
            user_enst = input("Please provide another ENST: (type exit to exit) ")
            if user_enst == "exit":
                exit(1)
            if user_enst == "next":
                return -1
            ccds_id = ccds_convert[ccds_convert["nucleotide_ID"] == user_enst]

            self._logger.info(f"Updated ENST: {user_enst}")
            self._enst_can = user_enst


        self._ccds_id = ccds_id["#ccds"].values.item()
        ccds_df = ccds_convert[ccds_convert["#ccds"] == self._ccds_id]
        # get NM id and NP id from CCDS current file
        self._nm = ccds_df.loc[(ccds_df["status_in_CCDS"] == "Accepted") & (ccds_df["source"] == "NCBI")]["nucleotide_ID"].values
        self._np = ccds_df.loc[(ccds_df["status_in_CCDS"] == "Accepted") & (ccds_df["source"] == "NCBI")]["protein_ID"].values

        # print(self._nm, self._np)
        if len(self._nm) != 1:
            self._logger.info("cannot find NM")
            self._logger.info(ccds_df[["#ccds", "nucleotide_ID", "protein_ID", "status_in_CCDS"]])
            self._logger.info(self._nm)
            user_enst = input("Please provide another NM: (type exit to exit) ")
            if user_enst == "exit":
                exit(1)
            self._nm = user_enst
        if len(self._np) != 1:
            self._logger.info("cannot find NP")
            self._logger.info(ccds_df[["#ccds", "nucleotide_ID", "protein_ID", "status_in_CCDS"]])
            self._logger.info(self._np)
            user_enst = input("Please provide another NP: (type exit to exit) ")
            if user_enst == "exit":
                exit(1)
            self._np = user_enst

        curCCDS = os.path.join(self._data_path, "CCDS.current.txt")
        if not os.path.isfile(curCCDS):
        #     # print(f"ccds file not found, downloading from cCDS..", file=log_file)
            download = "https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt"
            cmd = f"wget -O {curCCDS} {download}"
            os.system(cmd)
        self._curCCDS = curCCDS

        # read ccds into df
        ccds_df = pd.read_csv(self._curCCDS, sep="\t")
        get_cds = ccds_df[ccds_df["ccds_id"] == self._ccds_id]

        if get_cds.empty:
            raise ValueError(f"{self._ccds_id} not found in CCDS.current.")

        self._cds_loc = get_cds["cds_locations"].tolist()[0]
        self._cds_loc = self._cds_loc.strip("][").split(", ")
        self._nc = get_cds["nc_accession"].values.item()

        return 0

    def _get_dna(self, chrom, start, end):
        """
        get dna sequences from ucsc genome browser based on
        chrom, start, end
        """
        server="https://api.genome.ucsc.edu/"
        ext = f"getData/sequence?genome=hg38;chrom={chrom};start={start};end={int(end)+1}"
        r = requests.get(server+ext, headers={"Content-Type": "application/json"})
        #
        while not r.ok:
              #r.raise_for_status()
            time.sleep(200)
            r = requests.get(server+ext, headers={"Content-Type": "application/json"})
        decoded = r.json()
        dna_seq = Seq(decoded["dna"])
        #
        return dna_seq