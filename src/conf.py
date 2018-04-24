from __future__ import division
import matplotlib.pyplot as plt
import sys
import os
import vcf
import re
import subprocess
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
import seaborn as sns
import numpy as np
import argparse
import matplotlib.patches as mpatches
from time import gmtime, strftime
from scipy.stats.stats import pearsonr
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

######### FILE PATHs ##########

# path to the folder that contains fastq files
fastq =  "/Users/roujia/Desktop/20171219_PPS/"
r1 = None
r2 = None

# path to reference files
all_reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/fasta/missing_orf/missing_orf"

ref_fasta = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/fasta/missing_orf/missing_orf.fasta"

# path to output directory the last`/` is required!!!
output = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/P3_PPS_output/"

# pattern in file name that can be used to identify R1 and R2
# file_name_pattern = "scORFeome-HIP-[0-9]+_S[0 -9]+_L[0-9]+_"

######### VARIABLES ##########

# # if the fastq files are paired
# PAIRED = False
#
# # if you want to merge r1 and r2
# MERGE = True

# setting for alignment
ALIGNMENT_SETTING = "SENSITIVE"
# ALIGNMENT_SETTING = "DEFAULT"

# If the sequence need to be aligned
# change this variable to False to avoid alignment if the sequences are already aligned
ALIGN = True

# If variant call is already done, you can set this to False
# Then the program will only do analysis
VARIANT_CALL = True

# If you want to remove all the synomuous SNPs
# set this to True
remove_syn = True


# if you want to compare the variants with variants on Clinvar database
# set this to True
clinvar = True
orf_id_convert = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/new_orf_id_convert.txt"
clinvar_db = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/180213_clinvar_variant_summary.txt"


