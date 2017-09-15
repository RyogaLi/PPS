from __future__ import division
import matplotlib.pyplot as plt
import sys
import os
import vcf
import re
import subprocess
import plotly.plotly as py
import plotly.graph_objs as go
import seaborn as sns
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.patches as mpatches

# path to the folder that contains fastq files
fastq_path = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/fastq/"
# path to reference files
all_reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/combined_ref/ORF_combined_ref"
subset_reference = "" # future
# path to output directory the last`/` is required!!!
output = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_combined_ref/"

# pattern in file name that can be used to identify R1 and R2
file_name_pattern = "scORFeome-HIP-[0-9]+_S[0-9]+_L[0-9]+_"


# variables

# if the fastq files are paired
PAIRED = False

# setting for alignment
ALIGNMENT_SETTING = "SENSITIVE"
# ALIGNMENT_SETTING = "DEFAULT"

# If the sequence need to be aligned
# change this variable to False to avoid alignment if the sequences are already aligned
ALIGN = False

# If variant call is already done, you can set this to False
# Then the program will only do analysis
VARIANT_CALL = False


