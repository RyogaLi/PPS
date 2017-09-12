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
reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/old_ref/ORF_with_pDONR"
# path to output directory the last`/` is required!!!
output = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_old_ref/"

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

