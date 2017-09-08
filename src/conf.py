import matplotlib.pyplot as plt
import sys
import os
import vcf
import re
import plotly.plotly as py
import plotly.graph_objs as go
import seaborn as sns
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.patches as mpatches

# path to the folder that contains fastq files
fastq_path = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/fastq/"
# path to reference file
reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/ref/ORF_reference_pDONOR"
# path to output directory the last`/` is required!!!
output = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output/"

# pattern in file name that can be used to identify R1 and R2
file_name_pattern = "scORFeome-HIP-[0-9]+_S[0-9]+_L[0-9]+_"