from conf import *
from alignment import *
from variant_call import *
import os
import sys

def main():
	# step 1
	# check reference files
	# check fastq files
	# log these information

	# step 2
	# alignment
	alignment_obj = Alignment(reference, fastq_path)
	alignment_obj._main()

	# step 3
	# variant call
	variant_caller = VariantCall(reference + ".fasta")
	variant_caller._main()


if __name__ == "__main__":
	main()