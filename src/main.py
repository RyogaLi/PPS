from alignment import *
from variant_call import *

def main():
	# step 1
	# check reference files
	# check fastq files
	# log these information

	# step 2
	# alignment
	# if user want the sequnece file to be aligned first
	# if only want to call variants with existing reference and bam file, please set the variable to False in conf.py
	if ALIGN:
		alignment_obj = Alignment(all_reference, fastq_path, ALIGNMENT_SETTING)
		alignment_obj._main()

	# step 3
	# variant call
	if VARIANT_CALL:
		variant_caller = VariantCall(all_reference + ".fasta")
		variant_caller._main()

	# step 4
	# analysis




if __name__ == "__main__":
	main()