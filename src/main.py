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
	if ALIGN == True:
		alignment_obj = Alignment(reference, fastq_path, ALIGNMENT_SETTING)
		alignment_obj._main()

	# step 3
	# variant call and analysis
	variant_caller = VariantCall(reference + ".fasta")
	variant_caller._main()


if __name__ == "__main__":
	main()