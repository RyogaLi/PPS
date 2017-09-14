from alignment import *
from variant_call import *
from analysis import *

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
	# after variant calling
	# filter vcf file and get fully aligned genes with their avg read depth

	# go thru all output directory
	dir_list = os.listdir(output)
	for dir in dir_list:
		if not os.path.isdir(output + "/" + dir): continue
		os.chdir(os.path.join(output,dir))
		for file in os.listdir("."):
			if file.endswith(".log") and os.stat(file).st_size == 0:
				os.remove(file)
			if file.endswith(".raw.vcf"):
				# full_cover is a dictionary contains key=gene ids for genes that are fully covered
				# value = [fully covered gene length, average read depth for this gene]
				full_cover, total_gene_count = get_full_cover(file)

				snp, indel, read_depth = filter_vcf(file, full_cover)
				# take top 5 genes based on read depth
				plot_top_n(snp, indel, read_depth, 3)


if __name__ == "__main__":
	main()