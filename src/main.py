from __future__ import division
from alignment import *
from variant_call import *
from analysis import *
from plot import *
from sup import *
import logging.config


def write_full_cover(plate_name, all_genes, full_cover_genes, snp, indel, ref_dict, output_file):
	"""
	
	:param plate_name: str
	:param all_genes: dictionary 
	:param full_cover_genes: 
	:param snp: 
	:param indel: 
	:param ref_dict: 
	:param output_file: 
	:return: 
	"""
	with open(output_file, "a") as output_file:
		output_file.write("plate_name,gene_name,gene_length,aligned_length,alignment_rate,total_read_count,average_read_depth,number_of_SNP,number_of_INDEL\n")
		for gene in all_genes.keys():
			if gene in full_cover_genes:
				all_genes = full_cover_genes
			if gene in snp.keys():
				line = [plate_name,
						gene,
						str(ref_dict[gene]),
						str(all_genes[gene][0]),
						str(all_genes[gene][0]/ref_dict[gene]),
						str(all_genes[gene][1]),
						str(all_genes[gene][2]),
						str(len(snp[gene])),
						"0"]
			else:
				line = [plate_name,
						gene,
						str(ref_dict[gene]),
						str(all_genes[gene][0]),
						str(all_genes[gene][0]/ref_dict[gene]),
						str(all_genes[gene][1]),
						str(all_genes[gene][2]),
						"0",
						"0"]
			if gene in indel.keys():
				line[-1] = str(indel[gene])
			output_file.write(",".join(line)+"\n")


def main():
	# step 1
	# check reference files
	# check fastq files
	# log these information
	# create a dir for each alignment
	if os.path.exists("log"):
		shutil.rmtree("log")
	os.makedirs("log")
	logging.config.fileConfig("./src/logging.conf")
	logger = logging.getLogger("main")

	# step 2
	# alignment
	# if user want the sequnece file to be aligned first
	# if only want to call variants with existing reference and bam file, please set the variable to False in conf.py
	if ALIGN:
		alignment_obj = Alignment(all_reference, subset_reference, fastq_path, ALIGNMENT_SETTING)
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
				logger.info("Analyzing %s ...", file)
				full_cover, total_gene_count, all_gene_dict, ref_dict = get_full_cover(file)
				logger.info("Got fully aligned genes in %s", file)
				snp, indel, read_depth = filter_vcf(file, all_gene_dict)
				logger.info("Filtered %s", file)

				if remove_syn:
					# print "Here"
					dna_seq = get_dna_ref(ref_fasta)
					# print dna_seq
					filterd_snp = remove_synonymous(snp, dna_seq)
					# write this information to file
					# full_covered.csv
					# plate gene snp indel
					out_file = output + "filtered_synonymous_summary.csv"
					write_full_cover(file, all_gene_dict, full_cover, filterd_snp, indel, ref_dict, out_file)

				else:
					out_file = output + "summary.csv"
					write_full_cover(file, all_gene_dict, full_cover, snp, indel, ref_dict, out_file)

				# take top 5 genes based on read depth
				# plot_top_n(snp, indel, read_depth, 3)
				# logger.info("Plot generated")
				logger.info("Analaysis done for %s \n", file)

if __name__ == "__main__":
	main()