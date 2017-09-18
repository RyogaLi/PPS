from alignment import *
from variant_call import *
from analysis import *
from plot import *
import logging.config


def write_full_cover(plate_name, all_genes, full_cover_genes, snp, indel, ref_dict, output_file):
	with open(output_file, "a") as output_file:
		for gene in all_genes.keys():
			if gene in full_cover_genes:
				all_genes = full_cover_genes
			if gene in snp.keys():
				line = [plate_name, gene, str(ref_dict[gene]), str(all_genes[gene][0]),str(all_genes[gene][1]),str(snp[gene]),"N/A"]
			else:
				line = [plate_name, gene, str(ref_dict[gene]),str(all_genes[gene][0]),str(all_genes[gene][1]), "N/A", "N/A"]
			if gene in indel.keys():
				line[-1] = str(indel[gene])
			output_file.write(",".join(line)+"\n")


def main():
	# step 1
	# check reference files
	# check fastq files
	# log these information
	# create a dir for each alignment
	if os.path.exists("./log"):
		shutil.rmtree("./log")
	os.makedirs("./log")
	logging.config.fileConfig("./src/logging.conf")
	logger = logging.getLogger("main")

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
			logger.info("Removed empty log files in %s", dir)
			if file.endswith(".raw.vcf"):
				# full_cover is a dictionary contains key=gene ids for genes that are fully covered
				# value = [fully covered gene length, average read depth for this gene]
				logger.info("Analyzing %s ...", file)
				full_cover, total_gene_count, all_gene_dict, ref_dict = get_full_cover(file)
				logger.info("Got fully aligned genes in %s", file)
				snp, indel, read_depth = filter_vcf(file, all_gene_dict)
				logger.info("Filtered %s", file)

				# write this information to file
				# full_covered.csv
				# plate gene snp indel
				out_file = open("/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/full_covered_gene.csv", "w")
				write_full_cover(file,all_gene_dict, full_cover, snp, indel, ref_dict, out_file)
				out_file.close()

				# take top 5 genes based on read depth
				# plot_top_n(snp, indel, read_depth, 3)
				# logger.info("Plot generated")
				logger.info("Analaysis done for %s \n", file)


if __name__ == "__main__":
	main()