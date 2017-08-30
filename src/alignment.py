from conf import *
import os
import sys
import logging.config
import shutil

class Alignment(object):

	def __init__(self, reference, fastq_path, setting="DEFAULT"):
		self._reference = reference
		self._fastq_path = fastq_path
		self._setting = setting

	def _align(self, fastq):
		self._basename = os.path.basename(fastq).split(".")[0]
		# create a dir for each alignment
		if os.path.exists(output+self._basename):
			shutil.rmtree(output+self._basename)
		os.makedirs(output+self._basename)
		os.chdir(output+self._basename)
		if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
			command = "bowtie2 -a " + " -x " + self._reference +" -U " + fastq  + " -S " + self._basename + ".sam " + "2> bowtie_DEFAULT_" + self._basename + ".log"
			os.system(command)
		elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
			command = "bowtie2 -a -p 20 -x " + self._reference +" --local --very-sensitive-local -U " + fastq + " -S " + self._basename + ".sam " + "2> bowtie_SENSITIVE_" + self._basename + ".log"
			os.system(command)
		else:
			command = "ERROR: please provide correct setting (DEFAULT/SENSITIVE)"
			sys.exit(command)

		# convert sam file to a sorted bam file out put from samtools are save in corresponding log files, sterr
		os.system("samtools view -bS "+ self._basename + ".sam > " + self._basename + ".bam 2> sam_to_bam.log")
		os.system("samtools sort " + self._basename + ".bam -o " + self._basename + "_sorted.bam 2> sort_bam.log")
		# creating a bam index file
		os.system("samtools index " + self._basename + "_sorted.bam " + self._basename + "_sorted.bai 2> generate_index.log")

		return command

	def _gene_count(self):
		gene_count = {}
		with open(self._basename+".sam", "r") as sam_file:
			for line in sam_file:
				if "@" in line: continue
				line = line.split("\t")
				if line[2] in gene_count.keys():
					gene_count[line[2]] += 1
				else:
					gene_count[line[2]] = 1
		with open("gene_count.txt", "w") as output:
			output.write("gene_name\tgene_count\n")
			for key in gene_count.keys():
				output.write(key+"\t"+str(gene_count[key])+"\n")

	def _remove_errlog(self):
		"""
		remove error log file if it's empty
		:return: 
		"""
		pass

	def _main(self):
		# create a dir for each alignment
		if os.path.exists("./log"):
			shutil.rmtree("./log")
		os.makedirs("./log")
		# init logging
		logging.config.fileConfig("./logging.conf")
		logger = logging.getLogger("alignment")
		# create log directory


		# fastq files list
		fastq_files = []
		for file in os.listdir(fastq_path):
			if file.endswith(".fastq"):
				fastq_files.append(os.path.join(fastq_path, file))

		# for each fastq file in the list
		# align to ref
		# put bowtie debug into separate file
		# alignment.log: if the alignment complete successfully
		for fastq in fastq_files:
			logger.info("started aligning "+os.path.basename(fastq))
			command = self._align(fastq)
			logger.info(command)
			logger.info("alignment finished for "+os.path.basename(fastq))
			self._gene_count()



if __name__ == "__main__":
	# get all the names of fastq file
	fastq_path = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_DK/"
	reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/ref/ORF_reference_pDONOR"
	alignment_obj = Alignment(reference, fastq_path)
	alignment_obj._main()
	# alignment_obj._gene_count()

