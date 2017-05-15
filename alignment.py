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

	def error_check(self):
		# check if fasta file exits
		# check if fastq file exits
		pass

	def _align(self, fastq):
		basename = os.path.basename(fastq).split(".")[0]
		# create a dir for each alignment
		if os.path.exists(output+basename):
			shutil.rmtree(output+basename)
		os.makedirs(output+basename)
		os.chdir(output+basename)
		if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
			command = "bowtie2 -a " + " -x " + self._reference +" -U " + fastq  + " -S " + basename + ".sam " + "2> bowtie_DEFAULT_" + basename + ".log"
			os.system(command)
		elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
			command = "bowtie2 -a -p 20 -x " + self._reference +" --local --very-sensitive-local -U " + fastq + "-S " + basename + ".sam " + "2> bowtie_SENSITIVE_" + basename + ".log"
			os.system(command)
		else:
			command = "ERROR: please provide correct setting (DEFAULT/SENSITIVE)"
			sys.exit(command)

		# convert sam file to a sorted bam file out put from samtools are save in corresponding log files
		os.system("samtools view -bS "+ basename + ".sam > " + basename + ".bam 2> sam_to_bam.log")
		os.system("samtools sort " + basename + ".bam -o " + basename + "_sorted.bam 2> sort_bam.log")
		# creating a bam index file
		os.system("samtools index " + basename + "_sorted.bam " + basename + "_sorted.bai 2> generate_index.log")

		return command

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
			break # test


if __name__ == "__main__":
	# get all the names of fastq file
	fastq_path = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_DK/"
	reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/ref/ORF_reference_pDONOR"
	alignment_obj = Alignment(reference, fastq_path)
	alignment_obj._main()

