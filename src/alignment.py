from conf import *
from sge import *
import os
import sys
import logging.config
import shutil

# todo change all the fastq path to fastq file as it is for one file
# todo r1 and r2 should be separated in sge file

class Alignment(object):

	def __init__(self, all_reference, fastq_file, setting):
		self._reference = all_reference
		# self._sub_reference = subset_reference
		self._fastq_file= fastq_file
		self._setting = setting

	def _align(self, r1, r2=None):
		self._basename = os.path.basename(r1).split(".")[0]
		# create a dir for each alignment
		if os.path.exists(output+self._basename):
			shutil.rmtree(output+self._basename)
		os.makedirs(output+self._basename)
		os.chdir(output+self._basename)

		if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
			if r2 != None: # it's paired
				command = "bowtie2 -a " + " -x " + self._reference +" -1 " + r1 + " -2 " + r2 + " -S " + self._basename + ".sam " + "2> bowtie_DEFAULT_" + self._basename + ".log"
				os.system(command)
			else:
				command = "bowtie2 -a " + " -x " + self._reference + " -U " + r1 + " -S " + self._basename + ".sam " + "2> bowtie_DEFAULT_" + self._basename + ".log"
				os.system(command)
		elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
			if r2 != None:
				command = "bowtie2 -a -p 20 -x " + self._reference +" --local --very-sensitive-local -1 " + r1 + " -2 " + r2 + " -S " + self._basename + ".sam " + "2> bowtie_SENSITIVE_" + self._basename + ".log"
				os.system(command)
			else:
				command = "bowtie2 -a -p 20 -x " + self._reference + " --local --very-sensitive-local -U " + r1 + " -S " + self._basename + ".sam " + "2> bowtie_SENSITIVE_" + self._basename + ".log"
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

	def _main(self):

		# init logging
		logging.config.fileConfig("./src/logging.conf")
		logger = logging.getLogger("alignment")
		# create log directory

		# IF IT'S PAIRED (PAIRED == TRUE)
		# get corresponding R1 and R2 files
		if PAIRED == True:
			# fastq files list
			r1_files = []
			r2_files = []
			for file in os.listdir(fastq_path):
				if file.endswith(".fastq"):
					# fetch R1 in file name and add to r1
					if "R1" in file:
						r1_files.append(file)
					elif "R2" in file:
						r2_files.append(file)

			for r1 in r1_files:
				# find corresponding r2
				identifier = os.path.basename(r1).split("R1")[0]
				r2 = [i for i in r2_files if identifier in i][0]
				logger.info("started aligning %s and %s", os.path.basename(r1), os.path.basename(r2))
				command = self._align(fastq_path+r1, fastq_path+r2)
				logger.info(command)
				logger.info("alignment finished for %s and %s", os.path.basename(r1), os.path.basename(r2))
				# self._gene_count()
		else: # if it's not paired, align R1 and R2 separately
			# for file in os.listdir(fastq_path):
			if self._fastq_file.endswith(".fastq"):
				logger.info("started aligning %s ", os.path.basename(self._fastq_file))
				command = self._align(fastq_path + self._fastq_file)
				logger.info(command)
				logger.info("alignment finished for %s ", os.path.basename(self._fastq_file))

if __name__ == "__main__":
	# get all the names of fastq file
# 	fastq_path = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_DK/"
# 	reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/ref/ORF_reference_pDONOR"
	alignment_obj = Alignment(all_reference,fastq_path, ALIGNMENT_SETTING)
	alignment_obj._main()
	# alignment_obj._gene_count()

