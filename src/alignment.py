from conf import *
from sge import *
import os
import sys
import logging.config
import shutil

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
			self._align(r1, r2)
		else: # if it's not paired, align R1 and R2 separately
			self._align(self._fastq_file)

if __name__ == "__main__":
	# get all the names of fastq file
# 	fastq_path = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_DK/"
# 	reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/ref/ORF_reference_pDONOR"
	alignment_obj = Alignment(all_reference,fastq, ALIGNMENT_SETTING)
	alignment_obj._main()
	# alignment_obj._gene_count()

