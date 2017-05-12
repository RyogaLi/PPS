import os
import sys
import logging.config

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
		if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
			command = "bowtie2 -a " + " -x " + self._reference +" -U " + fastq  + " -S " + basename + ".sam " + " > bowtie_DEFAULT_" + basename + ".log"
			os.system(command)
		elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
			command = "bowtie2 -a -p 20 -x ORF_with_pDONR --local --very-sensitive-local -U "+fastq+ "-S "+self._reference + " > bowtie_SENSITIVE_" + basename + ".log"
			os.system(command)
		else:
			command = "ERROR: please provide correct setting (DEFAULT/SENSITIVE)"
			sys.exit(command)
		return command

	def _main(self):
		# init logging
		logging.config.fileConfig("./logging.conf")
		logger = logging.getLogger("alignment")

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
			break



if __name__ == "__main__":
	# get all the names of fastq file
	fastq_path = "../03_PPS_DK/"
	reference = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/ref/ORF_reference_pDONOR"
	alignment_obj = Alignment(reference, fastq_path)
	alignment_obj._main()

