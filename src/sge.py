from conf import *

class Sge(object):

	def __init__(self, input_dir, max_queue, sge_log):
		self._input_dir = input_dir
		self._max_queue = max_queue
		self._sge_log = sge_log

	def _submit(self):
		"""
		submit all the files in input_dir
		"""
		for file in self._input_dir:
			fastq_file = file
			cmd = "qsub -V -l mem=16G -cwd main.py -o " + self._sge_log
			os.system(cmd)

	def _monitor(self):
		"""
		monitor all the jobs in the queue
		wait untill all the jobs finished
		submit another round
		:return: 
		"""
		cmd = "qstat"
		stat = os.popen(cmd)
		output = stat.read()
		if output == "":
			self._submit()

