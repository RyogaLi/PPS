from conf import *

def get_alignment_rate(directory):
	"""
	get alignment rate from .log file in directory
	:param directory: 
	:return: 
	"""
	percentages = []
	file_list = os.listdir(directory)
	for dirpath, dirnames, filenames in os.walk(directory):
		for filename in filenames:
			if filename.startswith("bowtie_"):
				line = subprocess.check_output(['tail', '-1', os.path.join(dirpath, filename)])
				per = float(line.split("%")[0])
				percentages.append(per)
	return percentages

if __name__ == "__main__":
	print get_alignment_rate("/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_old_ref/")
