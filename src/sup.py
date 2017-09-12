from conf import *
from plot import *

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
	# todo compare alignment rate of different reference files
	old_ref_rate = get_alignment_rate("/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_old_ref/")
	combined_ref_rate = get_alignment_rate("/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_combined_ref/")

	plt.plot(range(len(old_ref_rate)), old_ref_rate, '.')
	plt.plot(range(len(combined_ref_rate)), combined_ref_rate, '.')
	plt.title("Compare alignment rate")
	plt.xlabel("plate")
	plt.ylabel("% aligned")
	plt.savefig("alignment_rate.png")