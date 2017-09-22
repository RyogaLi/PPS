# some supplimentary functions
# todo remove smorf in ref and align
from conf import *

def combine_fastq(f1, f2):
	"""
	combine r1 and r2
	change the id in r2 to somthing else
	:param f1: read 1
	:param f2: read 2
	:return: 
	"""
	# change all the ids in f2
	new_r2 = f1.split("_R1")[0] + ".fastq"
	with open(f2, "r") as r2, open(f1, "r") as r1:
		with open(new_r2, "w") as new:
			for line in r1:
				if line.startswith("@M00730"):
					line = line.split(" ")
					# print line[1:]
					line[1] = "1"+line[1][1:]
					id = "_".join(line)
					new.write(id)
				else:
					new.write(line)
			for line in r2:
				if line.startswith("@M00730"):
					line = line.split(" ")
					# print line[1:]
					line[1] = "1"+line[1][1:]
					id = "_".join(line)
					new.write(id)
				else:
					new.write(line)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Supplementary functions')
	parser.add_argument('combine',nargs="?",type=str, help='Combine R1 and R2')
	args = parser.parse_args()

	if args.combine:
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
			combine_fastq(fastq_path+r1, fastq_path+r2)
