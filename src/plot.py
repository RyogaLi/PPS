from __future__ import division
from conf import *

def plot_readdepth_genecount(output_dir):
	dirnames = os.listdir(output_dir)
	for dir in dirnames:
		rd = []
		gc = []
		if not os.path.exists(output_dir+dir+"/read_depth.txt"): continue
		os.chdir(output_dir+dir)
		read_depth = np.loadtxt("./read_depth.txt", dtype="str")
		gene_count = np.loadtxt("./gene_count.txt", dtype="str")
		for i in read_depth:
			find = gene_count[np.where(gene_count[:, 0] == i[0])]
			if len(find) != 0 and "gene_name" not in find and float(find[0][1]) < 500:
				rd.append(float(i[1]))
				gc.append(float(find[0][1]))
		corr = pearsonr(rd, gc)[0]
		print corr
		plt.plot(rd, gc, ".", label='The red data')
		plt.xlabel("read depth")
		plt.ylabel("gene count")
		plt.savefig("read_depth_vs_gene_count.png")
		plt.close()

def plot_recover_rate(output_dir):
	dirnames = os.listdir(output_dir)
	ref_dict = {}
	recover_rate = []
	for dir in dirnames:
		vcf = output_dir + dir + "/"+dir+"_sorted.raw.vcf"

		if not os.path.exists(vcf): continue
		variant_call_dict = {}

		with open(vcf, "r") as input_file:
			for line in input_file:
				id_line = re.search("<ID=(.+?),length=(.+?)>", line)
				if id_line:
					ref_dict[id_line.group(1)] = int(id_line.group(2))
				if "#" not in line:
					line = line.split("\t")
					gene_name = line[0]
					if gene_name not in variant_call_dict.keys():
						variant_call_dict[gene_name] = 1
					else:
						variant_call_dict[gene_name]+=1
			recover_rate.append(len(variant_call_dict.keys())/len(ref_dict.keys()))
	plt.plot(range(len(recover_rate)), recover_rate, ".")
	plt.title("% covered in each plate")
	plt.xlabel("plate")
	plt.ylabel("% covered")
	plt.savefig("percent_recovered.png")


def plot_recover(f):
	"""
	plot %aligned VS %recovered
	:param f: 
	:return: 
	"""
	recovered = {}
	with open(f, "r") as f:
		for line in f:
			if "plate_name" in line: continue
			line = line.split(",")
			recovered[line[1]] = float(line[4])
	# if the gene is aligned twice, take the one with higher value
	temp = {} # dict with unique gene names
	for key in recovered.keys():
		if key not in temp:
			temp[key]=recovered[key]
		else:
			if recovered[key]>temp[key]:
				temp[key] = recovered[key]
	total = len(temp.values())

	out = []
	aligned = temp.values()
	aligned.sort()
	for value in aligned:
		recovered_percent = len([i for i in aligned if i >value])/total
		out.append(recovered_percent)


	plt.plot(aligned, out, "-")

	plt.xlabel("percent aligned")
	plt.ylabel("percent recovered")
	plt.savefig("./aligned_recovered_plot.png")



if __name__ == "__main__":
	# plot_readdepth_genecount(output)
	plot_recover("./full_covered_gene.csv")