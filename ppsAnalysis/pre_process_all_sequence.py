import numpy as np
import os

def process_seq_file(seq_file):
    """
    remove all the dashes in the gene name
    remove SUP sequence if the same gene is in HIP
    save the duplicated SUP sequence to another file
    :param seq_file: all_sequence.txt
    :return:
    """
    HIP = {}
    PROTGEN = {}
    SGD = {}
    with open(seq_file, "r") as seq:

        # store sequences into dictionary based on their class
        for line in seq:

            if "orf_name" in line: continue
            line = line.split()

            if "-" in line[0]:
                line[0] = "".join(line[0].split("-"))
            if line[1] == "HIP":
                HIP[line[0]] = line[2]
            elif line[1] == "SGD":
                SGD[line[0]] = line[2]
            elif line[1] == "PROTGEN":
                PROTGEN[line[0]] = line[2]

    combined = HIP
    for key in PROTGEN.keys():
        if key in HIP.keys():
            continue
        else:
            combined[key] = PROTGEN[key]
    for key in SGD.keys():
        if key in HIP.keys():
            continue
        else:
            combined[key] = SGD[key]
    return combined

def make_fasta(input_dict, output_fasta):
    """
    make a reference (fasta file)
    """
    with open(output_fasta, "w") as fa:
        for key in input_dict.keys():
            fa.write(">"+key+"\n")
            fa.write(input_dict[key]+"\n")

def sep_sequences(fasta):
    """
    sep seq into HIP and SUP
    :param fasta:
    :return:
    """
    HIP = {}
    other = {}
    with open(seq_file, "r") as seq:
        # store sequences into dictionary based on their class
        for line in seq:
            if "orf_name" in line: continue
            line = line.split()
            if "-" in line[0]:
                line[0] = "".join(line[0].split("-"))
            if line[1] == "HIP":
                HIP[line[0]] = line[2]
            else:
                other[line[0]] = line[2]
    return HIP, other

if __name__ == "__main__":
    seq_file = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/02_data/all_sequence.txt"
    # combined_seqs = process_seq_file(seq_file)
    # make_fasta(combined_seqs, "./ORF_ref.fasta")
    HIP, other = sep_sequences(seq_file)
    make_fasta(HIP, "./ORF_hip.fasta")
    make_fasta(other, "./ORF_other.fasta")
