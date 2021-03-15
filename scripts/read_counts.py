#!/usr/bin/env python

from tqdm import tqdm
from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from sample_info import sample_info, rev_primer, proQ


def parse_reads(in_file):
    """Collect and collapse reads from fastq"""
    reads = {}
    for read in SeqIO.parse(in_file, "fastq"):
            seq  = str(read.seq)
            if seq not in reads:
                reads[seq]  = 1
            else:
                reads[seq] += 1
    return reads

def filter_reads(parsed_reads):
    """Filter reads and keep only the ones of correct length"""
    stat_label = ["tot", "correct_length"]
    stat_numbers = [0.0, 0.0, 0.0, 0.0]
    filt_reads = {}
    for read in tqdm(parsed_reads):
        stat_numbers[0] += parsed_reads[read]
        if len(read) == len(proQ):
            stat_numbers[1] += parsed_reads[read]
            filt_reads[read] = parsed_reads[read]

    return filt_reads, stat_label, stat_numbers

def compare_seq(seq, proQ):
    """ Check total number of mutations and position for an individual seq"""
    if seq == proQ: return"proQ"
    else:
        tot = 0
        for i in range(len(proQ)):
            if seq[i] != proQ[i]:
                tot += 1
                pos = i
                wt_aa = proQ[i]
                mut_aa = seq[i]
                read_name = "proQ_{}{}{}".format(wt_aa, i+1, mut_aa)
        if tot > 1: read_name = "filter"

    return read_name


def read_counts(reads, proQ_aa, read_table, fastq):
    """ Count reads """
    read_table = read_table
    total_counts = 0
    kept_counts = 0
    for read in tqdm(reads):
        total_counts += reads[read]
        protein = str(Seq(read).translate())
        read_name = compare_seq(protein, proQ_aa)
        if read_name != "filter":
            kept_counts += reads[read]
            read_table[read_name][fastq] += reads[read]
    print("Counts for {}:\n Reads left after filtering: {}\n Reads kept after filtering multiple mutations: {}". format(fastq, total_counts, kept_counts))
    return read_table

aa_dict = {"Alanine": "A",
           "Arginine": "R",
           "Asparagine": "N",
           "Aspartic acid": "D",
           "Cysteine": "C",
           "Glutamine": "Q",
           "Glutamic acid": "E",
           "Glycine": "G",
           "Histidine": "H",
           "Isoleucine": "I",
           "Leucine": "L",
           "Lysine":  "K",
           "Methionine": "M",
           "Phenylalanine": "F",
           "Proline": "P",
           "Serine": "S",
           "Threonine": "T",
           "Tryptophan": "W",
           "Tyrosine": "Y",
           "Valine": "V",
           "Stop": "*"}


# Create a dictionary with each possible aa mutations for each position in ProQ
proQ_aa = str(Seq(proQ).translate())
possible_reads = {"proQ": 0}
for index in range(len(proQ_aa)):
    for res in aa_dict.values():
        if proQ_aa[index] == res: continue
        mutation_name = "proQ_{}{}{}".format(proQ_aa[index], index+1, res)
        possible_reads[mutation_name] = 0

# Parse and count reads from trimmed fastq files
in_fastqs = argv[1] # , separated list of fastq
read_table = {entry: {in_fastq: 0 for in_fastq in in_fastqs.split(",")} for entry in possible_reads}
for fastq in in_fastqs.split(","):
    print("Parsing reads for {}..".format(fastq))
    reads = parse_reads(fastq)
    print("Filtering reads for {}..".format(fastq))
    reads, header, data = filter_reads(reads)
    print("Counting reads..")
    read_table = read_counts(reads, proQ_aa, read_table, fastq)

# Create count table to be used for DeSeq2 analyses
with open("read_counts.txt", "w") as fout:
    files =  in_fastqs.split(",")
    fout.write("\t".join(files) + "\n")
    for entry in read_table:
        counts = [read_table[entry][x] for x in files]
        if sum(counts) == 0: continue
        counts = [str(x) for x in counts]
        fout.write(entry + "\t" + "\t".join(counts) + "\n")



