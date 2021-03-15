#!/usr/bin/env python

from tqdm import tqdm
from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from sample_info import sample_info, rev_primer, proQ
from matplotlib import pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import os
import sys

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.sans-serif'] = "Arial"
rcParams['font.sans-serif'] = "Arial"



def parse_reads(file_list):
    """ Collect reads from fastq and return dict with {sample: {readseq: count}} """
    data = {}
    with open(file_list) as fin:
        for fastq in fin:
            fastq = fastq.strip("\n")
            index = int(fastq[fastq.find("bc10")+4:fastq.find("bc10")+6])
            sample = sample_info[index]["name"]
            data[sample] = {}
            for read in SeqIO.parse(fastq, "fastq"):
                    seq  = str(read.seq)
                    if seq not in data[sample]:
                        data[sample][seq]  = 1
                    else:
                        data[sample][seq] += 1
    return data

def filter_reads(parsed_reads):
    """Filter reads to only include correct length """
    stat_numbers = {fastq: [0.0, 0.0, 0.0, 0.0] for fastq in parsed_reads}
    filt_reads = {fastq: {} for fastq in parsed_reads}
    for fastq in parsed_reads:
        for read in tqdm(parsed_reads[fastq]):
            stat_numbers[fastq][0] += parsed_reads[fastq][read]
            if len(read) == len(proQ):
                stat_numbers[fastq][1] += parsed_reads[fastq][read]
                filt_reads[fastq][read] = parsed_reads[fastq][read]
            elif len(read) > len(proQ):
                stat_numbers[fastq][2] += parsed_reads[fastq][read]
            elif len(read) < len(proQ):
                stat_numbers[fastq][3] += parsed_reads[fastq][read]
    return filt_reads, stat_numbers

def print_read_stats(stats):
    """ Print basic read stats to read_summary.txt """
    labels = ["Type", "Rep", "Norm", "Stat", "Read count"]
    with open("read_stats/read_summary.txt", "w") as fout:
        fout.write("\t".join(labels) + "\n")
        for fastq in stats:
            typ, rep = fastq.split("-")
            line = "{}\t{}\t{}\t{}\t{}\n".format(typ, rep, {}, {}, {})
            fout.write(line.format("raw", "Total", stats[fastq][0]))
            fout.write(line.format("raw", "Correct length", stats[fastq][1]))
            fout.write(line.format("raw", "Longer", stats[fastq][2]))
            fout.write(line.format("raw", "Shorter", stats[fastq][3]))
            fout.write(line.format("%", "Correct length", stats[fastq][1]/stats[fastq][0]*100))
            fout.write(line.format("%", "Longer", stats[fastq][2]/stats[fastq][0]*100))
            fout.write(line.format("%", "Shorter", stats[fastq][3]/stats[fastq][0]*100))
    return

def plot_read_stats():
    """ Create a barplot (pdf) of read stats """
    data = pd.read_csv("read_stats/read_summary.txt", sep="\t")
    ax = sns.barplot("Stat", "Read count", data=data[data["Norm"] == "raw"], hue="Type")
    ax.set(xlabel = "")
    plt.savefig("read_stats/read_plot.pdf")
    plt.close()
    ax = sns.barplot("Stat", "Read count", data=data[data["Norm"] == "%"], hue="Type")
    ax.set(ylabel = "Reads (%)", xlabel = "")
    plt.savefig("read_stats/read_plot_percent.pdf")


def create_out_folder(name):
    """ Create folder for output files """
    out = os.getcwd() + "/" + name
    if os.path.exists(out):
        print(out + ' : exists')
        if os.path.isdir(out):
            print(out + ' : is a directory')
            if input("Use this folder? (Files may be overwritten) [y]") == "y":
                return
            else:
                sys.exit()
    else:
        os.mkdir(out)
    return

def print_parsed_reads(parsed_data):
    """ Save parsed reads as a python dictionary in a seperate file """
    out = os.getcwd() + "/scripts/parsed_data.py"
    if os.path.exists(out):
        print(out + " : exists")
        while os.path.exists(out):
            out = os.getcwd() + "/" + input("Enter new filename: ")
    with open(out, "w") as fout:
        fout.write("reads = {}".format(parsed_data))
    return


create_out_folder("read_stats")
fastq_list = argv[1] # File with list of paths to fastq files to be analyzed
print("Parsing reads...")
reads = parse_reads(fastq_list)
print("Filtering reads..")
reads, stats = filter_reads(reads)
print_read_stats(stats)
plot_read_stats()
print_parsed_reads(reads)




