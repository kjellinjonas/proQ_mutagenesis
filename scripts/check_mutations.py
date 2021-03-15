#!/usr/bin/env python

from tqdm import tqdm
from Bio.Seq import Seq
from sample_info import sample_info, rev_primer, proQ
from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np
import os
import sys
from parsed_data import reads 
import pandas as pd
import seaborn as sns

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.sans-serif'] = "Arial"
rcParams['font.sans-serif'] = "Arial"


### This script requires that reads_per_sample.py have been run first ###

def compare_seq(seq, proQ):
    """ Check total number of mutations and position for an individual seq"""
    tot = 0
    pos = 0
    wt  = []
    mut = []
    if seq != proQ:
        for i in range(len(proQ)):
            if seq[i] != proQ[i]:
                tot += 1
                pos = i # Can only be used for seqs with single mutation
                wt.append(proQ[i])
                mut.append(seq[i])
        if tot > 3:
            tot = 4 
    return tot, pos, wt, mut

def check_mutations(filtered_reads, type, outfolder):
    """ Compare read seq to ProQ nt sequence and quantify mutations and their position """
    mutations = {sample: [0.0, 0.0, 0.0, 0.0, 0.0] for sample in filtered_reads}
    mut_pos   = {sample: [0] * len(proQ) for sample in filtered_reads}
    for sample in filtered_reads:
        for read in tqdm(filtered_reads[sample]):
            tot, pos, wt, mut = compare_seq(read, proQ)
            mutations[sample][tot] += filtered_reads[sample][read]
            if tot == 1:
                mut_pos[sample][pos] += filtered_reads[sample][read]
    print_mutations(mutations, type, outfolder)
    print_mut_pos(mut_pos, type, outfolder)
    return mutations, mut_pos


def print_mut_pos(mutations, type, outfolder):
    """ Create outfile with mutations per position """
    if type == "nt":
        ref = proQ
    if type == "aa" or type == "sm":
        ref = proQ_aa
    with open(outfolder + "/mutated_positions_{}.txt".format(type), "w") as fout:
        fout.write("\t".join([str(x+1) for x in range(len(ref))]) + "\n")
        combined = {"Func": [0] * len(ref),
                    "NonFunc": [0] * len(ref),
                    "Mix": [0] * len(ref),
                    "Wt": [0] * len(ref)}
        for sample in mutations:
            name, rep = sample.split("-")
            for pos in range(len(ref)):
                combined[name][pos] += mutations[sample][pos]
        line = "Func\t" + "\t".join([str(x) for x in combined["Func"]]) + "\n"
        fout.write(line)
        line = "NonFunc\t" + "\t".join([str(x) for x in combined["NonFunc"]]) + "\n"
        fout.write(line)
        line = "Mix\t" + "\t".join([str(x) for x in combined["Mix"]]) + "\n"
        fout.write(line)
        line = "Wt\t" + "\t".join([str(x) for x in combined["Wt"]]) + "\n"
        fout.write(line)
    return

def print_mutations(mutations, type, outfolder):
    """ Create outfile with nr of mutations per sample """
    with open(outfolder + "/nr_mutations_per_sample_{}.txt".format(type), "w") as fout:
        header = ["Type", "Rep", "Norm", "Nr_of_mutations", "Value"]
        fout.write("\t".join(header) + "\n")
        for sample in mutations:
            name, rep = sample.split("-")
            line = "{}\t{}\t{}\t{}\t{}\n".format(name, rep, {}, {}, {})
            fout.write(line.format("raw", "0", mutations[sample][0]))
            fout.write(line.format("raw", "1", mutations[sample][1]))
            fout.write(line.format("raw", "2", mutations[sample][2]))
            fout.write(line.format("raw", "3", mutations[sample][3]))
            fout.write(line.format("raw", "≥4", mutations[sample][4]))
            fout.write(line.format("%", "0", mutations[sample][0]/sum(mutations[sample]) * 100))
            fout.write(line.format("%", "1", mutations[sample][1]/sum(mutations[sample]) * 100))
            fout.write(line.format("%", "2", mutations[sample][2]/sum(mutations[sample]) * 100))
            fout.write(line.format("%", "3", mutations[sample][3]/sum(mutations[sample]) * 100))
            fout.write(line.format("%", "≥4", mutations[sample][4]/sum(mutations[sample]) * 100))
    return

def plot_nr_mutations(type, outfolder):
    """ Create barplot (pdf) with nr of mutations per sample """
    data = pd.read_csv(outfolder +  "/nr_mutations_per_sample_{}.txt".format(type), sep="\t")
    ax = sns.barplot("Nr_of_mutations", "Value", data=data[data["Norm"] == "%"], hue="Type")
    ax.set(ylabel = "Reads (%)", xlabel= "Nr mutations ({})".format(type))
    plt.savefig(outfolder + "/nr_mutations_plot_{}.pdf".format(type))
    plt.close()

def plot_nr_mut_pos(type, outfolder):
    """ Create barplot (pdf) with nr of mutations per conditions (pooled sample counts) """
    data = pd.read_csv(outfolder + "/mutated_positions_{}.txt".format(type), sep = "\t")
    data = data.transpose()
    
    
    data["Position (aa)"] = [x +1 for x in range(229)]

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, sharey=True, gridspec_kw={'hspace': 0}) 
    sns.barplot(x = data["Position (aa)"], y = data["Func"], ax=ax1, color="darkblue") 
    ax1.text(3, 150, str(round(1-len(data[data["Func"] == 0])/len(data), 2) * 100) + " %")
    sns.barplot(x = data["Position (aa)"], y = data["NonFunc"], ax=ax2, color="darkblue") 
    ax2.text(3, 150, str(round(1-len(data[data["NonFunc"] == 0])/len(data), 2) * 100) + " %")
    sns.barplot(x = data["Position (aa)"], y = data["Mix"], ax=ax3, color="darkblue") 
    ax3.text(3, 150, str(round(1-len(data[data["Mix"] == 0])/len(data), 2) * 100) + " %")
    sns.barplot(x = data["Position (aa)"], y = data["Wt"], ax=ax4, color="darkblue") 
    ax4.text(3, 150, str(round(1-len(data[data["Wt"] == 0])/len(data), 2) * 100) + " %")
    plt.xticks(np.arange(0, 230, 20), [x for x in range(0, 230, 20)])
    plt.yticks(np.arange(0, 160, 50), [0, 50, 100, 150])
    plt.savefig(outfolder + "/mutation_distribution_{}.pdf".format(type))
    plt.close()

def check_aa_mutations(filtered_reads, outfolder):
    """ Compare read seq to ProQ aa sequence and quantify mutations and their position """
    mutations = {sample: [0.0, 0.0, 0.0, 0.0, 0.0] for sample in filtered_reads}
    mut_pos   = {sample: [0] * len(proQ_aa) for sample in filtered_reads}
    reads_sm  = {sample: {} for sample in filtered_reads}
    for sample in filtered_reads:
        for read in tqdm(filtered_reads[sample]):
            protein = str(Seq(read).translate())
            tot, pos, wt, mut = compare_seq(protein, proQ_aa)
            mutations[sample][tot] += filtered_reads[sample][read]
            if tot == 1:
                reads_sm[sample][read] = filtered_reads[sample][read]
                mut_pos[sample][pos] += filtered_reads[sample][read]

    print_mutations(mutations, "aa", outfolder)
    print_mut_pos(mut_pos, "aa", outfolder)

    return mutations, mut_pos

def create_out_folder(name):
    """ Create outfolder, warns if folder already exists """
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
    return name




outfolder = create_out_folder("mutation_stats")
proQ_aa = str(Seq(proQ).translate())
reads = reads
nr_mutations, mut_pos = check_mutations(reads, "nt", outfolder)
plot_nr_mutations("nt", outfolder)
nr_aa_mutations, mut_aa_pos = check_aa_mutations(reads, outfolder)
plot_nr_mutations("aa", outfolder)
plot_nr_mut_pos("aa", outfolder)

