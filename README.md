# proQ_mutagenesis


Collection of python scripts used to analyze ProQ saturation mutagenesis PacBio data.
Sample info summarized in scripts/sample_info.py

# Generate read stats and mutation distributions

scripts/reads_per_sample.py fastq_list.txt

scripts/check_mutations.py # Requires reads_per_sample.py to be run first

# Generate count table to be used for DESeq2 analyses

scripts/read_counts.py file1.fastq,file2,fastq...
