"""
File for basic utilities.
"""
from math import ceil


def get_fasta_length(fasta):
    tlen = 0
    with open(fasta) as infile:
        for line in infile:
            if line.startswith(">"):
                continue

            tlen+=len(line.rstrip())

    print("Fasta size is " + str(tlen) + "bp")
    return tlen


def compute_number_of_reads_to_simulate(coverage, fasta_len, mean_rl):
    nreads = int(ceil(coverage * fasta_len / mean_rl))
    print("Number of reads to simulate for {}x coverage: {}".format(str(coverage), str(nreads)))
    return nreads
