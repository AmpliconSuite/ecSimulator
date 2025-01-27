#!/usr/bin/env python3

import argparse
from math import ceil
from subprocess import call
import os
import sys

from utilities import compute_number_of_reads_to_simulate, get_fasta_length

# this script will run the nanosim pipeline to simulate nanopore reads from ecDNA and background genome.
# it assumes nanosim is installed and its scripts are on the system path (true if installed by conda).

model_mean_read_length = 5266.8  # this summarizes the mean read length in the colo320dm training dataset


def extract_nanosim_model(model_path):
    if not os.path.isdir(model_path):
        print("Extracting model files...")
        cmd = "tar -xzf {}".format(model_path + ".tar.gz")
        print(cmd) 
        call(cmd, shell=True)


def run_nanosim(fasta, model_pre, output_pre, coverage, circ_or_linear, read_length_tuple, nthreads, seed=None):
    fasta_length = get_fasta_length(fasta)
    num_reads = str(compute_number_of_reads_to_simulate(coverage, fasta_length, model_mean_read_length))

    cmd = "simulator.py genome -rg {} -c {} -o {} -n {} -dna_type {} -t {} -b guppy --fastq".format(
        fasta, model_pre, output_pre, num_reads, circ_or_linear, str(nthreads))

    if read_length_tuple:
        if read_length_tuple[0]:
            cmd+=(" -min " + str(read_length_tuple[0]))

        if read_length_tuple[1]:
            cmd+=(" -max " + str(read_length_tuple[1]))

        if read_length_tuple[2]:
            cmd+=(" -med " + str(read_length_tuple[2]))

    if seed:
        cmd+=(" --seed " + str(seed))

    print(cmd)
    call(cmd, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Call Nanosim to simulate nanopore reads from amplicon genome "
                                                 "structures and background genome regions.")
    parser.add_argument("-o", "--output_prefix", help="Prefix for output filenames.", required=True)
    parser.add_argument("--amplicon_fasta", help="Fasta file of simulated amplicon structure.", required=True)
    parser.add_argument("--amplicon_coverage", type=float, help="Coverage for amplicon region - scale up to simulate "
                                                                "higher amplicon CNs", required=True)

    parser.add_argument("--background_fasta", help="Fasta file of background genome regions. Do not provide if you want"
                                                   " reads from amplicon only")
    parser.add_argument("--background_coverage", type=float, help="Coverage for background region - scale up to "
                                                                  "simulate lower amplicon CNs")
    parser.add_argument("-t", "--num_threads", type=int, help="Number of threads to use in all stages (default=1).",
                        default=1)
    parser.add_argument("--seed", help="Manually seed the pseudo-random number generator for the simulation.")
    parser.add_argument("-min", "--min_len", type=int, help="Minimum length for simulated reads (default=50)",
                        default=None)
    parser.add_argument("-max", "--max_len", type=int, help="Maximum length for simulated reads (default Infinity)",
                        default=None)
    parser.add_argument("-med", "--median_len", type=int, help="Median lenght for simulated reads (will override the "
                                                               "values used in the training set).", default=None)
    parser.add_argument("--amp_is_linear", action='store_true', help="Set if the simulated amplicon is linear instead"
                                                                     " of circular.")
    args = parser.parse_args()
    if args.background_fasta and not args.background_coverage:
        print("--background_coverage required if --background_fasta is provided!")
        sys.exit(1)

    read_length_tuple = (args.min_len, args.max_len, args.median_len)
    if args.median_len:
        model_mean_read_length = args.median_len  # this is not perfect, because nanosim does not allow for control over
        # the mean read length, rather the median read length. It appears those two are fairly close after inspecting
        # the nanosim source code.

    circular_or_linear = "linear" if args.amp_is_linear else "circular"

    # set the nanopore read simulation model
    src_path = os.path.dirname(os.path.realpath(__file__))
    model_path = src_path.rsplit("/src")[0] + "/nanosim_models/promethion_10_4_1_guppy_model"
    model_prefix = "/colo320dm_promethion_r10_4_1"
    extract_nanosim_model(model_path)
    model_path+=model_prefix

    # simulate reads from the amplicon
    print("Simulating reads from amplicon...")
    amplicon_prefix = args.output_prefix + "_amplicon_reads"
    run_nanosim(args.amplicon_fasta, model_path, amplicon_prefix, args.amplicon_coverage, circular_or_linear,
                read_length_tuple, args.num_threads, args.seed)

    # simulate reads from the background
    if args.background_fasta:
        print("\nSimulating reads from background...")
        background_prefix = args.output_prefix + "_background_reads"
        run_nanosim(args.background_fasta, model_path, background_prefix, args.background_coverage, "linear",
                    read_length_tuple, args.num_threads, args.seed)

