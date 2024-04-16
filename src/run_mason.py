#!/usr/bin/env python3
"""
Simulates short-read Illumina data from amplicon sequences using Mason.
"""

import argparse
from subprocess import call
import os
import sys

from utilities import compute_number_of_reads_to_simulate, get_fasta_length


def check_mason_path(mason_path=None):
    if not os.path.isfile(mason_path):

        raise Exception("Cannot find Mason. Specify an existing file path for " 
                        " Mason.")

def run_mason(fasta, mason_path, output_prefix, read_length, coverage, circ_or_linear, nthreads, seed=None):
    """Simulates short-read data with Mason.

    Args:
        fasta: Path to FASTA file containing amplicon.
        mason_path: Path to Mason executable.
        output_prefix: Prefix for file path to write out simulated read data.
        read_length: Read length for simulated Illumina data (e.g., 150).
        coverage: Target coverage for amplicon.
        circ_or_linear: Whether or not to simulate reads from a circular or
            linear amplicon. (Note used currently.)
        nthreads: Number of threads to utilize.
        seed: Random seed to utilize.

    Returns:
        None. Call Mason to simulate data.
    """
    fasta_length = get_fasta_length(fasta)
    num_reads = str(compute_number_of_reads_to_simulate(coverage, fasta_length, 2*read_length))

    R1 = f"{output_prefix}_R1.fastq.gz"
    R2 = f"{output_prefix}_R2.fastq.gz"
    cmd = "{} -ir {} -n {} --illumina-read-length {} --seq-technology illumina --num-threads {} -o {} -or {}".format(
        mason_path, fasta, num_reads, read_length, nthreads, R1, R2
    )

    if seed:
        cmd+=(" --seed " + str(seed))

    print(cmd)
    call(cmd, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Call Mason to simulate nanopore reads from amplicon genome "
                                                 "structures and background genome regions.")
    parser.add_argument("-o", "--output_prefix", help="Prefix for output filenames.", required=True)
    parser.add_argument("--amplicon_fasta", help="Fasta file of simulated amplicon structure.", required=True)
    parser.add_argument("--amplicon_coverage", type=float, help="Coverage for amplicon region - scale up to simulate "
                                                                "higher amplicon CNs", required=True)

    parser.add_argument("-l", "--read_length", help="Read length", type=int, default=150)
    parser.add_argument("--background_fasta", help="Fasta file of background genome regions. Do not provide if you want"
                                                   " reads from amplicon only")
    parser.add_argument("--background_coverage", type=float, help="Coverage for background region - scale up to "
                                                                  "simulate lower amplicon CNs")
    parser.add_argument("-t", "--num_threads", type=int, help="Number of threads to use in all stages (default=1).",
                        default=1)
    parser.add_argument("--seed", help="Manually seed the pseudo-random number generator for the simulation.")
    parser.add_argument("--mason_path", help="Custom path to Mason executable (mason_simulator).", default='mason_simulator')

    
    # Note: this parameter is not used for now.
    parser.add_argument("--amp_is_linear", action='store_true', help="Set if the simulated amplicon is linear instead"
                                                                     " of circular.")
    args = parser.parse_args()
    if args.background_fasta and not args.background_coverage:
        print("--background_coverage required if --background_fasta is provided!")
        sys.exit(1)

    circular_or_linear = "linear" if args.amp_is_linear else "circular"

    # set the path of the mason read simulation tool
    if not args.mason_path.endswith('mason_simulator') and args.mason_path != 'mason_simulator':
        args.mason_path += "/mason_simulator"

    check_mason_path(args.mason_path)

    # simulate reads from the amplicon
    print("Simulating reads from amplicon...")
    amplicon_prefix = args.output_prefix + "_amplicon_reads"
    background_prefix = args.output_prefix + "_background_reads"
    run_mason(
        args.amplicon_fasta,
        args.mason_path,
        amplicon_prefix,
        args.read_length,
        args.amplicon_coverage,
        circular_or_linear,
        args.num_threads,
        args.seed
    )

    # simulate reads from the background
    if args.background_fasta:
        print("\nSimulating reads from background...")
        background_prefix = args.output_prefix + "_background_reads"
        run_mason(
            args.background_fasta,
            args.mason_path,
            background_prefix,
            args.read_length,
            args.background_coverage,
            "linear",
            args.num_threads,
            args.seed
        )