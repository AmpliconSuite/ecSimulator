#!/usr/bin/env python3

__version__ = "0.4"
__author__ = "Jens Luebeck (jluebeck [ at] ucsd.edu)"

import argparse
import os
import sys
import yaml

if sys.version_info[0] < 3:
    raise Exception("Must use Python 3")

from ecSim_IO import *
from ecSim_sv_gen import *

min_segment_size = 1000
min_interval_size = 10000
max_interval_size = 20000000
flanking_length = 100000
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROG_DIR = SRC_DIR[:SRC_DIR.rindex('/')]
LC_DIR = PROG_DIR + "/low_comp_regions/"
VIR_DIR = PROG_DIR + "/oncoviruses/"


def run_sim(LC_DIR, ref_name, ref_fasta, sim_config, num_amplicons, output_prefix):
    # get the list of chromosomes for this sample
    chrSet = readChromSet(LC_DIR + ref_name + "_chroms.txt")

    # read the low complexity list
    excIT = read_excludedRegions(LC_DIR + ref_name + "_centromere.bed",
                                 LC_DIR + ref_name + "_wgMapabilityExcludable.bed")

    # set the random seed if requested
    if sim_config["random_seed"] is not None:
        r.seed(sim_config["random_seed"])
        logging.info("Setting random seed to " + str(sim_config["random_seed"]))

    # read the reference genome
    print("Reading ref")
    seqD, seqStartInds, ref_gsize = readFasta(ref_fasta, chrSet)
    if sim_config["viral_insertion"]:
        # read the viral genome
        logging.info("Viral insertion amplicon: " + sim_config["viral_strain"])
        viralSeqD, viralSeqStartInds, vir_gsize = readFasta(VIR_DIR + sim_config["viral_strain"],
                                                            sim_config["viral_strain"].rsplit(".")[0])
        if len(viralSeqD) > 1:
            sys.stderr.write("Viral genome has more than one fasta entry. Only first will be used.\n")

        viralName = viralSeqStartInds[0][0]
        viralSeq = viralSeqD[viralName]

    else:
        viralName, viralSeq = "", ""

    # set target size
    isCircular = True
    target_size = sim_config["target_size"]
    if "mean_segment_size" in sim_config:
        mean_seg_size = sim_config["mean_segment_size"]
    else:
        mean_seg_size = 150000.0  # refers to the average distance between breakpoints.

    if sim_config["viral_insertion"] and sim_config["num_breakpoints"] == "auto":
        mean_seg_size = target_size / 13.0

    sameChrom = sim_config["same_chromosome"] is True
    origin = sim_config["origin"].lower()
    allowed_origins = ["episome", "chromothripsis", "tst", "bfb"]
    if not origin in allowed_origins:
        logging.error("Argument for 'origin' in yaml file must be one of " + str(allowed_origins))
        sys.exit(1)

    if sim_config["num_intervals"] == "auto":
        nIntervals = 1  # episomal or tst
        if origin == "chromothripsis":
            nIntervals = 10
            sameChrom = True
        #
        # elif origin == "bfb":
        #     sameChrom = True

    else:
        nIntervals = sim_config["num_intervals"]

    used_intervals = defaultdict(IntervalTree)
    if sim_config["overlap_bed"]:
        overlap_regions = read_overlap_regions(sim_config["overlap_bed"])
    else:
        overlap_regions = None

    # Do the simulations
    logging.info("Amplicons to be simulated: " + str(num_amplicons))
    all_padded_intervals = []
    for x in range(1, num_amplicons + 1):
        amp_num = str(x)
        print("\nGenerating amplicon " + amp_num)
        logging.info('')
        logging.info("--- Generating amplicon " + amp_num + " ---")

        if sim_config["num_breakpoints"] == "auto":
            num_breakpoints = max(0, r.randint(1, int(math.ceil(target_size / mean_seg_size))) - nIntervals)
        else:
            num_breakpoints = int(sim_config["num_breakpoints"])

        interval_sizes = compute_interval_sizes(nIntervals, target_size, min_interval_size, max_interval_size)
        logging.info("Num breakpoints: " + str(num_breakpoints))
        logging.info("Interval sizes: " + str(interval_sizes))
        raw_intervals = select_interval_regions(interval_sizes, ref_gsize, seqStartInds, seqD, excIT,
                            used_intervals, sameChrom, origin, overlap_regions, viralName, viralSeq)

        logging.info("Selected intervals:")
        for y in raw_intervals:
            logging.info(str(y))

        # write a fasta of the background intervals (intervals + padding)
        padded_raw_intervals = get_padded_intervals(raw_intervals, seqD, flanking_length)
        all_padded_intervals.extend(padded_raw_intervals)


        if sim_config["allow_interval_reuse"]:
            used_intervals.clear()

        bp_intervals = assign_bps(raw_intervals, min_segment_size, num_breakpoints)

        logging.info("Doing simulation")
        rearranged_amplicon = conduct_EC_SV(bp_intervals, num_breakpoints, sim_config["sv_probs"])

        cycle_fname = output_prefix + "_amplicon" + amp_num + "_cycles.txt"
        graph_fname = output_prefix + "_amplicon" + amp_num + "_graph.txt"
        write_cycles_file(raw_intervals, bp_intervals, rearranged_amplicon, cycle_fname, isCircular)
        write_bpg_file(bp_intervals, rearranged_amplicon, graph_fname, isCircular)
        write_amplicon_fasta(rearranged_amplicon, output_prefix + "_amplicon" + amp_num + ".fasta", amp_num)

    merged_padded_raw_intervals = merge_gIntervals(all_padded_intervals, seqD)
    merged_intervals_fname = output_prefix + "_background_intervals.fasta"
    write_interval_fasta(merged_intervals_fname, merged_padded_raw_intervals)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate ecDNA genome structures.")
    parser.add_argument("--ref_name", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38"],
                        required=True)
    parser.add_argument("--ref_fasta", help="Reference genome fasta.", required=True)
    # parser.add_argument("--mode", help="Type of amplicon to simulate.", choices=["ecDNA", "BFB"], default="ecDNA")
    parser.add_argument("--config_file", help="Path to config file for run. "
                                              "Will use mode's default params if not provided.")
    parser.add_argument("-o", "--output_prefix", help="Prefix for output filenames.", required=True)
    parser.add_argument("-n", "--num_amplicons", help="Number of amplicons to generate.", type=int, default=1)
    parser.add_argument("-v", "--version", help="Print version and exit.", action='store_true')
    args = parser.parse_args()

    if args.version:
        print("ecSimulator: version " + __version__)
        print("author: " + __author__)

    logging.basicConfig(filename=args.output_prefix + '.log', format='%(asctime)s %(levelname)-8s %(message)s',
                        level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S', filemode='w')

    logging.info("ecSimulator " + __version__)
    logging.info("Using command line args: ")
    logging.info(" ".join(sys.argv[1:]))

    if not args.config_file:
        args.config_file = SRC_DIR + "/default_config.yaml"

    logging.info("Loading YAML config file from " + args.config_file)
    with open(args.config_file) as f:
        sim_config = yaml.safe_load(f)

    run_sim(LC_DIR, args.ref_name, args.ref_fasta, sim_config, args.num_amplicons, args.output_prefix)

    logging.info('')
    logging.info("Finished")
    logging.shutdown()
    sys.exit()
