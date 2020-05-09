#!/usr/bin/env python

__version__ = "0.2"
__author__ = "Jens Luebeck"

import os
import yaml
import argparse

from ecSim_IO import *
from ecSim_sv_gen import *

min_segment_size = 1000
min_interval_size = 10000
max_interval_size = 20000000
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROG_DIR = SRC_DIR[:SRC_DIR.rindex('/')]
LC_DIR = PROG_DIR + "/low_comp_regions/"
VIR_DIR = PROG_DIR + "/oncoviruses/"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate focal amplifications such as ecDNA "
                                                 "or virally-mediated amplicons.")
    parser.add_argument("--ref_name", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38"], required=True)
    parser.add_argument("--ref_fasta", help="Reference genome fasta.", required=True)
    parser.add_argument("--mode", help="Type of amplicon to simulate.", choices=["ecDNA", "BFB"], default="ecDNA")
    parser.add_argument("--config_file", help="Path to config file for run. "
                                              "Will use mode's default params if not provided.")
    parser.add_argument("-o", "--output_prefix", help="Prefix for output filenames.", required=True)
    parser.add_argument("-n", "--num_amplicons", help="Number of amplicons to generate.", type=int, default=1)
    parser.add_argument("-v", "--version", help="Print version and exit.", action='store_true')
    args = parser.parse_args()

    if args.version:
        print("ecSimulator: version " + __version__)
        print("author: " + __author__)

    logging.basicConfig(filename=args.output_prefix + '.log',
        format = '%(asctime)s %(levelname)-8s %(message)s',
        level = logging.DEBUG,
        datefmt = '%Y-%m-%d %H:%M:%S',
        filemode = 'w')

    logging.info("ecSimulator " + __version__)
    logging.info("Using command line args: ")
    logging.info(" ".join(sys.argv[1:]))

    if not args.config_file:
        args.config_file = SRC_DIR + "/default_config.yaml"

    logging.info("Loading YAML config file from " + args.config_file)
    with open(args.config_file) as f:
        sim_config = yaml.safe_load(f)

    # get the list of chromosomes for this sample
    chrSet = readChromSet(LC_DIR + args.ref_name + "_chroms.txt")

    # read the low complexity list
    excIT = read_excludedRegions(LC_DIR + args.ref_name + "_centromere.bed",
                                 LC_DIR + args.ref_name + "_wgMapabilityExcludable.bed")

    # set the random seed if requested
    if sim_config["random_seed"] is not None:
        r.seed(sim_config["random_seed"])
        logging.info("Setting random seed to " + str(sim_config["random_seed"]))

    # read the reference genome
    print("Reading ref")
    seqD, seqStartInds, ref_gsize = readFasta(args.ref_fasta, chrSet)
    if sim_config["viral_insertion"]:
        # read the viral genome
        logging.info("Viral insertion amplicon: " + sim_config["viral_strain"])
        viralSeqD, viralSeqStartInds, vir_gsize = readFasta(VIR_DIR + sim_config["viral_strain"])
        if len(viralSeqD) > 1:
            sys.stderr.write("Viral genome has more than one fasta entry. Only first will be used.\n")

        viralName = viralSeqStartInds[0][0]
        viralSeq = viralSeqD[viralName]

    else:
        viralName, viralSeq = "", ""

    # set target size
    if args.mode == "ecDNA":
        granularity = 150000.0
        isCircular = True
        if not sim_config["viral_insertion"]:
            target_size = sim_config["ecdna_target_size"]
        else:
            target_size = sim_config["viral_amp_target_size"]
            granularity = target_size/13.0

        nIntervals = sim_config["num_intervals"]

    elif args.mode == "BFB":
        granularity = 2500000.0
        isCircular = False
        if sim_config["viral_insertion"]:
            warnMess = "BFB mode does not support viral integration. Ignoring viral insertion flag."
            sys.stderr.write(warnMess + "\n")
            logging.warning(warnMess)

        if sim_config["num_intervals"] > 1:
            warnMess = "BFB mode uses one interval only. YAML file specified " + str(sim_config["num_intervals"]) + \
                       " intervals. Setting num_intervals to 1."
            logging.warning(warnMess)

        target_size = sim_config["bfb_target_size"]
        nIntervals = 1

    else:
        isCircular = False
        granularity = 999999999.0
        sys.stderr.write("Unspecified mode\n")
        sys.exit(1)

    # check if user specified number of breakpoints
    if sim_config["num_breakpoints"] != "auto":
        num_breakpoints = int(sim_config["num_breakpoints"])

    num_amplicons = args.num_amplicons if args.num_amplicons > 1 else sim_config["num_amplicons"]
    used_intervals = defaultdict(IntervalTree)

    # Do the simulations
    logging.info("Amplicons to be simulated: " + str(num_amplicons))
    for x in range(1, num_amplicons + 1):
        amp_num = str(x)
        print("\nGenerating amplicon " + amp_num)
        logging.info('')
        logging.info("--- Generating amplicon " + amp_num + " ---")

        if sim_config["num_breakpoints"] == "auto":
            num_breakpoints = r.randint(1, int(math.ceil(target_size/granularity)))

        interval_sizes = compute_interval_sizes(nIntervals, target_size, min_interval_size, max_interval_size)

        logging.info("Num breakpoints: " + str(num_breakpoints))
        logging.info("Interval sizes: " + str(interval_sizes))

        if sim_config["viral_insertion"]:
            viralName = sim_config["viral_strain"].rsplit(".")[0]

        raw_intervals = compute_ec_interval_regions(interval_sizes, ref_gsize, seqStartInds, seqD, excIT,
                                                    used_intervals, viralName, viralSeq)

        if sim_config["allow_overlapping_intervals"]:
            used_intervals.clear()

        bp_intervals = assign_bps(raw_intervals, min_segment_size, num_breakpoints)

        logging.info("Doing simulation")
        if args.mode == "ecDNA":
            rearranged_amplicon = conduct_EC_SV(bp_intervals, num_breakpoints, sim_config["sv_probs"])

        cycleFname = args.output_prefix + "_amplicon" + amp_num + "_cycles.txt"
        graphFname = args.output_prefix + "_amplicon" + amp_num + "_graph.txt"
        write_cycles_file(raw_intervals, bp_intervals, rearranged_amplicon, cycleFname, isCircular)
        write_bpg_file(bp_intervals,rearranged_amplicon, graphFname, isCircular)
        write_amplicon_fasta(rearranged_amplicon,args.output_prefix + "_amplicon" + amp_num + ".fasta", amp_num)

    logging.info('')
    logging.info("Finished")
    logging.shutdown()
    sys.exit()
