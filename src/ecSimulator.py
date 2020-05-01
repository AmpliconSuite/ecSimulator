#!/usr/bin/env python3
import sys
import os
import copy
import yaml
import bisect
import argparse
import numpy as np
from itertools import groupby
from collections import defaultdict
from intervaltree import Interval, IntervalTree

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
PROG_DIR = SRC_DIR[:SRC_DIR.rindex('/')]
LC_DIR = PROG_DIR + "/low_comp_regions/"

print(SRC_DIR,PROG_DIR,LC_DIR)

class gInterval(object):
    def __init__(self,chrom,start,end,seq):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq
        self.size = abs(self.end - self.start)

    def to_string(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + " | " + str(self.size)

#--------------------------------------------------------------------
#IO METHODS

#read the excludable region database
#return a dictionary of interval trees
def read_excludedRegions(cent_file,map_exc_file):
    excIT = defaultdict(IntervalTree)
    with open(exc_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fields[1],fields[2] = int(fields[1]),int(fields[2])
            excIT[fields[0]].add(Interval(fields[1],fields[2]))

    return excIT


#fasta parsing, derived from brent pedersen
#https://www.biostars.org/p/710/
def readFasta(fname,chrSet):
    seqD = {}
    seqStartIndices = []
    with open(fname) as infile:
        currPos = 0
        faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">"))
        for header in faiter:
            name = header.__next__()[1:].strip()
            if name in chrSet:
                seq = "".join(s.strip() for s in faiter.__next__())
                seqD[name] = seq
                seqStartIndices.append((name,currPos))
                currPos+=len(seq)

    return seqD,seqStartIndices,currPos

def readChromSet(fname):
    return set(line.strip() for line in open(fname))

#---------------------------------------------------------------------
def compute_interval_sizes(num_segs,target_size,min_seg_size = 100000):
    interval_sizes = []
    df = num_segs - 1
    remaining_space = target_size
    while df > 0:
        if (remaining_space - min_seg_size < 2*min_seg_size):
            print("Remaining space too small. Increasing target amplicon size.")
            remaining_space+=(3*min_seg_size)

        loc = np.random.randint(min_seg_size,remaining_space - min_seg_size)
        currSegSize = min(loc, remaining_space - loc)
        interval_sizes.append(currSegSize)
        remaining_space -= currSegSize
        df-=1

    interval_sizes.append(remaining_space)
    return interval_sizes


def compute_interval_regions(interval_sizes, ref_gsize, seqStartInds, seqD, excIT):
    intervals = []
    chrom_snames, chrom_sposns = zip(*seqStartInds)
    for s in interval_sizes:
        foundInt = False
        iters = 0
        while not foundInt and iters < 10000:
            iters+=1
            #pull a random number from the ref.
            spos = np.random.randint(0,ref_gsize - s)
            sposn_index = bisect.bisect_left(chrom_sposns,spos) - 1
            schrom = chrom_snames[sposn_index]
            epos = spos + s
            if epos < chrom_sposns[sposn_index + 1]:
                #lookup the sequence
                normStart = spos - chrom_sposns[sposn_index]
                normEnd = epos - chrom_sposns[sposn_index]
                if not excIT[schrom][normStart:normEnd]:
                    currSeq = seqD[schrom][normStart:normEnd]
                    if not "N" in currSeq:
                        foundInt = True
                        intervals.append(gInterval(schrom,normStart,normEnd,currSeq))

    return intervals



#---------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate focal amplifications such as ecDNA "
                                                 "or virally-mediated amplicons.")
    parser.add_argument("--ref_name", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38"], required=True)
    parser.add_argument("--ref_fasta", help="Reference genome fasta.", required=True)
    parser.add_argument("--mode", help="Type of amplicon to simulate.", choices=["ecDNA","viralAmp"], default="ecDNA")
    parser.add_argument("--config_file", help="Path to config file for run. "
                                              "Will use mode's default params if not provided.")
    parser.add_argument("-o", "--output_prefix", help="Prefix for output filenames.", required=True)
    args = parser.parse_args()

    if not args.config_file:
        args.config_file = SRC_DIR + "/default_config.yaml"

    with open(args.config_file) as f:
        run_config = yaml.safe_load(f)

    #get the list of chromosomes for this sample
    with open(LC_DIR + args.ref_name + "_chroms.txt") as f:
        chrSet = set(f.read().splitlines())

    print("Reading ref")
    #read the reference genome
    seqD, seqStartInds, ref_gsize = readFasta(args.ref_fasta,chrSet)
    if args.mode == "viralAmp":
        #read the viral genome
        pass

    print(chrSet,len(seqD),seqStartInds,ref_gsize)

    #read the low complexity list
    excIT = read_excludedRegions(LC_DIR + args.ref_name + "_wgMapabilityExcludable.bed")

    print("Selecting intervals")
    target_size = run_config["ecdna_target_size"] if args.mode == "ecDNA" else run_config["viral_amp_target_size"]
    interval_sizes = compute_interval_sizes(run_config["num_intervals"],target_size)
    interval_regions = compute_interval_regions(interval_sizes,ref_gsize,seqStartInds,seqD,excIT)
    for i in interval_regions:
        print(i.to_string())














