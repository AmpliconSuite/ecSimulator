from collections import Counter, defaultdict, namedtuple
from datetime import datetime
from itertools import groupby
import os

from intervaltree import Interval, IntervalTree


Bpoint = namedtuple("Bpoint", "chrom pos")

# --------------------------------------------------------------------
# Input methods


# read the excludable region database
# return a dictionary of interval trees
def read_excludedRegions(cent_file, map_exc_file):
    excIT = defaultdict(IntervalTree)
    for f in [cent_file, map_exc_file]:
        with open(f) as infile:
            for line in infile:
                fields = line.rstrip().rsplit("\t")
                fields[1], fields[2] = int(fields[1]), int(fields[2])
                excIT[fields[0]].add(Interval(fields[1], fields[2]))

    return excIT


# Reads a bed file into an intervaltree dict
def read_overlap_regions(reqRegBedF):
    reqRegIvalD = defaultdict(IntervalTree)
    with open(reqRegBedF) as infile:
        for line in infile:
            line = line.rstrip()
            if line:
                fields = line.rsplit("\t")
                reqRegIvalD[fields[0]].addi(int(fields[1]), int(fields[2]))

    return reqRegIvalD


# fasta parsing, derived from brent pedersen, https://www.biostars.org/p/710/
def readFasta(fname, ref_name, chrSet=None):
    if chrSet is None:
        pre = "chr" if ref_name != "GRCh37" else ""
        chrSet = set([pre + str(x) for x in range(1, 23)] + [pre + "X", pre + "Y"])

    seqD = {}
    seqStartIndices = []
    with open(fname) as infile:
        currPos = 0
        faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">"))
        for header in faiter:
            name = next(header)[1:].rstrip().rsplit()[0]
            seq = "".join(s.rstrip() for s in next(faiter))
            if name in chrSet:
                # print(name)
                seqD[name] = seq
                seqStartIndices.append((name, currPos))
                currPos+=len(seq)

    return seqD, seqStartIndices, currPos

# read a list of chromosomes names we can use
def readChromSet(fname):
    return set(line.strip() for line in open(fname))


# --------------------------------------------------------------------
# Output methods

def write_outputs(output_prefix, amp_num, raw_intervals, bp_intervals, all_amplicons, isCircular):
    final_amplicon = all_amplicons[-1]
    cycle_fname = output_prefix + "_amplicon" + amp_num + "_cycles.txt"
    graph_fname = output_prefix + "_amplicon" + amp_num + "_graph.txt"
    write_dt = datetime.now().strftime("%B %d, %Y %H:%M:%S")
    write_cycles_file(raw_intervals, bp_intervals, final_amplicon, cycle_fname, isCircular, write_dt)
    write_bpg_file(bp_intervals, final_amplicon, graph_fname, isCircular, write_dt)
    write_amplicon_fasta(final_amplicon, output_prefix + "_amplicon" + amp_num + ".fasta", amp_num)

    # make a directory for the intermediate structures
    intermed_outdir = os.path.dirname(output_prefix) + "/intermediate_structures/"
    if not os.path.exists(intermed_outdir):
        os.makedirs(intermed_outdir)

    intermed_basepre = intermed_outdir + os.path.basename(output_prefix)
    for ind, intermed_amplicon in enumerate(all_amplicons[1:-1]):  # index 0 is reference genome segs, index -1 is the final genome
        intermed_prefix = intermed_basepre + "_intermediate" + str(ind+1)
        cycle_fname = intermed_prefix + "_amplicon" + amp_num + "_cycles.txt"
        graph_fname = intermed_prefix + "_amplicon" + amp_num + "_graph.txt"
        write_cycles_file(raw_intervals, bp_intervals, intermed_amplicon, cycle_fname, isCircular, write_dt)
        write_bpg_file(bp_intervals, intermed_amplicon, graph_fname, isCircular, write_dt)
        write_amplicon_fasta(intermed_amplicon, output_prefix + "_amplicon" + amp_num + ".fasta", amp_num)


def write_cycles_file(raw_intervals, bp_intervals, amplicon, outname, isCircular, write_dt):
    with open(outname,'w') as outfile:
        outfile.write("# ecSimulator " + write_dt + "\n")
        for ival in raw_intervals:
            outline = "\t".join(["Interval", str(ival.seg_id), ival.chrom, str(ival.start), str(ival.end)]) + "\n"
            outfile.write(outline)

        outfile.write("List of cycle segments\n")
        for ival in bp_intervals:
            outline = "\t".join(["Segment", str(ival.seg_id), ival.chrom, str(ival.start), str(ival.end)]) + "\n"
            outfile.write(outline)

        segSeq = ",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in amplicon])
        ampl = str(int(sum([x.size for x in amplicon])))
        if not isCircular:
            segSeq = "0+," + segSeq + ",0+"
        outline = "Cycle=1;Copy_count=1.0;Length={};Segments={}\n".format(ampl, segSeq)
        outfile.write(outline)


def write_bpg_file(bp_intervals, amplicon, outname, isCircular, write_dt):
    dirToChar = {-1: "-", 1: "+"}
    with open(outname, 'w') as outfile:
        outfile.write("# ecSimulator " + write_dt + "\n")
        outfile.write("SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, "
                      "NumberReadsMapped\n")

        segCounts = Counter([i.seg_id for i in amplicon])
        for ival in bp_intervals:
            startString = ival.chrom + ":" + str(ival.start) + "-"
            endString = ival.chrom + ":" + str(ival.end) + "+"
            occ = segCounts[ival.seg_id]
            outline = "\t".join(["sequence", startString, endString, str(2.0*occ), str(occ), str(ival.size), "0"]) + "\n"
            outfile.write(outline)

        outfile.write("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, "
                      "HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n")

        # construct the set of edges, such that reverse oriented breakpoint and forward oriented breakpoint are the same
        edgeCounts = defaultdict(int)
        pairedElems = zip(amplicon, amplicon[1:]) if not isCircular else zip(amplicon, amplicon[1:] + [amplicon[0]])
        for a, b in pairedElems:
            left_end = a.end if a.direction == 1 else a.start
            right_end = b.start if b.direction == 1 else b.end
            bp1 = Bpoint(a.chrom, left_end)
            bp2 = Bpoint(b.chrom, right_end)

            if bp1 <= bp2:
                opair = (dirToChar[a.direction], dirToChar[b.direction * -1])
                edgeCounts[(bp1, bp2, opair)]+=1

            else:
                opair = (dirToChar[b.direction * -1], dirToChar[a.direction])
                edgeCounts[(bp2, bp1, opair)]+=1

        for e,count in edgeCounts.items():
            # concordant edge
            if e[1].pos - e[0].pos == 1 and e[1].chrom == e[0].chrom:
                outfields = ["concordant", ]

            # discordant edge
            else:
                outfields = ["discordant", ]

            posString1 = e[0].chrom + ":" + str(e[0].pos) + e[2][0]
            posString2 = e[1].chrom + ":" + str(e[1].pos) + e[2][1]
            posPair = posString1 + "->" + posString2 + "\t"
            outfields+=[posPair, str(float(count)), str(count), "None", "None"]
            outstring = "\t".join(outfields) + "\n"
            outfile.write(outstring)


def write_amplicon_fasta(amplicon, outname, amp_num):
    with open(outname, 'w') as outfile:
        ampSeq = ">amplicon_" + amp_num + "\n"
        for ival in amplicon:
            ampSeq+=ival.seq

        outfile.write(ampSeq)


def write_interval_fasta(outname, intervals):
    with open(outname, 'w') as outfile:
        outputSeq = ""
        for ind, ival in enumerate(intervals):
            ivalEntry = ">" + str(ival).rsplit(" | ")[0] + "\n"
            ivalEntry+=ival.seq
            ivalEntry+="\n"
            outputSeq+=ivalEntry

        outfile.write(ivalEntry)
