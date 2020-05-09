from itertools import groupby
from collections import Counter, defaultdict, namedtuple

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


# fasta parsing, derived from brent pedersen, https://www.biostars.org/p/710/
def readFasta(fname, chrSet=None):
    if chrSet is None:
        chrSet = set()

    seqD = {}
    seqStartIndices = []
    with open(fname) as infile:
        currPos = 0
        faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">"))
        for header in faiter:
            name = header.__next__()[1:].strip()
            if not chrSet or name in chrSet:
                seq = "".join(s.strip() for s in faiter.__next__())
                seqD[name] = seq
                seqStartIndices.append((name,currPos))
                currPos+=len(seq)

    return seqD, seqStartIndices, currPos

# read a list of chromosomes names we can use
def readChromSet(fname):
    return set(line.strip() for line in open(fname))

# --------------------------------------------------------------------
# Output methods

def write_cycles_file(raw_intervals, bp_intervals, amplicon, outname, isCircular):
    with open(outname,'w') as outfile:
        for ival in raw_intervals:
            outline = "\t".join(["Interval", str(ival.seg_id), ival.chrom, str(ival.start), str(ival.end)]) + "\n"
            outfile.write(outline)

        outfile.write("List of cycle segments\n")
        for ival in bp_intervals:
            outline = "\t".join(["Segment", str(ival.seg_id), ival.chrom, str(ival.start), str(ival.end)]) + "\n"
            outfile.write(outline)

        segSeq = ",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in amplicon])
        if not isCircular:
            segSeq = "0+," + segSeq + ",0+"
        outline = "Cycle=1;Copy_count=2.0;Segments=" + segSeq + "\n"
        outfile.write(outline)


def write_bpg_file(bp_intervals, amplicon, outname, isCircular):
    dirToChar = {-1: "-", 1: "+"}

    with open(outname, 'w') as outfile:
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
        for a,b in pairedElems:
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
            outfields+=[posPair, str(2.0*count), str(count), "None", "None"]
            outstring = "\t".join(outfields) + "\n"
            outfile.write(outstring)


def write_amplicon_fasta(amplicon, outname, amp_num):
    with open(outname, 'w') as outfile:
        ampSeq = ">amplicon_" + amp_num + "\n"
        for ival in amplicon:
            ampSeq+=ival.seq

        outfile.write(ampSeq)
