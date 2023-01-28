from collections import defaultdict
import warnings

from intervaltree import IntervalTree

lookup = str.maketrans("ACGT", "TGCA")


# take a list of gIntervals and merge them
def merge_gIntervals(intervals, seqD):
    ivald = defaultdict(IntervalTree)
    merged_intervals = []
    for ival in intervals:
        ivald[ival.chrom].addi(ival.start, ival.end)

    for c, itree in ivald.items():
        itree.merge_overlaps()
        for ival in itree:
            curr_seq = seqD[c][ival.begin:ival.end]
            merged_intervals.append(gInterval(c, ival.begin, ival.end, curr_seq, -1))

    return merged_intervals


class gInterval(object):
    def __init__(self, chrom, start, end, seq, seg_id):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seq = seq
        self.size = abs(self.end - self.start) + 1
        self.seg_id = seg_id
        self.direction = 1
        self.preserve = False
        self.is_chromosomal = True

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + " | " + str(self.size)

    def reverse_complement(self):
        return self.seq.translate(lookup)[::-1]

    def break_interval(self, relative_location):
        if self.preserve:
            warnings.warn("break_interval called on preserved segment")
            return self, gInterval("", 0, 0, "", -1)

        a_s, a_e = self.start, self.start + relative_location - 1
        b_s, b_e = self.start+relative_location, self.end
        a_seq, b_seq = self.seq[:relative_location], self.seq[relative_location:]
        a = gInterval(self.chrom, a_s, a_e, a_seq, -1)
        b = gInterval(self.chrom, b_s, b_e, b_seq, -1)

        return a, b
