import bisect
import copy
import logging
import math
import sys

from numpy import cumsum
from numpy import random as r


from amplicon_interval import *


# Interval size selection
def compute_interval_sizes(num_segs, target_size, min_interval_size, max_interval_size):
    interval_sizes = []
    df = num_segs - 1
    remaining_space = target_size
    while df > 0:
        if remaining_space - min_interval_size < 2*min_interval_size:
            logging.warning("Remaining space too small. Increasing target amplicon size.")
            remaining_space+=(3*min_interval_size)

        loc = r.randint(min_interval_size,remaining_space - min_interval_size)
        currSegSize = min(loc, remaining_space - loc)
        if currSegSize > max_interval_size:
            currSegSize = max_interval_size
            logging.warning("Interval size exceeded maximumum allowable interval size. \
            Set interval to maximum size instead. This may change the target amplicon size.")
        interval_sizes.append(currSegSize)
        remaining_space-=currSegSize
        df-=1

    interval_sizes.append(remaining_space)
    return interval_sizes


# take a gInterval and make a new gInterval that includes some flanking regions for the interval
def get_padded_intervals(intervals, seqD, flanking_length):
    padded_intervals = []
    for ival in intervals:
        if not ival.is_chromosomal:
            continue

        flanked_start = max(0, ival.start - flanking_length)
        flanked_end = min(ival.end + flanking_length, len(seqD[ival.chrom]) - 1)
        flanked_seq = seqD[ival.chrom][flanked_start:flanked_end + 1]
        padded_intervals.append(gInterval(ival.chrom, flanked_start, flanked_end, flanked_seq, ival.seg_id))

    return padded_intervals

# select the ecDNA intervals
def select_interval_regions(interval_sizes, ref_gsize, seqStartInds, seqD, excIT, used_intervals, sameChrom,
                                origin, overlap_ivald, viralName="", viralSeq=""):
    intervals = []
    if overlap_ivald:
        flat_overlapl = [(x, y.begin, y.end) for x in overlap_ivald.keys() for y in overlap_ivald[x]]
        if len(set([x[0] for x in flat_overlapl])) > 0 and sameChrom:
            logging.warning("multichromosomal is set to false but overlap bed contains regions from multiple "
                            "chromosomes. Overlap bed takes precedence over multichromosomal = False")

    chrom_snames, chrom_sposns = map(list, zip(*seqStartInds))
    chrom_sposns.append(ref_gsize)
    for unadj_s in interval_sizes:
        s = unadj_s - 1
        foundInt = False
        iters = 0
        while not foundInt and iters < 20000:
            iters += 1
            # if there are regions desired to be in the amplicon, use those
            if overlap_ivald:
                elem = flat_overlapl[r.randint(0, len(flat_overlapl))]
                schrom, init_start, init_end = elem
                # check if interval wide enough to pull from within the interval
                elem_len = init_end - init_start
                if elem_len >= s:
                    d = elem_len - s
                    lclip = r.randint(1, d)
                    rclip = d - lclip
                    normStart = init_start + lclip
                    normEnd = init_end - rclip

                # try to flank the interval
                else:
                    d = s - elem_len
                    lpad = r.randint(1, d)
                    rpad = d - lpad
                    normStart = init_start - lpad
                    normEnd = init_end + rpad

                if not used_intervals[schrom].overlaps(normStart, normEnd):
                    currSeq = seqD[schrom][normStart:normEnd+1]
                    if "N" not in currSeq:
                        foundInt = True
                        logging.info("Identified an interval in " + str(iters) + " iterations")
                        seg_id = len(intervals) + 1
                        intervals.append(gInterval(schrom, normStart, normEnd, currSeq, seg_id))
                        if not overlap_ivald:
                            used_intervals[schrom].addi(normStart, normEnd)

                    if excIT[schrom][normStart:normEnd] or "N" in currSeq:
                        logging.warning("Interval {}:{}-{} contained some low-complexity regions".format(
                            schrom, str(normStart), str(normEnd)))

            # else pull a random number from the ref.
            else:
                spos = r.randint(0, ref_gsize - s - 1)
                sposn_index = bisect.bisect_right(chrom_sposns,spos) - 1
                epos = spos + s
                if epos < chrom_sposns[sposn_index + 1]:
                    # lookup the sequence
                    schrom = chrom_snames[sposn_index]
                    if not sameChrom or not intervals or schrom == intervals[0].chrom:
                        normStart = spos - chrom_sposns[sposn_index]
                        normEnd = epos - chrom_sposns[sposn_index]
                        if not excIT[schrom][normStart:normEnd] and not used_intervals[schrom].overlaps(normStart, normEnd):
                            currSeq = seqD[schrom][normStart:normEnd+1]
                            if "N" not in currSeq:
                                foundInt = True
                                logging.info("Identified an interval in " + str(iters) + " iterations")
                                seg_id = len(intervals)+1
                                intervals.append(gInterval(schrom, normStart, normEnd, currSeq, seg_id))
                                used_intervals[schrom].addi(normStart, normEnd)

        if not foundInt:
            errorM = "Could not generate " + str(len(interval_sizes)) + " valid random intervals of size " + str(s) + \
                     " in a reasonable number of iterations!"
            logging.error(errorM)
            sys.stdout.write(errorM + "\n")
            sys.exit(1)

    if viralSeq and viralName:
        vInt = gInterval(viralName, 1, len(viralSeq), viralSeq, len(intervals)+1)
        vInt.preserve = True
        vInt.is_chromosomal = False
        intervals.append(vInt)

    if origin == "chromothripsis" or origin == "two-foldback":
        intervals = sorted(intervals, key=lambda x: (x.chrom, x.start))
        for ind, x in enumerate(intervals):
            x.seg_id = ind + 1

    logging.info("Selected intervals: ")
    for ival in intervals:
        logging.info(str(ival))

    return intervals


def assign_bps(amp_intervals, flankingLength, num_breakpoints):
    bp_intervals = copy.deepcopy(amp_intervals)
    i = 0
    while i < num_breakpoints:
        # get the length of everything
        sizeList = [s.size for s in bp_intervals]
        cumulative_starts = [0, ] + list(cumsum(sizeList))
        origLength = float(cumulative_starts[-1])

        gotBP = False
        while not gotBP:
            # pick a random breakpoint
            loc = r.randint(flankingLength, origLength-flankingLength)
            sposn_index = bisect.bisect(cumulative_starts, loc) - 1
            relStart = loc - cumulative_starts[sposn_index]
            # print(loc, sposn_index, relStart, origLength, cumulative_starts)
            d2E = cumulative_starts[sposn_index+1] - loc

            # check if it's allowable
            if relStart > flankingLength and d2E > flankingLength and not bp_intervals[sposn_index].preserve:
                a, b = bp_intervals[sposn_index].break_interval(relStart)
                bp_intervals = bp_intervals[:sposn_index] + [a, b] + bp_intervals[sposn_index+1:]
                i+=1
                gotBP = True

    # update the ids
    for ind, ival in enumerate(bp_intervals):
        ival.seg_id = ind + 1

    return bp_intervals


def conduct_EC_SV(segL, num_breakpoints, sv_probs, origin):
    intermediates = [copy.deepcopy(segL), ]  # this stores the intermediate structures from each step
    newSegL = copy.deepcopy(segL)
    origUniqueSegs = len(segL)
    if origUniqueSegs == 1:
        return [newSegL, ]

    origLength = float(sum(s.size for s in segL))
    delProb, dupProb, invProb, transProb, fbackProb = sv_probs["del"], sv_probs["dup"], sv_probs["inv"], sv_probs["trans"], sv_probs["fback"]
    safeIds = set([x.seg_id for x in newSegL if x.preserve])

    i = 0
    if origin == "two-foldback":
        rev_ints = copy.deepcopy(newSegL[::-1])
        for x in rev_ints:
            x.direction = -1
            x.seq = str(x.reverse_complement())

        newSegL = newSegL + rev_ints
        intermediates.append(copy.deepcopy(newSegL))
        origLength*=2
        logging.info("initialized two-foldback case")
        i+=1

    while i < int(math.ceil(num_breakpoints/2.0)):
        logging.info(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in newSegL]) + "\n")
        currUniqueSegs = len(set(newSegL))
        currLength = float(sum(s.size for s in newSegL))
        lenDiff = (currLength - origLength)
        logging.debug("lenDiff: " + str(lenDiff))
        meanSegsTogether = int(math.ceil(math.log2(len(newSegL))))

        zeroLen = True
        manipSegs, unManipSegs, manipLen, numSegs, loopedEndI = [], [], 0, 0, 0
        while zeroLen:
            # select random uniform startpoint
            strtP = r.randint(0, len(newSegL))
            # poisson seg nums
            numSegs = r.poisson(meanSegsTogether - 1) + 1
            if numSegs >= len(newSegL):
                numSegs = len(newSegL)
            loopedEndI = (strtP + numSegs) % len(newSegL)
            startI = strtP
            if startI != loopedEndI:
                zeroLen = False
                manipSegs = (copy.deepcopy(newSegL) + copy.deepcopy(newSegL))[strtP:strtP + numSegs]
                manipLen = float(sum(s.size for s in manipSegs))
                unManipSegs = (copy.deepcopy(newSegL) + copy.deepcopy(newSegL))[strtP + numSegs:len(newSegL) + strtP]
                logging.debug("len newSegL, len manipSegs, len unManipSegs: " + str((len(newSegL), len(manipSegs), len(unManipSegs))))

        currSafeIds = set([x.seg_id for x in unManipSegs if x.preserve])
        delAllowed = (currSafeIds == safeIds)

        # deletion
        # print(delAllowed,delProb,(manipLen - lenDiff) / origLength, currUniqueSegs / origUniqueSegs,
        # float(len(manipSegs)) / len(newSegL))
        if delAllowed and r.random() < delProb and (manipLen - lenDiff) / origLength < 0.2 and \
                currUniqueSegs / origUniqueSegs > 0.4 and float(len(manipSegs)) / len(newSegL) <= 0.4:
            newSegL = unManipSegs
            i += 1
            logging.info(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in manipSegs]) + " | " + str(manipLen))
            logging.debug("del")

        # other SV
        else:
            trans, dup, inv = False, False, False

            if r.random() < fbackProb:
                dup = True
                inv = True

            else:
                # is translocated
                if r.random() < transProb:
                    trans = True

                # is duplicated
                if r.random() < dupProb and (currLength + manipLen) / origLength < 1.3:
                    logging.debug("length info: " + str((manipLen, lenDiff, origLength, (currLength + manipLen) / origLength)))
                    dup = True

                # is inverted
                if r.random() < invProb:
                    inv = True

            logging.info(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in manipSegs]) + " | " + str(manipLen))
            logging.debug("trans,dup,inv: " + str((trans, dup, inv)))
            # do inversion now
            if inv:
                manipSegs = manipSegs[::-1]
                for ind, j in enumerate(manipSegs):
                    j.seq = str(j.reverse_complement())
                    j.direction*=-1

                i += 1

            if trans:
                if dup:
                    newSegL = unManipSegs + manipSegs
                    i+=1

                else:
                    newSegL = unManipSegs

                insP = r.randint(0, len(newSegL))
                newSegL[insP:insP] = manipSegs
                i+=1

            else:
                if dup:
                    newSegL[loopedEndI:loopedEndI] = manipSegs
                    i+=1

                else:
                    newSegL = unManipSegs + manipSegs

        intermediates.append(copy.deepcopy(newSegL))
        # logging.debug(
        #     ",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in newSegL]) + "\n")
        # for x in newSegL:
        #     print(str(x), len(x.seq))

    intermediates.append(copy.deepcopy(newSegL))
    logging.info(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in newSegL]) + "\n")
    logging.info("Final size: " + str(sum(s.size for s in newSegL)))
    return intermediates


# def restrictedExponential(mean, minV, maxV):
#     if mean < minV or mean > maxV:
#         sys.stderr.write("Improper args, distribution mean outside bounds\n")
#         sys.exit()
#
#     val = int(round(r.exponential(mean)))
#     while val > maxV or val < minV:
#         val = int(round(r.exponential(mean)))
#
#     return val


# # takes a list of 2-tuples containing numeric values and checks for overlap in the tuples
# def checkOverlappingTuples(tupList, tupToAdd):
#     for i in tupList:
#         if tupToAdd[0] >= i[0] and tupToAdd[0] <= i[1]:
#             return True
#         elif tupToAdd[1] >= i[0] and tupToAdd[1] <= i[1]:
#             return True
#
#         return False
