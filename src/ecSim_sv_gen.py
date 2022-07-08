import sys
import copy
import math
import bisect
import logging

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


# select the ecDNA intervals
def compute_ec_interval_regions(interval_sizes, ref_gsize, seqStartInds, seqD, excIT, used_intervals, sameChrom,
                                origin, viralName="", viralSeq=""):
    intervals = []
    chrom_snames, chrom_sposns = map(list, zip(*seqStartInds))
    chrom_sposns.append(ref_gsize)
    for s in interval_sizes:
        foundInt = False
        iters = 0
        while not foundInt and iters < 10000:
            iters += 1
            # pull a random number from the ref.
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
                        currSeq = seqD[schrom][normStart:normEnd]
                        if "N" not in currSeq:
                            foundInt = True
                            logging.info("Identified an interval in " + str(iters) + " iterations")
                            seg_id = len(intervals)+1
                            intervals.append(gInterval(schrom, normStart + 1, normEnd + 1, currSeq, seg_id))
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
        intervals.append(vInt)

    if origin == "chromothripsis":
        intervals = sorted(intervals, key=lambda x: (x.chrom, x.start))
        for ind, x in enumerate(intervals):
            x.seg_id = ind + 1

    elif origin == "tst":
        intervals = sorted(intervals, key=lambda x: (x.chrom, x.start))
        for ind, x in enumerate(intervals):
            x.seg_id = ind + 1
        rev_ints = copy.deepcopy(intervals)[::-1]
        for x in rev_ints:
            x.direction = -1
            x.seq = str(x.reverse_complement())

        intervals = intervals + rev_ints

    logging.info("Intervals: ")
    for ival in intervals:
        logging.info(ival.to_string())

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


def conduct_EC_SV(segL, num_breakpoints, sv_probs):
    newSegL = copy.deepcopy(segL)
    origUniqueSegs = len(segL)
    if origUniqueSegs == 1:
        return newSegL

    origLength = float(sum(s.size for s in segL))
    delProb, dupProb, invProb, transProb, fbackProb = sv_probs["del"], sv_probs["dup"], sv_probs["inv"], sv_probs["trans"], sv_probs["fback"]
    safeIds = set([x.seg_id for x in newSegL if x.preserve])

    i = 0
    while i < int(math.ceil(num_breakpoints/2.0)):
        print(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in newSegL]) + "\n")
        currUniqueSegs = len(set(newSegL))
        currLength = float(sum(s.size for s in newSegL))
        lenDiff = (currLength - origLength)
        print(lenDiff)
        meanSegsTogether = int(math.ceil(math.log2(len(newSegL))))

        # select runif strtpoint
        zeroLen = True
        manipSegs, unManipSegs, manipLen, numSegs, loopedEndI = [], [], 0, 0, 0
        while zeroLen:
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
                print(len(newSegL), len(manipSegs), len(unManipSegs))

        currSafeIds = set([x.seg_id for x in unManipSegs if x.preserve])
        delAllowed = (currSafeIds == safeIds)

        # deletion
        # print(delAllowed,delProb,(manipLen - lenDiff) / origLength, currUniqueSegs / origUniqueSegs, float(len(manipSegs)) / len(newSegL))
        if delAllowed and r.random() < delProb and (manipLen - lenDiff) / origLength < 0.2 and \
                currUniqueSegs / origUniqueSegs > 0.4 and float(len(manipSegs)) / len(newSegL) <= 0.4:
            newSegL = unManipSegs
            i += 1
            print("del")

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
                    print(manipLen,lenDiff,origLength,(currLength + manipLen) / origLength)
                    dup = True

                # is inverted
                if r.random() < invProb:
                    inv = True

            print(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in manipSegs]) + " | " + str(manipLen))
            print("trans,dup,inv: " + str((trans, dup, inv)))
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

    print(",".join([str(x.seg_id) + "+" if x.direction == 1 else str(x.seg_id) + "-" for x in newSegL]) + "\n")
    logging.info("Final size: " + str(sum(s.size for s in newSegL)))
    return newSegL


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
