#!/usr/bin/env python3
#
# This parses a "Stats.json" file generated from Illumina's bcl2fastq (version 2)
#
#

import sys
import json

if len(sys.argv) < 2:
    print("Give me a Stats.json file, please", file=sys.stderr)
    exit(1)

files = sys.argv[1:]

def getSum(dic, tag):
    """
    The tags we want are in some dictionary (which is in some nested list of dictionaries and/or lists)
    This extracts out the SUM of some tag
    """
    ret = 0
    if type(dic) == dict:
        if tag in dic:
            ret += int(dic[tag])
        for key, val in dic.items():
            ret += getSum(val, tag)
    elif type(dic) == list: # list of dictionaries
        for item in dic:
            ret += getSum(item, tag)
            
    return ret

print("RunNumber", "Lane", "TotalClustersPass", "PercentClustersPass", "YieldQ30", "PercentYieldQ30", sep="\t")

for fn in files:

    f = open(fn)
    
    dat = json.load(f)

    runNum = dat["RunNumber"]

    resultsByLane = dat["ConversionResults"]

    for laneResults in resultsByLane:
        laneNumber = laneResults["LaneNumber"]
        tc= laneResults["TotalClustersRaw"]
        tcPass= laneResults["TotalClustersPF"]
        demux = laneResults["DemuxResults"]

        yld = getSum(demux, "Yield")
        yldq30 = getSum(demux, "YieldQ30")

        print(runNum, laneNumber, tcPass, tcPass/tc, yldq30, yldq30/yld, sep="\t")

            



