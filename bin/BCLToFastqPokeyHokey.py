#!/usr/bin/env python3
'''Unmangle sample names.'''

import os
import sys

tgtFold = sys.argv[1]

hexBytes = "0123456789ABCDEF"
possExts = [".fq",".fastq",".fastq.gzip","fastq.gz","fq.gzip","fq.gz"]

noneSeen = True
for sfn in os.listdir(tgtFold):
    origP = os.path.join(tgtFold, sfn)
    # only want files
    if not os.path.isfile(origP):
        continue
    # with an underscore
    underInd = sfn.find("_")
    if underInd < 0:
        continue
    # that are fastqs
    origExt = None
    for pe in possExts:
        if sfn.endswith(pe):
            origExt = pe
            break
    if origExt is None:
        continue
    # that start with only hexadecimal characters
    origHex = sfn[0:underInd]
    if len(origHex) % 2 != 0:
        continue
    postHex = bytearray()
    for i in range(0, len(origHex), 2):
        nib0 = hexBytes.find(origHex[i])
        nib1 = hexBytes.find(origHex[i+1])
        if (nib0 < 0) or (nib1 < 0):
            postHex = None
            break
        postHex.append((nib0 << 4) + nib1)
    if postHex is None:
        continue
    postHex = str(postHex, "utf-8")
    # rename
    postP = os.path.join(tgtFold, postHex + sfn[underInd:])
    os.rename(origP, postP)
    noneSeen = False

if noneSeen:
    raise ValueError("No output fastqs encountered")

