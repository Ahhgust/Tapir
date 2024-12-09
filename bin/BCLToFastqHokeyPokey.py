#!/usr/bin/env python3
'''Mangle a sample sheet.'''

import sys

hexBytes = "0123456789ABCDEF"

workF = sys.stdin
if len(sys.argv) == 2:
    workF = open(sys.argv[1], "r")
elif len(sys.argv) > 2:
    raise ValueError("Too many files provided...")

def hexifyThing(toHex):
    if "_" in toHex:
        raise ValueError("Underscore in value: " + toHex)
    hexSampN = bytearray(toHex, "utf-8")
    hexSampN = "".join([hexBytes[15 & (cv >> 4)] + hexBytes[15 & cv] for cv in hexSampN])
    return hexSampN


preData = True
dataHead = True
noneSeen = True
dataIndex = -1
nameIndex = -1
for line in workF:
    if line.strip() == "":
        continue
    lineS = [cv.strip() for cv in line.split(",")]
    if preData:
        if lineS[0] == "[Data]":
            preData = False
    elif dataHead:
        if "Sample_ID" in lineS:
            dataIndex = lineS.index("Sample_ID")
        if "Sample_Name" in lineS:
            nameIndex = lineS.index("Sample_Name")
        if (dataIndex < 0) and (nameIndex < 0):
            raise IOError("Both Sample_ID and Sample_Name are not present.")
        dataHead = False
    else:
        if dataIndex >= 0:
            lineS[dataIndex] = hexifyThing(lineS[dataIndex])
        if nameIndex >= 0:
            lineS[nameIndex] = hexifyThing(lineS[nameIndex])
        noneSeen = False
    print(",".join(lineS))
if noneSeen:
    raise ValueError("No sample IDs present in sample sheet.")


