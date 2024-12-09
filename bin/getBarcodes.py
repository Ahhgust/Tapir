#!/usr/bin/python3
# modified from: https://gencore.bio.nyu.edu/how-to-find-out-what-barcodes-are-in-your-undetermined-reads/
from operator import itemgetter
import sys, gzip

barcodes = {}

fh = sys.stdin

if len(sys.argv) > 1:
        f = sys.argv[1]
        if f.endswith(".gz"):
                fh = gzip.open(f, "rt")
        else:
                fh = open(f)

i = 0
for line in fh:
        if i % 4 == 0:
                
                bc = line.split(':')[-1].strip()
                if bc not in barcodes:
                        barcodes[bc] = 1
                else:
                        barcodes[bc]+=1

        i += 1

if fh != sys.stdin:
        fh.close()
        
total = sum(barcodes.values())
for k, v in sorted(barcodes.items(), key=itemgetter(1)):
        print(k, v, round(v/total*100, 2))

        
