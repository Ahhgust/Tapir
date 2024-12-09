#!/usr/bin/env python3

import subprocess
import argparse
import sys
import os

weights = [None] * 11

# these are filtering rules compatible with bcftools
# note that the "GQ" (from the paper) is the second-highest (also second lowest) PL score.
# The PL should look like: 5,0,100
# which consider an order of genotypes (AA,AB,BB, where A is the reference, and B is the alternative allele)
# the most likely genotype should have a PL of 0 (first sanity check)
#
# note that the sum of the above PLs (105) - the max(100) gives you 5 (the value you want).
# this expression only works when there are three genotype possibilities. the second boolean caters to this possibility (more than three GTs; assuming a single PL of 0, the sum must be the max + the second-from-max)

weights[1] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35)) || (GT="hom"&(QUAL<62))'
weights[2] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35 | SUM(FMT/PL)-MAX(FMT/PL)<55)) || (GT="hom"&(QUAL<100))'
weights[3] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35 | SUM(FMT/PL)-MAX(FMT/PL)<55)) || (GT="hom"&(QUAL<131))'
weights[4] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35 | SUM(FMT/PL)-MAX(FMT/PL)<55)) || (GT="hom"&(QUAL<155))'
weights[5] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35 | SUM(FMT/PL)-MAX(FMT/PL)<55)) || (GT="hom"&(QUAL<156))'
weights[6] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35 | SUM(FMT/PL)-MAX(FMT/PL)<55 | MIN(FMT/AD)/MAX(FMT/AD) < 0.51 & SUM(FMT/PL)-MAX(FMT/PL)<77)) || (GT="hom"&(QUAL<178))'
weights[7] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35 | SUM(FMT/PL)-MAX(FMT/PL)<57 | MIN(FMT/AD)/MAX(FMT/AD) < 0.51 & SUM(FMT/PL)-MAX(FMT/PL)<77)) || (GT="hom"&(QUAL<178 | SUM(FMT/PL)-MAX(FMT/PL)<24))'
weights[8] = None # the decision tree learned a nn-monotonic relationship. let's just ignore this class
weights[9] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.41)) || (GT="hom"&(QUAL<179 | SUM(FMT/PL)-MAX(FMT/PL)<24))'
weights[10] = 'MIN(FMT/PL) > 0 || STRLEN(ALT[1]) > 0|| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.41)) || (GT="hom"&(QUAL<199 | SUM(FMT/PL)-MAX(FMT/PL)<27))'



parser = argparse.ArgumentParser(description="Let's +setGT that b/vcf file!")

parser.add_argument('-d', '--dry_run', action='store_true', dest='D', help="Just prints the bcftools command to standard out")
parser.add_argument('-v', '--vcf_file', dest='V', default="-", help="The name of the b/vcf file. - (default) is standard input")
parser.add_argument('-w', '--weight', dest='W', default=10, type=int, help="The weight (stringency) of filtering. 1-10; no 8")
parser.add_argument('-o', '--output_file', dest='O', default="", help="Writes data to the named file. Default is standard out")
parser.add_argument('-i', '--index_file', dest='I', action='store_true', help="Automagically indexes the result")
parser.add_argument('-p', '--print_rules', dest='P', action='store_true', help="Prints out the filtering rules and their index (weight, -w value)")


(results, args) = parser.parse_known_args(sys.argv[1:])[0:2]


if results.W < 1 or results.W >= len(weights) or weights[ results.W ] is None:
    print("Index out of bounds; you asked for weight:" , results.W, " but no such weight exists!", file=sys.stderr)
    exit(1)

if results.P:
    print("Weight" , " : " , "Rules")
    for i in range(len(weights)):
        print(i, " : " , weights[i])
    exit(0)
    
# we're emulating the following:
#bcftools +setGT -Ov Mix1to4-1ng.la.md.13X.bcf -- -t q -n . -i'MIN(PL) > 0| (GT="het"&(MIN(FMT/AD)/MAX(FMT/AD) < 0.35)) | (GT="hom"&(QUAL<62))'

if results.O:
    outputFormat='v'
    if len(results.O) >3 and results.O[-3:] == 'bcf': # a bit crude, but if you ask for a bcf, let's have it make one!
        outputFormat='b'

    # the if/else setup is a bit clunky here; the gist of the problem is that ' is not interpreted when we use subprocess, but it is when we use os.system.
    # and this is a rather janky way to deal with it. Setting shell=True might also work...?
    if results.V == '-' or results.V == '/dev/stdin/':
        command = "bcftools +setGT -O%s -o %s %s -- -t q -n . -i%s" % (outputFormat, results.O, results.V, weights[ results.W ].replace(" ", "") )
    else:
        command = "bcftools +setGT -O%s -o %s %s -- -t q -n . -i'%s'" % (outputFormat, results.O, results.V, weights[ results.W ].replace(" ", "") )
        
    if results.I and outputFormat=='b':
        command += " && bcftools index " + results.O 
    
else:
    if results.V == '-' or results.V == '/dev/stdin/':
        command = "bcftools +setGT -Ov %s -- -t q -n . -i%s" % ( results.V, weights[ results.W ].replace(" ", "") )
    else:
        command = "bcftools +setGT -Ov %s -- -t q -n . -i'%s'" % ( results.V, weights[ results.W ].replace(" ", "") )

#command = "echo " + command # remove me later...

if results.D:
    print(command)

elif results.V == '-' or results.V == '/dev/stdin/':
    #    print("Not implemented yet", file=sys.stderr)
    #subprocess.run( ["cat"], input=sys.stdin.read().encode())
    subprocess.run(  command.split(), input=sys.stdin.read().encode())

else:
    os.system(command)


