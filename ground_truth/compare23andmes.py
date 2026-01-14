#!/usr/bin/python3
#
# Written by August Woerner
# Compares genotypes; takes in a ground truth file (-t), and a comparison file (-c)
# both need to be "23andme" format.
# Optionally, the comparison can be constrained to some SNP panel (vcf or vcf.gz)
# 

import sys
import os
import gzip
import argparse
from collections import namedtuple

Genotype = namedtuple("Genotype", "a b")

REF_AL=3
ALT_AL=4
INFO_COL=8

header=[]
indels=0



parser = argparse.ArgumentParser(description="Let's compare some VCFs!")
parser.add_argument('-s', '--snp_panel', dest='S', help="A snp panel file; reduce comparison to SNPs in this set; vcf or vcf.gz", default="")
parser.add_argument('-t', '--truth', dest='T', help="Ground truth; 23andme format (gzipped okay)", default="")
parser.add_argument('-c', '--comparison', dest='C', help="Comparison (to ground truth); 23andme format (gzipped okay)", default="")

results, args = parser.parse_known_args(sys.argv[1:])[0:2]


if len(args):
    print("Extra arguments detected. Not good!\n", file=sys.stderr)
    exit(1)

if results.T=="":
    print("No truth data?!\nThat's not good!", file=sys.stderr)
    exit(1)

if results.C=="":
    print("No comparison data?!\nThat's not good!", file=sys.stderr)
    exit(1)
    
    
def getGT(alleles):
    """
    s is a character string; must be two letters (DNA alphabet)...
    """
    if len(alleles) != 2:
        print("Parsing failure with ", alleles, file=sys.stderr)
        exit(1)

    if alleles[0] < alleles[1]:
        return Genotype( alleles[0], alleles[1] )
    return Genotype( alleles[1], alleles[0] )

def parseLine(line):
    sp = line.rstrip().split("\t")
    gt = getGT(sp[3])
    return sp[1], sp[2], gt

def getFH(f):
    if f.endswith(".gz"):
        fh = gzip.open(f, "rt")
    else:
        fh=open(f)
    return fh


def parsePanel(f):
    
    fh = getFH(f)

    d = {}

    for line in fh:
        if line.startswith("#"):
            continue
        line = line.rstrip()
        s = line.split("\t")
        # chromosome + position + ref allele
        
        if s[0][0] == 'c': # stripping off the 'chr' prefix
            key = s[0][3:] + ":" + s[1]
        else:
            key = s[0] + ":" + s[1]
        
        # one or more alt alleles
        val = set(s[4].split(","))
        val.add(s[3])

        # panel can be norm (-m-)
        # ie, the same position may have multiple variants on multiple lines
        if key in d:
            d[key] = d[key].union( val )
        else:
            d[key]=val
        
    fh.close()
    return d


def getCat(gt1,gt2):
    '''
    Takes a pair of gt objects (from getGt)
    And returns a nominal category
    gt1 is taken from ground truth
    '''

    # handle the cases of missing data (in one or the other)
    if gt1.a is None:
        cat1="Missing"
        if gt2.a is None:
            return "Missing,Missing"
        if gt2.a==gt2.b:
            return "Missing,Hom"
        return "Missing,Het"


    type1='Het'
    if gt1.a==gt1.b:
        type1="Hom"
    
    if gt2.a is None:
        return type1 + ",Missing"

    type2='Het'
    if gt2.a==gt2.b:
        type2="Hom"

    # common case; match
    if gt1 == gt2:
        return type1 + "," + type1

    if type1=='Hom':
        if gt1.a == gt2.a or gt1.a== gt2.b: # 1 allele matches
            return "Hom,Het"
        elif type2 == 'Hom': # homozygous and inconsistent
            return "Hom,AltHom"
        return "Hom,AltHet" # and crazytown

    # type1==Het
    if type2=='Hom':
        if gt2.a == gt1.a or gt2.a==gt1.b:
            return "Het,Hom" # 1 allele matches, the other direction
        return "Het,AltHom"

    # types 1 and2 == Het
    if gt1.a in gt2 or gt1.b in gt2: # both het, one allele matches; eg, A/C vs A/G)
        return "Het,AltHet"
    return "Het,AltAltHet" # neither allele matches (eg, A/C vs G/T)
        


catcounts = dict.fromkeys(
        [\
            "Hom,Hom", "Het,Het","Hom,Het", "Het,Hom", \
            "Hom,AltHom", "Hom,AltHet", "Het,AltHet", "Het,AltAltHet", "Het,AltHom",\
            "Missing,Missing","Missing,Hom", "Missing,Het",\
            "Hom,Missing", "Het,Missing", "Indels"],\
        0 )
        

byCompartment = {}

panel={}
if results.S != "":
    panel = parsePanel(results.S)


truth = {}
fh = getFH(results.T)
for line in fh:
    # header
    if line.startswith("#"):
        continue
    
    (chrom, pos, gt2) = parseLine(line)
    
    key = chrom + ":" + pos
    if results.S != "":
        if not key in panel:
            continue

    if key not in truth:
        truth[key] = gt2
    else:
        truth[key] = Genotype(None, None)
        print("Duplicate truth records:", key, gt2, file=sys.stderr)

fh.close()

fh = getFH(results.C)
for line in fh:
    # header
    if line.startswith("#"):
        continue

    (chrom, pos, gt1) = parseLine(line)
    
    key = chrom + ":" + pos
    if key not in truth or truth[key].a is None:
        continue

    gt2 = truth[key]

    if results.S != "":
        if not key in panel:
            continue
        
        if gt1.a not in panel[key] or gt1.b not in panel[key]:
            continue

        if gt2.a not in panel[key] or gt2.b not in panel[key]:
            continue

    chromcat = 'Autos'
    if chrom == 'X':
        chromcat = 'chrX'

    if chromcat not in byCompartment:
        byCompartment[chromcat] = catcounts.copy()

    inner = byCompartment[chromcat]
    
    cat = getCat(gt1,gt2)
    inner[cat] += 1

        

keys = catcounts.keys()
print("DepthThreshold\tQualityThreshold\tPosteriorThreshold\tCompartment", "\t".join(keys), sep="\t")

for compartment, inner in byCompartment.items():
    print("--", "--", "--", compartment, "\t".join([ str(inner[k]) for k in keys]), sep="\t")


