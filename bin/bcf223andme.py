#!/usr/bin/env python3
# Author: August Woerner
# This is a little script to convert a VCF file (stdin)
# into a "23andme" file. Optionally gzipped.
# for VCFs that have been lifted over to hg19, it is highly recommended that you
# remove records from inconsistent chromosomes, a la:
# bcftools view -e 'INFO/OriginalContig!=CHROM' ...
# (and that you have GATK add the OrginalContig tag)
# One, it makes sense to do, and two, you can get (haploid) X calls
# showing up on the autosomes, which also breaks many bcftools routines.

import sys
import argparse
import gzip

EPSILON=0.0001


def getInfo(infoStr, query, cast=str):
    """
    Takes an INFO string and a query; returns the value of the string.
    eg, AF=0;INFO=1;OriginalContig=chr1;OriginalStart=14487;RAF=0.00147929
    query=RAF will return [0.00147929]
    (set cast=float and it'll be a number)
    """
    
    s = infoStr.split(";")
    for item in s:
        
        sp = item.split("=")
        if len(sp)<2:
            continue
        if sp[0] == query:
            return [cast(f) for f in sp[1].split(",")]
    
    return None
    
    
def getTagPosition(tagString, query):
    """
    Takes: GT:PL:DP (tag string)
    and a query (PL)
    and returns the 0-based index of the tag string (1 in the example)
    or -1 if not found.
    """
    tags = tagString.split(":")
    for i in range(len(tags)):
        if tags[i]==query:
            return i
    return -1


def probToOdds(val, eps=EPSILON):
    """
    Converts a probability into an odds
    values at 0 and 1 are adjusted by epsilon to make the odds non-infinite and non-zero
    """
    if val >= 1.-eps:
        val = 1.-eps
    elif val < eps:
        val = eps

    return val/(1.0-val)


parser = argparse.ArgumentParser(
    prog=sys.argv[0],
    description='This takes a (single sample) VCF (standard input ONLUY), and writes a "23 and me" genotype file. Best used when piped from bcftools',
    epilog='Use with:\nbcftools view -e \'INFO/OriginalContig!=CHROM\' ...\nfor VCFs that have been lifted over using GATK/Picard')

AF_TAG='RAF' # for GLIMPSE2 only; used with the Bayes factor 
parser.add_argument('-a', '--autosomes', dest='A', action='store_true', help="Just the autosomes, please!")
parser.add_argument('-x', '--x_and_autosomes', dest='X', action='store_true', help="Just the autosomes and the X, please!")
parser.add_argument('-s', '--sample_index', dest='S', type=int, default=1, help="The index (1-based) of the individual in the vcf file. Defaults to the 1st person.")
parser.add_argument('-z', '--gzip', dest='Z', action='store_true', help="Gzip the output")
parser.add_argument('-l', '--liftover_sanity', dest='L', action='store_true', help="Sanity test for a liftover VCF; ensures that the chromosomes (before/after) match", default=False)
parser.add_argument('-o', '--out', dest='O', action='store_true', help="Writes to standard out...")

parser.add_argument('-q', '--min_quality', dest='Q', type=int, help="The minimum genotype quality (second-lowest PL tag value); use with BCFtools", default=-1)
parser.add_argument('-p', '--min_posterior', dest='P', type=float, help="The minimum genotype posterior probability (taken as the max GQ value); use w/ GLIMPSE2", default=-1.)
parser.add_argument('-b', '--min_bayes_factor', dest='B', type=float, help="The minimum Bayes factor (requires allele frequency estimates; taken as max posterior odds / corresponding genotype's prior odds); use w/ GLIMPSE2", default=-1.)
parser.add_argument('-k', '--keep_indels', dest='K', action='store_true', help="Reports indels, as well as SNPs.")
parser.add_argument('-B', '--bayes_factor_af_tag', dest='AF_TAG', type=str, help="Use w/ imputation solvers other than GLIMPSE2; what tag to use for the alternative allele frequency; defaults to RAF (glimpse default)", default=AF_TAG) 



flags = parser.parse_args()
autosOnly=flags.A or flags.X
alsoX = flags.X
calculateBF=False
AF_TAG = flags.AF_TAG

whichIndex=flags.S + 8 
f = None

tag=""
if flags.Q>-1:
    tag='PL'
elif flags.P>-1:
    tag='GP'
    if flags.Q>-1:
        print("-p and -q are set! Pick one or the other, folks.", file=sys.stderr)
        exit(1)

if flags.B > 0:
    calculateBF=True
    # we filter on a Bayes factor AND a genotype posterior.
    if tag == "": # trivial; must be > 0, if not set, making a trivial AND
        tag='GP'
        flags.P=0. # add a trivial value for GP
        
    elif tag != 'GP':
        print(tag, " is not allowed. Filtering on a Bayes factor only makes sense for GLIMPSE", file=sys.stderr)
        exit(1)
        
ntot=0
nhet=0
nhalt=0
        
for line in sys.stdin:
    if line.startswith("##"):
        continue
    elif line.startswith("#"):
        sp = line.rstrip().split("\t")
        if whichIndex >= len(sp) or whichIndex < 9:
            print("That's a problem. The specified index " , (whichIndex+1) , " doesn't make a whole lot of sense!", file=sys.stderr)
            exit(1)
        who = sp[whichIndex]
        
        if flags.Z:
            f= gzip.open(who + ".23andme.tsv.gz", "wt")
        elif flags.O:
            f = sys.stdout
        else:
            f = open( who + ".23andme.tsv", "w")
            
        print("# rsid\tchromosome\tposition\tgenotype", file=f)
        
    else:
        sp = line.rstrip().split("\t")
        rec = sp[whichIndex]

        rs = sp[2]
        chrom = sp[0]
        if chrom.startswith("chr"):
            chrom = chrom[3:]

        # VCF assumed to be in genomic order (1-22, X,Y,M,...)
        if autosOnly and (chrom[-1] < '0' or chrom[-1]>'9'):
            if not alsoX:
                break
            elif chrom[-1] != "X":
                break
                
        # get the alleles at a locus
        alleles = [sp[3]]
        alleles.extend( sp[4].split(",") )
        gt = rec.split(":")[0]
        if gt.startswith("."):
            continue


        # biallelic SNPs only.
        if len(alleles)>2:
            continue
        if not flags.K: #NOT recommended to keep indels...
            if len(alleles[0])>1:
                continue
            if len(alleles)==2 and len(alleles[1])>1:
                continue       


        # sanity test.
        if flags.L:
            af=getInfo(sp[7], "OriginalContig", str)
            if af is not None:
                # chromosome annotations do not match
                # probably don't want to use these for IBD-segment inference...
                if sp[0] != af[0]:
                    continue

        # phase or unphased, we do not care
        gts = gt.split("/")
                
        if len(gts)==1:
            gts = gt.split("|")

        isHemi=False
        # "1" to 1
        gts = [int(x) for x in gts]
        # hemizygous calls get "diplotyzed"
        if len(gts)==1:
            gts.append(gts[0])
            isHemi=True
        
        finalgts = [ alleles[x] for x in gts ]

        tagLabels=sp[8]
        tagIndex=-1
        
        bf=prior=None
        if calculateBF:
            af=getInfo(sp[7], AF_TAG, float)
            if af is None:
                print("Should never happen", sp[7], file=sys.stderr)
                continue
            if len(af)>1:
                continue
            
            #p^2 
            prior=(1.0-af[0]) * (1.0-af[0])
            if finalgts[0] != finalgts[1]: # 2pq
                prior= (1.0-af[0])*(af[0])*2
            # from here, we only observe one allele
            elif not isHemi:
                if gts[0] != 0: # and it's the alternative allele
                    prior = af[0]*af[0]
            elif gts[0]==0: # hemizygous reference
                prior = 1.0-af[0]
            else: # hemizygous alternative
                prior=af[0]

            prior=probToOdds(prior)
            
            
        if tag != "":
            tagIndex=getTagPosition(sp[8], tag)

            if tagIndex < 0:
                continue
            tagVals = rec.split(":")[tagIndex].split(",")
            
            # used w/ bcftools
            if tag == "PL":
                tagVals = [int(i) for i in tagVals if not i.startswith(".")] # need integers
                tagVals.sort()
                if len(tagVals)< 2:
                    continue
                if tagVals[1] < flags.Q: # second-smallest likelihood is taken as the genotype quality. (this is the legacy GATK definition of genotype quality)
                    continue
            elif tag == 'GP': # used w/ GLIMPSE
                tagVals = [float(i) for i in tagVals] # need reals
                if len(tagVals) < 1:
                    continue
                elif  max(tagVals) < flags.P: # posterior probability is too small.
                    continue
                    
                
                if calculateBF:
                    maxProb=max(tagVals)
                    maxProb = probToOdds(maxProb)
                    
                    if prior is not None and prior > 0:
                        bf = maxProb/prior                        
                        if bf < flags.B:
                            continue
                    else: # prior is 0 (can't happen) or none (shouldn't happen); let's just skip those...
                        print("Unexpected prior", line, file=sys.stderr)
                        continue


        # note this will not play nice w/ how the genotypes are printed (all letters concatenated)...
        if not flags.K:
                # skipping insertions.
            if max( [ len(x) for x in finalgts ] ) > 1:
                continue

            # and deletions
            if min( [ len(x) for x in finalgts ] ) < 1:
                continue


        ntot+=1
        if finalgts[0] != finalgts[1]:
            nhet+=1
        elif finalgts[0] != alleles[0]: # homozygous (else) and not equal to the reference
            nhalt+=1

        print(rs, chrom, sp[1], "".join(finalgts), sep="\t", file=f)
        
if f is not None and f is not sys.stdout:
    f.close()

# a little janky, but let's print some summaries to stderr
# note auto/X call
print("Ntotal:", ntot, "NHet:",  nhet, "NHomAlt:", nhalt, sep="\t", file=sys.stderr)