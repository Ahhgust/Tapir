#!/usr/bin/env python3

import pysam
import sys
import os
import argparse

MAX_DEPTH=512
MAX_READLEN=351
MAX_TLEN=1023




def countAlignedBases(cigar):
    '''
    This takes in a cigar string and counts the number of aligned bases
    Taken as the number of Ms and Is
    '''

    nchar = len(cigar)
    num=0
    mtag=0

    for i in range(nchar):
        if cigar[i] >= '0' and cigar[i] <= '9':
            num *= 10
            num += int(cigar[i])
        elif cigar[i] == 'M' or cigar[i]=='I':
            mtag += num
            num=0
        else:
            num=0

    return mtag



def main(args):

    parser = argparse.ArgumentParser(description="Let's summarize bams!")
    parser.add_argument('-v', '--verbose', dest='V', help="Verbose option; prints raw stats, not marginals", action='store_true')
    parser.add_argument('-r', '--reference', dest='R', help="The reference sequence (fasta). Necessary for CRAM files", default="")
    parser.add_argument('-V', '--very_verbose', dest='VV', help="Very Verbose option; prints raw stats, not marginals; includes duplicates", action='store_true')
    parser.add_argument('-m', '--mapping_quality', dest='M', help="Only retains reads with a minimum mapping quality of Q in the final output", type=int, default=20)

    results = parser.parse_known_args(args[1:])[0]
    argv = parser.parse_known_args(args)[1]


    if results.VV:
        results.V=True
    
    if len(argv) < 2:
        print("I need a bam file and an optional bed file, yo!")
        return 1
    
    bam = argv[1]

    if len(argv) > 2:
        bed = argv[2]
        if bed.upper() == "GSA":
            bed = os.path.join( os.path.dirname(os.path.realpath(__file__)), "hg38.gsa_autos_filtered.bed") # shortcut for the GSA
    else:
        bed = os.path.join( os.path.dirname(os.path.realpath(__file__)), "hg38.10ksamp.bed") # default bed file-- keep this in my src/
    

    regions = list()
    with open(bed, "r") as fh:
        for line in fh:
            sp = line.rstrip().split("\t")
            if len(sp) < 3:
                if len(sp):
                    print("Unexpected line: ", line, file=sys.stderr)

                continue
            
            regions.append( (sp[0], int(sp[1]), int(sp[2])))

    if bam.endswith(".cram"):
        if results.R == "":
            print("CRAM detected; I need a reference filename to work!", file=sys.stderr)
            exit(1)
            
        bammy = pysam.AlignmentFile(bam, "rc", reference_filename=results.R)
    else:
        bammy = pysam.AlignmentFile(bam, "rb")

    pileDepths = [0] * (MAX_DEPTH+1)
    pileDepthsWithDups = [0] * (MAX_DEPTH+1)
    readLens = [0] * (MAX_READLEN+1)
    templateLengths = [0] * (MAX_TLEN+1)
    
    # in the sample, the maxima
    maxReadLen=0
    maxDepth=0
    maxDepthWithDups=0
    maxTLen=0
    if results.V:
        print("Chrom\tStop\tRefStart\tReadLen\tTlen\tIsDuplicate")
    for locus in regions:
        it = bammy.fetch( locus[0], locus[1], locus[2])
        
        nreads=0
        nreadsWithDups=0
        readnames = set()
        for c in it:


            if c.mapping_quality<results.M:
                continue
            
            if c.cigarstring is None:
                continue
            
            if c.query_name in readnames:
                continue

            readnames.add(c.query_name) #ignore overlapping reads

            nreadsWithDups +=1
            
            if not c.is_proper_pair or c.query_name is None:
                continue


            if not c.is_duplicate:
                nreads += 1
               
                nbases = countAlignedBases(c.cigarstring)
                
                if nbases > MAX_READLEN:
                    maxReadLen= nbases = MAX_READLEN
                elif nbases > maxReadLen:
                    maxReadLen=nbases

                readLens[nbases]+=1

                tlen = abs( c.template_length )
                                
                if tlen > MAX_TLEN:
                    maxTLen= tlen = MAX_TLEN
                elif tlen > maxTLen:
                    maxTLen=tlen

                templateLengths[tlen]+=1

                
            if results.V:
                if c.is_duplicate:
                    if not results.VV:
                        continue
                    
                    tlen = abs( c.template_length )
                    nbases = countAlignedBases(c.cigarstring)
                
                print(locus[0], locus[2], c.reference_start, nbases, tlen, c.is_duplicate, sep="\t")

        
        if nreads > MAX_DEPTH:
            nreads = maxDepth = MAX_DEPTH
        elif nreads > maxDepth:
            maxDepth=nreads

        pileDepths[ nreads ] += 1
        
        if nreadsWithDups > MAX_DEPTH:
            nreadsWithDups = maxDepthWithDups = MAX_DEPTH
        elif nreadsWithDups > maxDepthWithDups:
            maxDepthWithDups=nreadsWithDups

        pileDepthsWithDups[ nreadsWithDups ] += 1


    if results.V:
        return 1
    
    print("Index\tCount\tLabel")
    print("0", len(regions), "Nsites", sep="\t")
    
    for i in range(maxReadLen+1):
        print(i , readLens[i], "ReadLength", sep="\t")

    for i in range(maxDepth+1):
        print(i , pileDepths[i], "Depth", sep="\t")

    for i in range(maxDepthWithDups+1):
        print(i , pileDepthsWithDups[i], "DepthWithDups", sep="\t")

    for i in range(maxTLen+1):
        print(i , templateLengths[i], "TemplateLengths", sep="\t")

    return 1
    
if __name__ == "__main__":
    main(sys.argv)

    
    
    
    


