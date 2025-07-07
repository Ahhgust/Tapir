#!/usr/bin/env python3

import os
import sys
import argparse

VERSION=0.001
EPSILON=0.0001

class AllPairs:
    """
    This is a poor mans itertools cross
    It lets you iterate through allpairs of values in two lists
    (ivals, jvals)
    """
    def __init__(self, ivals, jvals):

        self.imax=len(ivals)
        self.jmax=len(jvals)
        
        self.ivals=ivals
        self.jvals=jvals
        
    def __iter__(self):
        self.i=0 # indexes
        self.j=0
        return self

        
    def __next__(self):

        if self.j>= self.jmax:
            self.i=self.j=0        # note that the iterator resets itself for the next use
            raise StopIteration

        # the current iteration value
        x=self.ivals[self.i], self.jvals[self.j], self.i, self.j

        # and set up the next indices
        self.i += 1
        if self.i < self.imax:
            return x

        self.i= 0
        self.j+= 1        

        return x

def getInfo(tag, dat, callback=max, cast=float):
    dat = dat.split(";")
    for d in dat:
        sp = d.split("=")
        if len(sp)==2 and sp[0]==tag:
            v = [cast(f) for f in sp[1].split(",") if f != '.']
            return callback(v)
        
    return None

# used to parse the GP tag.
# returns the max
def getFormat(tag, names, vals, callback=max, cast=float):
    names=names.split(":")
    vals = vals.split(":")
    if len(names) != len(vals):
        print(names, vals, file=sys.stderr)
        return None
        
    for i in range(len(names)):
        if names[i] == tag:
            v = [cast(f) for f in vals[i].split(",") if f != '.']
            if callback is None:
                return v
            return callback(v)
        
    return None

def prob2odds(prob, eps=EPSILON):
    
    # avoids division by 0
    if prob < eps:
        prob = eps
    elif prob>=(1.0-eps):
        prob=1.0-eps

    return prob/(1.0-prob)
        
if __name__ == "__main__":


    # Room to grow
    # We can onsider multiple genotypers. For now, glimpse is it.
    formatFilter='GP'
    fLevels = [0.50, 0.75, 0.90, 0.95, 0.99]

    infoFilter='RAF'
    iLevels= [1,5,10,50,100,200, 500, 1000] # note that the filter is not on the RAF, but on the Bayes Factor (which is a function of the RAF)


    allPairs = AllPairs(fLevels, iLevels)

    
    sampleid=""
    for line in sys.stdin:
        if line.startswith("##"):
            print(line, end="")
        elif line.startswith("#"):
            print("##", formatFilter, "levels=", ",".join([str(f) for f in fLevels]), sep="")
            print("##", infoFilter, "levels=", ",".join([str(i) for i in iLevels]), sep="")
            
            
            s = line.rstrip().split("\t")
            if len(s) != 10:
                print("Unexpected input. I need a a single-sample VCF to work...", file=sys.stderr)
                exit(1)
                
            sampleid=s[-1]
            
            print("\t".join(s), end="")

            for p in allPairs:
                print("\t", p[-2] , "-" , p[-1], sep="", end="")
            print()
            
        else:
            s = line.rstrip().split("\t")
            gtall=s[-1].split(":")
            
            gt=gtall[0]

            #no missing data
            if gt[0]=='.':
                continue

            # can either be diploid or haploid; 
            if len(gt) != 1 and len(gt) != 3:
                continue

            # and strictly biallelic
            if gt[0] != '0' and gt[0] != '1':
                continue

            if gt[-1] != '0' and gt[-1] != '1':
                continue


            raf = getInfo(infoFilter, s[7]) # alternative allele frequency

            if raf is None:
                continue
            
            gps = getFormat(formatFilter, s[8], s[9], None)
            if gps is None:
                continue
            
            gp = max(gps)
            if len(gps)>3:
                continue
            
            postodds= prob2odds(gp)

            # hemizygous
            if len(gt)==1:
                if gt=="0":
                    priorodds=prob2odds( 1.0-raf)
                else:
                    priorodds=prob2odds(raf)
            elif gt[0] != gt[-1]: #0/1
                priorodds=prob2odds( 2*raf*(1.0-raf) ) #2pq
            elif gt[0] == '0':
                priorodds=prob2odds( (1.0-raf)*(1.0-raf)) #p2
            else:
                priorodds=prob2odds( raf*raf) # q2
            

            bf = postodds/priorodds
                
            print("\t".join(s), end="")

            
            for p in allPairs:
                gpFilt=p[0]
                bfFilt=p[1]

                if gp > gpFilt and bf > bfFilt:
                    print("\t", s[-1], sep="", end="")
                else:
                    print("\t./.", end="")

            print("")
            
            
    
        
    


