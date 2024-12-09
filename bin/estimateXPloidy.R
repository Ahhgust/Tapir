#!/usr/local/bin/Rscript
# Written by August Woerner
# This script estimates the ploidy of the X chromosome relative to the autosomes (either 2:1, AA:X(possibly Y) or 2:2, AA:XX (also possibly XXY))
# It does so using a likelihood ratio test using the total number of reads aligned to the autosomes and the X
# See usage statement
# In principle, you can do the same thing for the Y, but you'd need to define how much of the "Y" you assayed in females...
# this script neglects mapping error;

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==2 || length(args)==4) {

    # default probabiliities, assuming that the regions on the X and the autosomes are equal in size
    xxprob=0.5
    xprob=2./3

    if (length(args)==4) {
        atot <- as.numeric(args[[3]])
        xtot <- as.numeric(args[[4]])

        xxprob <- atot/(atot+xtot) # expected number of successes if ploidies are equal, and autosomal reads are "success"
        xprob <-  (atot*2)/((atot*2)+xtot) # ditto for 2:1 ploidy
    }

    acount <- as.integer(args[[1]])
    xcount <- as.integer(args[[2]])
    

    xxlike <- dbinom(acount, acount+xcount, prob=xxprob, log=TRUE)
    xlike <- dbinom(acount, acount+xcount, prob=xprob, log=TRUE)

    cat(
        sprintf("%f\t%f\t%f\t%d\n", xxlike, xlike, xxlike-xlike, ifelse(xxlike>xlike, 2, 1))
    )

} else {
    stop("Usage: guessXPloidy.R aCount xCount (aTot xTot)\nBy default, the size of the autosomal and x chromosome regions are assumed to be equal; otherwise provide the sizes of the two regions (aTot, xTot)")
}





