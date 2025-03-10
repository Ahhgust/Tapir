#!/usr/local/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Hmisc))

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("Gimme an Rscript file!\n")
}

getTib <- function(f) {
    tib <- read_tsv(f, col_types=cols(), progress=FALSE)
    tib$File <- f
    NSites <- filter(tib, Label=="Nsites") %>% pull(Count)
    tib$Nsites <- as.integer(NSites[[1]])
    return(tib)
}


# (a refactored) Lander-Waterman equation to estimate the library size
f <- function(x, c, n) {
    return( c/x - 1.0+exp(-n/x) )
}

# an R version of:
#https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java#L115
estimateLibraryFractionScalar <- function(nreads, nuniqueReads) {
     # Estimates the size of a library based on the number of paired end molecules observed
     # and the number of unique pairs observed.
     # <p>
     # Based on the Lander-Waterman equation that states:
     # C/X = 1 - exp( -N/X )
     # where
     # X = number of distinct molecules in library
     # N = number of read pairs
     # C = number of distinct fragments observed in read pairs

    ndups <- nreads - nuniqueReads

    
    if (nreads==0 || nuniqueReads==0) {
        return(-1)
    }
    if (ndups < 0) {
                                        #stop("Cannot happen")
        return(-2)
    }

    if (f(nuniqueReads, nuniqueReads, nreads) < 0) {
                                        #stop("Also cannot happen")
        return(-3)
    }

    m=1.0
    M=100.

    while (f(M*nuniqueReads, nuniqueReads, nreads)>0) {
        M = M*10
    }

    # bisect away!
    for(i in 1:40) {
        r=(m+M)/2.
        u=f(r*nuniqueReads, nuniqueReads, nreads)
        if (u==0) {
            break
        } else if (u > 0) {
            m=r
        } else {
            M=r
        }
    }
    LibrarySize=nuniqueReads* (m+M)/2. # estimate of the number of unique templates in the library
    return(    
        nuniqueReads/LibrarySize     # what fraction of those things have we sampled?
        )
        
        
}

# and let's make a vector-friendly version
estimateLibraryFraction <- Vectorize(estimateLibraryFractionScalar)


tib <- lapply(args, getTib) %>% bind_rows()

nsites <- filter(tib, Label=='Nsites') %>% pull(Count) %>% unique()



filter(tib,
       startsWith(Label, "Depth") | grepl("Length", Label) ) %>%
    group_by(File, Label) %>%
    dplyr::summarize(
        # MeanDepth= sum(Count*Index)/nsites[[1]],
               Mean=Hmisc::wtd.mean(Index, weights=Count),
               Var=Hmisc::wtd.var(Index, weights=Count),
#               LB0.025=Hmisc::wtd.quantile(Index, weights=Count, probs=0.025)[[1]],
 #              UB0.975=Hmisc::wtd.quantile(Index, weights=Count, probs=0.975)[[1]],
                                        #             NSites=nsites[[1]],
#               NSites=Nsites[[1]],
               .groups='keep'
              ) %>%
    ungroup() -> depths

                                        #options(scipen=999)
                                        #options(tibble.width=Inf)

#filter(tib, grepl(Label, "Length")) %>% # template and read lengths
 #   group_by(File, Label) %>%
  #  dplyr::summarize(

filter(tib,
       Label == "Depth" ) %>%
    group_by(File) %>%
    dplyr::summarize(
               Breadth1X=sum( Count[ Index>=1])/sum(Count),
               Breadth5X=sum( Count[ Index>=5])/sum(Count),
               Breadth10X=sum( Count[ Index>=10])/sum(Count),
               Breadth20X=sum( Count[ Index>=20])/sum(Count)
           ) -> alsodepth




filter(tib,
       Label == "Depth" | Label=='DepthWithDups' ) %>%
    group_by(File) %>%
    dplyr::summarize(
               FracSampledOfLibrary=
               estimateLibraryFraction(
                   sum(Count[Label=='DepthWithDups']*Index[Label=='DepthWithDups' ]),
                   sum(Count[Label=='Depth']*Index[Label=='Depth'])
                   )
           ) -> complexity

options(dplyr.print_max = 1e9)

depths %>%
    pivot_wider(names_from=Label, values_from=c("Mean", "Var")) %>%
    left_join(alsodepth, by="File") %>%
    left_join(complexity, by="File") %>%
    mutate_if(is.double, ~sprintf(., fmt="%.5f")) %>%
    format_tsv() %>%
    cat()
#    arrange(Label, MeanDepth) %>%
#    print(n=nrow(depths))

