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


tib <- lapply(args, getTib) %>% bind_rows()

nsites <- filter(tib, Label=='Nsites') %>% pull(Count) %>% unique()

if (length(nsites) != 1) {
 #   print(nsites)
#    stop("Variable number of sites detected.")
}


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


options(dplyr.print_max = 1e9)

depths %>%
    pivot_wider(names_from=Label, values_from=c("Mean", "Var")) %>%
    left_join(alsodepth, by="File") %>%
    mutate_if(is.double, ~sprintf(., fmt="%.5f")) %>%
    format_tsv() %>%
    cat()
#    arrange(Label, MeanDepth) %>%
#    print(n=nrow(depths))

