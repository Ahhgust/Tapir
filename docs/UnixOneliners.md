# One liners!
Tapir is a unix utility. Here are some helpful "one-liners" (full disclosure; may not be just one line) to help pull useful bits out of Tapir.
All one liners print to screen; you can copy and paste data from the screen to Excel. <br>
When you do go to (in Excel), and select on column A (ie, where you pasted the data) <br>
Data -> Text to Columns -> Delimited -> (select Tab and Space, as well as "Treat consecutive delimiters as one")
Doing so will make you a nice tabular excel worksheet.

## Experiment-level stats

All instructions assume you are in the `Experiment` directory. I.e.,:
whatever experiment name you gave Tapir <br>

```
ls -1
```
should give you:
```
Run_Fastqs
Sample_Data
```
If the above directories are not reported by `ls`, please `cd` into the right directory.

### Do my UDIs have too much data?
```
$TAPIR/bin/fcat.pl -h Sample_Data/Offtargets/*/Reports/*cov | cut -f1,3
```
Grabs the filename and mean read depth (coverage)
<br>
If you see "large" values (for example, >0.01x; use your validation as a guide), that may indicate carry over. <br>

<br>
We sometimes also want to look at our reads <br>

```
fgrep properly Sample_Data/Offtargets/*/Reports/*flagstat | sed 's/+ 0 properly paired //' | tr ':()'' ' '
```
which reports the number (and %) of properly paired reads.


### How much data do I have?
Ignoring the "Offtarget" samples, you can get coverage a la:
```
$TAPIR/bin/fcat.pl -h Sample_Data/*/*/Reports/*cov | cut -f1,3- | fgrep -v Offtargets
```
As well as the %/number of properly paired reads as:

```
fgrep properly Sample_Data/*/*/Reports/*flagstat | sed 's/+ 0 properly paired //' | tr ':()' ' ' | fgrep -v Offtargets
```

### Are my samples mixtures?

```
$TAPIR/bin/fcat.pl -h Sample_Data/*/*/Reports/*demix.summary
```

See [QC.md](QC.md) for some ideas on how to interpret the output of *Demixtify*.



### What are my call rates?
Calculating the "call-rate" has some subjectivity. This is a genomic assay after all (so the divisor is... 3 billion? Or do you only care about SNPs? Or SNVs?). The raw number of genotype calls can be had a la:
For glimpse:
```
for file in Sample_Data/*/Uploads/glimpse2/*snps.tsv.gz; do echo -n $file ' '; zgrep -c -v '#' $file; done
```

Note that `*` is unix wildcard; it matches everything and anything. So in the above, the `/*/` will match every sample in the `Sample_Data/` directory.

For bcftools:
```
for file in Sample_Data/*/Uploads/bcftools/*snps.tsv.gz; do echo -n $file ' '; zgrep -c -v '#'  $file; done
```
This gives you the number of *snps* (i.e., SNVs with MAF>1%) typed by Tapir.

-  The total number of SNPs (Autos + X) 
   -  11115284
   -  from: `bcftools view -H $TAPIR/resources/snppanel_autos_x.vcf.gz | wc -l`
-  The total number of X SNPs:
   -  427912
   -  from: `bcftools view -H $TAPIR/resources/snppanel_autos_x.vcf.gz chrX | wc -l`
   

### How do I report genotype accuracy?
We provide "23andme" files for HG002, HG003 and HG004; these come from Genome in a Bottle [External Link](https://www.nist.gov/programs-projects/genome-bottle).
You can compare "23andme" files a la:
```
python3 $TAPIR/ground_truth/compare23andmes.py -t $TAPIR/ground_truth/HG004.23andme.tsv.gz -c $TAPIR/ground_truth/HG003.23andme.tsv.gz
```
(note that the above command is probably not what you want (sometimes teachers do this kind of thing, in this case, to force you not to just copy and paste your way through this document). If you haven't spotted the issue, the above compares HG003 and HG004's data. It maybe provides you an example of a 'null' accuracy'...?)
<br>
where `-t` (HG004 in the above) provides the "true" genotypes, and `-c` provides the comparison. The script provides the following:
-  DepthThreshold,QualityThreshold,PosteriorThreshold
   -  Ignore; we have other tools that fill in these values, they are meaningless here.
-  Compartment
   -  At present either Autos or chrX (for/if/when we have truth data on the X)
-  Hom,Hom
   -  In truth homozygous (first Hom), and estimated to be Homozygous and identical.
-  Het,Het
   -  In truth heterozygous (first Het) and estimated to be Heterozygous and identical
      -  The sum of `Hom,Hom` and `Het,Het` is the number of correct genotypes.
-  Hom,*
   -  For the remaining columns, the different ways a genotyping error can occur, given that the true genotype is homozygous
-  Het,*
   -  As above, but given that the true genotype is Heterozygous.
-  The number of genotype calls (intersected with the calls made in the ground truth) is the row-sum
   -  As an artifact of how the computation is performed on 23andme files, no missing data (Missing) are available.    
   -  For convenience, you can estimate call rates from the row-sum as well
      -  It will be slightly biased (too small). We say a genotype is correct if it matches ground truth. If ground truth has missing data, that implicitly becomes a "no call" (because we can't say if the call is right or not). A consequence is that perfectly completed genotyping will result in call-rates that are not 100%, but closer to 99.9%)
<br>

**We highly recommend using the SNP 23andme files (.snps.tsv.gz), otherwise your "accuracy" will be inflated (glimpse).**
To understand why, recall that the vast majority of polymorphisms in populations are rare (remember that new alleles arise in an individual, and the common case is that they wink out of existence in the next few generations; only very very few drift to an appreciable frequency). GLIMPSE estimates a posterior probability, which can have some unusual properties at times.<br><br>
To motivate us, let's imagine a polymorphic site that is rare (~2/3 of the sites in the imputation panel), and you have *no* data at the site (or in the vicinity). In this case, GLIMPSE will guess the most likely genotype (spoiler alert, the sample probably doesn't have the rare allele), and the odds that GLIMPSE is correct is staggeringly high. <br>
If 2/3 of our SNVs are rare (which is about what we observe), we probably got most all of these calls correct. Likewise, if we have very little data, the remaining "common" sites may also present a problem. If we require 99% confidence (a posterior probability >99% for us to call a genotype), then the remaining sites will probably mostly be no-calls. In this thought experiment, we just generated a profile that is ~99.9% accurate that has a call-rate of ~2/3, and will be treated by most IBD-segment algorithms as "matching" nearly every person in the population (as the sites needed to tell people apart have been converted to no-calls). <br>
If the above is highly unsettling, good. It should. Users need to be extremely careful of how they treat information.<br>
In case you were wondering, these issues tend to appear when the sample is low-coverage (<1x; this value should be determined internally). Additionally, if you use a Bayes Factor (in the right way), this issue largely disappears.
