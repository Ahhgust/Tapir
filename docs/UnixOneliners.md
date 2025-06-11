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
If you see "large" values (>0.01x), that may indicate carry over. <br>

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
Okay, the call-rate is a little sticky. But the raw number of genotype calls can be had a la:
For glimpse:
```
for file in Sample_Data/*/Uploads/glimpse2/*snps.tsv.gz; do echo -n $file ' '; zgrep -c -v '#' $file; done
```
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
(note that the above command is ... not so smart. It's comparing HG003 and HG004's data. It maybe provides you an example of a 'null' accuracy')
<br>
where `-t` (HG004 in the above) provides the "true" genotypes, and `-c` provides the comparison. The script provides the following:
-  DepthThreshold,QualityThreshold,PosteriorThreshold
   -  Ignore; these are historical.
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
      -  It will be slightly biased (too small, as call-rates in the ground truth are not 100%, but often nearly so)
<br>

**We highly recommend using the SNP 23andme files (.snps.tsv.gz), otherwise your "accuracy" will be inflated (glimpse).**
