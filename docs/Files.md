# Files and formats
Tapir makes a ridiculous number of files.
This is my attempt to document all non-obvious and/or non-trivial files.

-	[Standard genomic file formats](#Standard-Formats) Links to various standard genomic file formats.
-	[Suffixes](#Suffixes) File names tell you how the file was processed. Learn more!
-	[Run_Data](#Run-Data) Custom files from Tapir (Step 1)
-	[Sample_Data](#Sample-Data) Custom files from Tapir (Step 2)

Forgot what a term means, consult the [Glossary](Glossary.md)

## What is (and is not) kept
Tapir creates many files; not all are retained. Some are explicitly not retained; 
eg, when we genotype we break the genome into pieces, genotype (in parallel), and then stitch them back together again.
These pieces are not retained (and are written to /tmp by default).
We also create intermediate files (eg, a bam file exists before and after duplicates are marked).
We use the `temp` directive in snakemake to do so. These files can be (forcibly) retained. Just add `--notemp` to you snakemake command line.

## Standard Formats
Tapir uses many standard file formats. These include:

-  Fastq
   -  [Format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
   -  [Explained in plain English](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211)
-  SAM/BAM 
   -  [Format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
   -  [Wikipedia article](https://en.wikipedia.org/wiki/Binary_Alignment_Map)
-  VCF/BCF
   -  [Format specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
   -  [Wikipedia article](https://en.wikipedia.org/wiki/Variant_Call_Format)
-  BED
   -  [Format specification](https://genome.ucsc.edu/FAQ/FAQformat.html)
   -  [Wikipedia article](https://en.wikipedia.org/wiki/BED_(file_format))
   -  Note: the BED file format can have different columns (other than the first 3, which are always `chrom, start and stop`).

Note that SAM/VCF are text (uncompressed), while BAM, BCF and .vcf.gz are compressed, and often indexed (common extensions for the indexe files are: .bai, .csi, .tbi) so that the files can be queried quickly. 


### Suffixes 
Genomic tools read in information and then produce information. We also track these tools using the file name. Ideally, a good file name should tell you (at the level of 10,000 feet) what was done to that file.
In general, the steps that lead to some file can be read from left to right, meaning that tasks are added to the end of the file name.
<br>
The extension:
```
.la.md.bqsr.bam
```
Is read in plain English as:
-  Left-aligned (la) around indels (a genomic standard)
-  marked duplicates (md)
-  base quality score recalibrated (bqsr), by GATK
-  bam (SAM/BAM above).

which tells you both what was done, and the order (first `la`, then `md`, then `bqsr`).

# Custom Files

## Run Data 
(Step 1 of Tapir) <br>
Files in the `Run_Data` directory are run- (as opposed to sample-) specific.

### Files for reproducibility
-  `Run_Fastqs/*your_run*/SampleSheet.csv`
   -	This is the sample sheet you used when you ran `bcl2bams.smk`. It is written just after the FASTQs are extracted.
-  `Run_Fastqs/*your_run*/version.txt`
   -	The version of `bcl2bams.smk` that was used.
-  `Run_Fastqs/*your_run*/config.yaml`
   -	This is the config file used.


### Files related to FASTQ extraction
-  `Run_Fastqs/*your_run*/Fastqs/Stats/AdapterTrimming.json`
   -	Created by `bcl2fastq`. This gives trimming
-  `Run_Fastqs/*your_run*/Fastqs/bcl2fastq.outerr
   -	Also createrd by bcl2fastq. This is the standard output/standard error of the run
   -	Note it is not in a 'Logs' directory.

### FASTQs
Fastqs are stored here:
`Run_Fastqs/*your_run*/Fastqs/*fastq.gz`
Tapir (by default) retains all FASTQs. It also extracts your indexes (necessary input for `deML`). 
The "index" fastqs have `_I1_` or `_I2_` in their file name, while the read data have `_R1_` or `_R2_` 
The *size* of the fastqs can be a rough indicator of performance. <br>
In particular, FASTQs from the unclaimed indexes (typically UDI/UDP as a prefix) can be diagnostic--
 if these files are very large, you likely have an error in your sample sheet (or worse; contamination).

## Sample Data
From step 1 of Tapir. <br>
### Sample_Data/Undetermined/
This directory has all of the data related to the "Undetermined" fastqs. <br>
i.e., data that failed to demultiplex. Tapir provides limited reporting on these files.
Relevant data are in: 
```
Sample_Data/Undetermined/*your_run*/Reports/
```
Which includes samstats/flagstat files (read depth information may be helpful), as well as the fastqc report.
Tapir does not run BQSR on these files (BQSR is slow and serves little point here).
If you have poor run performance, it may pay to check if demultiplexing performed poorly. 
The configuration files do let you modify the stringency when demultiplexing (bcl2fastq; you can add `--barcode-mismatch 3` to the `params` directive).
Do so at your own risk, however. (eg, mismatch of 1 is the default; you likely won't be able to demultiplex 96 samples with anything larger than 1).
You can also check out the FASTQ files (the indexes themselves; eg, from the *Fastqs/* directory:
```
ls -lSr UD*_I[12]_*gz | tail
```
will show you the 10 largest barcode files; while:
```
zcat *biggestIndexFile* | head -n 40
``` 
will show you the first 10 FASTQ records in that file (each fastq record is 4 lines long; hence 40). Demultiplexing failure can happen, for example, if the sample indexes have a lot of Ns in them.
`deML` may be able to save some of the data, but note that `deML` has not been integrated into Tapir. (ie, this fix requires a fair amount of expertise)

<br>
You can also use tapir to see what indexes are appearing in the Undetermined files  a la:
```
zcat Undetermined_S0_L001_R1_001.fastq.gz | head -n 400000 | python3  $TAPIR/bin/getBarcodes.py|tail -n 30
```

Which reports the index (pair, observed), a raw count and a percentage of reads associated with that barcode pair.
<br>
Remember, if you're working with a NovaSeq, an index that is all "G" is really just missing data (in this two-dye system, G is no signal).


### Sample_Data/Offtargets/
We recommend extracting all indexes, regardless of their intended use. We call these data "offtargets". They are treated as quasi-samples. It is important to evaluate the throughput *in every run*.
You can get the mean read depth (from the `Sample_Data/Offtargets/*your_run*/Reports/` directory) using the below:
```
$TAPIR/bin/fcat.pl -h *cov | cut -f2-3
```
which concatenates (combines) all the coverage estimates. An offtarget file that has a lot of data in it probably indicates a mislabeled sample in your sample sheet.
It can also indicate carryover/contamination.

### Flagstat
See `Sample_Data`

### Samstats
See `Sample_Data`

### Samstats.cov
See `Sample_Data`

### BQSR report
Tapir runs GATK's BQSR [External Link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)
and generates a recalibration report (item 5. in the above link). The report is named:
```
Sample_Data/*your_sample*/*your_run*/Reports/*bqsr_summary.pdf
```
And it provides a high-level description on the extent and type of quality score recalibration that occurred.

### Fastqc reports
Tapir provides two fastqc reports; one for all of read1, and another for all of read2 (r1_fastqc.zip, r2_fastqc.zip). 
These are sample+run specific (so the same directory as the bqsr reports). <br>
Here's an [External link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to what these reports look like. <br>
Of note, we run fastqc on the aligned reads from the BAMs, which is a bit atypical (but generating 4x2=8 reports per sample; r1 and r2, 4 lanes) is crazy stupid and redundant. This is less information at least :)

## Sample Data
(Step 2 of Tapir) <br>
Files in the `Sample_Data` directory are sample specific.
-	The "Combined" bam file:
	-  *merged.md.bam
	-  Is either a symlink to the *single* bam file for this sample
		-  OR
	-  Is a the file has been merged and duplicates marked again ("remarked") (NOT a symlink)

### Xploidy
In the VCFs directory, you will find an "xploidy" file. Tapir uses it to estimate the ploidy of the X chromosome, it then uses that information to genotype the X (haploid or diploid).
<br>
Tapir does so by considering two hypotheses: the individual is AA:XX (eg, female), or the individual is AA:X (eg, male). The number of reads mapping to those chromosomes lets us assess these hypotheses.
The result from the hypothesis testing is written to:
```
cat Combined.la.md.bqsr.xploidy.txt
```
which was estimated using: `$TAPIR/bin/estimateXPloidy.R`

The file is a headerless tsv; it reports (in this order), 
-	The log-likelihood of AA:XX hypothesis
	-  (Equal ploidy between chromosome X and the autosomes)
-	The log-likelihood of the AA:X hypothesis	
	-  (2:1 ploidy of the autosomes vs X)
	-  Both likelihoods are from a binomial distribution, comparing read counts on the autosomes (chr7) and the X chromosome (non-par, with matching coordinates between the chromosomes).
-	The log likelihood ratio (LLR, estimated as column 1 - column 2)
-   And the estimated X ploidy (1 iff LLR < 0, 2 otherwise)

### Flagstat
This is the direct output of `samtools flagstat` on the `.la.md.bam` (prior to BQSR) bam file. In brief, flagstat gives you some summary information on the *reads* in your BAM file. (how many, what alignment properties). <br>
See the samtools documentation [here](https://www.htslib.org/doc/samtools-flagstat.html)

### Samstats
Flagstats files are about (how many) reads you generated. Samstats instead focuses on sites. It considers 10k autosomal sites (from the GSA), and looks at the properties of the reads there.
It estimates the following summary statistics:
-	Depth (ie, the read depth)
	-	At a site, how many reads do I see (after QC)
-	DepthWithDups (the read depth, including dups)
	-	At a site, how many reads do I see (prior to QC; this counts PCR/optical duplicates)
-	Template Length
	-	How long is the template (ie, the thing you sequenced). Taken as the TLEN field in the BAM file format.
-	ReadLength
	-	How long is the read. (IMO, the tlen field is much more meaningful)
-	Nsites
	-	Trivial; how many sites do we consider. By default, 10,0000
	-	Associated Index is always 0 (can be ignored).

For each summary statistic, Samstats reports the exact distribution (marginally; for each summary ignoring the others).<br>
The samstats file is a TSV, with the following headers:
```
Index\tCount\tLabel
```
Where Label is the summary statistic, and `Count` is how many times you observed a given `Index`. eg,

```
[augustw@hal Final_Reports]$ fgrep Depth Combined.la.md.bqsr.merged.md.samstats | head -n 3
0       5839    Depth
1       2490    Depth
2       877     Depth
```

Which means we have 5839 instances of a read depth (after QC) of 0, 2490 cases of a read depth of 1 and 877 instances of a read depth of 2...


### Samstats.cov
This is a digest of the `samstats` file. It estimates the mean and variance of the samstats summaries. It also estimates Breadth (1,5,10,20x). eg, Breadth5x is fraction of sites that have a read depth >= 5. <br>
We also provide "FracSampledOfLibrary", which looks at the duplication rate (taken as the Mean(DepthWithDups)/Mean(Depth)), and estimates the fraction of the library that was sequenced. <br>
This is a naive implementation of the `EstimateLibraryComplexity` tool of GATK/Picard.
See [link](https://gatk.broadinstitute.org/hc/en-us/articles/360037591931-EstimateLibraryComplexity-Picard)
<br>
Note, this implementation assumes that all duplicate reads arise from PCR, while the incredibly slow Picard version tries to separate those variables. <br>
In practice, the value provided shouldn't be taken at face value, but it does tell you a lot about if it's worth it to re-sequence that library again.

### Demix

Refers to Demixtify (v1). In brief, demixtify assesses two hypotheses: the sample is a DNA mixture (2 unknown persons) versus single source sample (1 unknown).
It estimates the two corresponding likelihoods using the biallelic SNVs found in: `$TAPIR/resources/GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.vcf.gz`) 
Demixtify is available on (my) github: [Link](https://github.com/Ahhgust/Demixtify)
The file formats are described [Here](https://github.com/Ahhgust/Demixtify/blob/main/MFfile.md)
<br>
Note, Tapir uses V1 of Demixtify. It is only used for mixture identification, not deconvolution.
<br>
In practice, we find that testing the two stated hypotheses is sufficient; eg, 3 person mixtures are still much more likely under the two (versus) one contributor hypothesis.
<br>
In addition, Demixtify is very efficient; even a small fraction of another person's reads (~0.5%, perhaps less) can provide enough evidence that a sample is a DNA mixture. (CE is far less sensitive)
We use the following interpretation scheme:
- 	The LLR < 6.635/2
	-	This is an application of a chi-square approximation to the likelihood ratio (LR) *test*. It is a significance test on the LR.
	-	We call the sample single source
		-	6.635 uses an alpha of 0.01
		-	Use a chi-square table (1 degree of freedom) to adjust
	- 	If the LLR is negative, then the single source hypothesis is more likely
	-	If it's positive, you might expect it to be positive just by chance alone
		-	That's why we use a likelhood ratio *test*.
-	The mixture proportion (mp) <= 5% (and the mixture hypothesis has significant support)
	- 	We call the sample *effectively* single source
		-	5% is an arbitrary cutoff, where the major contributor tends to be correctly genotyped even if it's truly a highly imbalanced mixture.
	-	If the mp is estimated to be 0.5% (0.005)
		- 	Don't worry about it. The true *mp* can be truly trivial `<<0.005` and what is observed may be explained by a few rogue reads.
		-   All that can be really said is that 0.5% is more likely than 0%. The true MF may be 0.9% or 0.000001%.
	-	Otherwise
		- 	Consider diagnosing the situation. You have alleles of an unexpected kind in your data. Bacterial contamination can (in principle) cause this, for example.
		-	In *some* hair samples we have noticed that the MP can be small (~1-2%) and have significant support.
		-	Consider checking the distribution of read depths; are there outliers?
			- 	For low coverage samples (<1x), you can check the "Poisson-ness" of the sample:
			-	`1.0 - dpois(0, 0.060)`
			-	(where the mean read depth is 0.060).
			-	In words, `dpois(0, 0.060)` gives the probability of observing 0 reads given a 0.060x coverage genome (under a Poisson model; not a bad model, especially in low coverage scenarios).
			-	1.0- then flips it; it's the probability of observing 1+ reads.
			-	This is an estimate of the Breadth at 1X. If the values are similar, then you're happy.
			-	If not, you likely have some form of drop-in
				- 	In essence, some places in the genome have too many reads.
-	The mixture proportion (mp) > 5% (and the mixture hypothesis has significant support)
	- 	We call this a mixture.
		- 	Perhaps you merged some biological replicates together?
			-	(mis-labeling?)
		-	Or perhaps not.

**Sequence error** <br>
Demixtify also estimates the sequence error rate based on the prevalence of alleles other than the two potential alleles (e.g., at C/T sites, Demixtify will tabulate the fraction of G and A alleles to estimate an overall rate of error).
<br>A high sequence error rate can be indicative of real problems; e.g., runs with poor phasing/prephasing. If there is *bias* in which alleles arise (ie, most tools assume that sequence error is uniform; when this assumption is violated that is a kind of bias), then some of your downstream tools may misbehave. BQSR fixes some of these issues, but not all. 

### snpsummary

The "snpsummary" file contains summary statistics on genotypes. Summaries can either be "unfiltered" (no QC, aka, no filtering on genotype quality (bcftools) or the genotype posterior (glimpse)).
<br>
Likewise, statistics can be filtered (*after* QC filtering) <br>
All summaries are autosomal (only). <br>
e.g., the number of autosomal genotypes:
```
zegrep -c -v '#|X' Uploads/glimpse2/Combined.la.md.bqsr.glimpse2.hg19.23andme.snps.tsv.gz
```
Should match the number of genotype calls <br>
(the 3rd column from: `perror_filt_mean`)

The summaries are as follows:
-  Heterozygosity (filtered)
   -  Observed (`observed_het`)
      -  This is the fraction of sites called as heterozygous (2nd column)
	  -  And the total number of sites considered (3rd column)
   -  Expected (e.g., `AF_nfe_ehet`)
      -  Expectations are per population
 	     -  See gnomad (v3.1.2) [External Link](https://gnomad.broadinstitute.org/news/2017-02-the-genome-aggregation-database/)
	     -  `AF_nfe_ehet` uses allele frequencies from Non-Finnish Europeans
	  -  Expectations are theta-corrected 
	     -  2pq-(1-fst) 
		 -  Fst is estimated per locus, using the (somewhat naive) estimator of Hudson [External Link](doi.org/10.1093/genetics/132.2.583)
-  Rare homozygotes (filtered)
   -  See heterozygotes
   -  Of the sites with a minor allele frequency (MAF) < 10%
      -  What fraction are homozygous minor?  
      -  And how many might be expected?
         -  (Note, also theta corrected)	  
   -  In general, GEDmatch's IBD segment tools are pretty insensitive to heterozygotes. 
      -  But they are quite sensitive to rare homozygotes
-  Error rate
   -  Unfiltered (`perror_unfilt_mean`) *AND*
   -  Filtered (`perror_filt_mean`)
      -  Uses either genotype quality or the genotype posterior
-  Quality bins (unfiltered)
   -  Using either the genotype quality (PL) *OR*
   -  The genotype posterior (Phred-encoded) (GQ)
      -  values clamped to Q50
   -  Make a histogram of counts
      - e.g., `phred_qbin_GP 0_5 100` would say that there are 100 instances of genotype quality scores between 0 and 5.
	  