# Files and formats
Tapir makes a ridiculous number of files.
This is my attempt to document all non-obvious and/or non-trivial files.

-	[**Standard genomic file formats](#Standard-Formats) Links to various standard genomic file formats.
-	[**Suffixes](#Suffixes) File names tell you how the file was processed. Learn more!
-	[**Run_Data](#Run-Data) Custom files from Tapir (Step 1)
-	[**Sample_Data](#Sample-Data) Custom files from Tapir (Step 2)

Forgot what a term means, consult the [Glossary](Glossary.md)

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
Different genomic tools read in information and produce information. In general, the steps that lead to some file can be read from left to right, meaning that tasks are added to the end of the file name.
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
Files in the `Run_Data` directory are run specific.

### Flagstat
See `Sample_Data`

### Samstats
See `Sample_Data`

### Samstats.cov
See `Sample_Data`

## Sample Data
(Step 2 of Tapir) <br>
Files in the `Sample_Data` directory are sample specific.
-	The "Combined" bam file:
	-  *merged.md.bam
	-  Is either a symlink to the *single* bam file for this sample
		-  OR
	-  Is a the file has been merged and duplicates marked again ("remarked") (NOT a symlink)

### Xploidy
In the VCFs directory, you will find an "xploidy" file. Tapir uses it to estimate the ploidy of the X chromosome, and then use that information to genotype the X (haploid or diploid).
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

Refers to Demixtify (v1). In brief, demixtify assesses two hypotheses: the sample is a DNA mixture (2 unknown persons) versus a single source sample (1 unknown).
It estimates the two corresponding likelihoods using the biallelic SNVs found in: `$TAPIR/resources/GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.vcf.gz`) 
Demixtify is available on (my) github: [Link](https://github.com/Ahhgust/Demixtify)
The file formats are described [Here](https://github.com/Ahhgust/Demixtify/blob/main/MFfile.md)
<br>
Note, Tapir uses V1 of Demixtify. It is only used for mixture identification, not deconvolution.
<br>
In practice, we find that testing the two stated hypotheses is sufficient; eg, 3 person mixtures are still much more likely under the two (versus) one contributor hypothesis.
<br>
In addition, Demixtify is very efficient; even a small fraction of another person's reads(~0.5%, perhaps less) provides enough evidence for a DNA mixture.
We use the following interpretation scheme:
- 	The LLR < 6.635/2
	-	We call the sample single source
		-	6.635 uses an alpha of 0.01
		-	Use a chi-square table (1 degree of freedom) to adjust
	- 	If the LLR is negative, then the single source hypothesis is more likely
	-	If it's positive, you might expect it to be positive just by chance alone
		-	That's why we use a likelhood ratio *test*.
-	The mixture proportion (mp) <= 5% (and the mixture hypothesis has significant support)
	- 	We call the sample *effectively* single source
		-	5% is an aribtrary cutoff, where the major contributor tends to be correctly genotyped even if it's truly a highly imbalanced mixture.
	-	If the mp is estimated to be 0.5% (0.005)
		- 	Don't worry about it. The true *mp* can be truly trivial `<<0.005`
	-	Otherwise
		- 	Maybe consider diagnosing the situation. You have alleles of an unexpected kind in your data. Bacterial contamination can (in principle) cause this, for example.
-	The mixture proportion (mp) > 5% (and the mixture hypothesis has significant support)
	- 	We call this a mixture.
		- 	Perhaps you merged some biological replicates together?
			-	(mis-labeling?)
		-	Or perhaps not.

**Sequence error** <br>
Demixtify also estimates the sequence error rate based on the prevalence of alleles other than the two potential alleles (e.g., at C/T sites, Demixtify will tabulate the fraction of G and A alleles to estimate an overall rate of error).
<br>A high sequence error rate can be indicative of real problems; e.g., runs with poor phasing/prephasing. If there is *bias* in which alleles arise (ie, most tools assume that sequence error is uniform; when this assumption is violated that is a kind of bias), then some of your downstream tools may misbehave. BQSR fixes some of these issues, but not all. 
