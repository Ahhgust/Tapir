# QC

For definitions, file formats and specifics on the software packages used by Tapir, please see the [Glossary](Glossary.md).

## Sample QC

### Alignment statistics
Tapir runs `samtools flagstat` [External link](http://www.htslib.org/doc/samtools-flagstat.html). <br>
Of note, `samtools flagstat` reports both the raw number of reads and the *number* of properly paired reads; both can be informative on sample performance. 
Note that having a low percentage/number of properly paired reads can be sympomatic of many things, 
included in which are bacterial contamination and adapter dimers (adapter dimers can generate fastq records that are blank; blanks are not properly paired).

### Site-level statistics
One of the most important questions of a sample is how much data do I have? Flagstat information does not (readily) provide that 
(what if your reads are short? or your reads are long and your templates are short?). 
Tapir includes custom software to estimate the breadth and depth of coverage, as well as the number/rate of duplicates.

### Mixture status
Step 2 of Tapir (genotyping) only applies to single-source samples. We use Demixtify to estimate the mixture proportion and the number of contributors (NoC). Although Demixtify only supports 1 or 2 contributor scenarios (really, hypotheses), 
in practice more than 2 contributors are identified as mixtures as well. <br>
Of note, genomic techniques can be very sensitive; Demixtify can detect mixtures below 1%. 
In practice, values below 1% are unlikely to be true mixtures (perhaps environmental DNAs, and/or bacterial DNAs are to blame). 
We recommend a three-tiered interpretation scheme:
- The mixture hypothesis does not have significant support
  - Treat as single source
- The mixture hypothesis has significant support, but the mixture proportion is quite small
  - Treat as single source.
  - If less than 1%, low-level "contamination" may be to blame (i.e., contamination present that is well below what is detectable by CE).
  - If less than ~5-10%, genotyping will be large unaffected, but the mixture hypothesis may be realistic. If the sample *should not* be a mixture, try to investigate what happened.
    - Perhaps evaluate sex-linked data (x/y/mt)
	- If you know ground truth, try pulling out the sites that are homozygous and rare. Do you see the common allele too (and more often than the other two alleles)? (if so, then this strong evidence for a mixture).
- The mixture hypothesis has significant support, but the mixture proportion is NOT small
  - Stop. Collaborate. And listen.
  - Really, figure out what happened. Odds are, you found yourself a mixture.
  
Additional notes:
Demixtify uses a likelihood-ratio *test* [Wikipedia article](https://en.wikipedia.org/wiki/Likelihood-ratio_test) to contrast the single source vs mixture hypotheses. <br>
In short, we treat (twice the) log likelihood ratio (LLR) as being chi-square distributed with 1 degree of freedom. Procedurally, if your data only weakly support the 
mixture hypothesis, that may not be enough evidence to support the hypothesis. Intuitively, the mixture hypothesis is more complicated, and you need to correct for the fact that a more complicated model may fit your data better simply by chance. <br>
We use an alpha of 0.01 (this helps control the false discovery rate), but know that any alpha that you choose will lead to cases where we incorrectly reject the single source hypothesis.
In practice, you should worry about this most when you have very little data.

### Additional reports

- Tapir creates a FastQC report. FastQC looks for likely problems in your Fastq records. Note that Tapir only uses aligned/mapped reads to inform Fastqc. Please see the FastQC website [External link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for instructions on how to interpret this report. I would note, that (IMO) some of the "problems" reported by FastQC (e.g., over-represented sequences) matter for *some* genomic applications, but likely matter very little for characterizing well-described SNPs/SNVs.
- Tapir also creates a BQSR report. As a reminder BQSR (base quality score recalibration) empirically adjusts base quality scores; base qualities tend to be a bit "optimistic", and cannot account for PCR error/DNA damage; likewise, poor pre-phasing results tend to lead to (highly) biased base quality. In low-pass scenarios, quality really matters.  Many of these issues are fixed by BQSR. Please see GATK's BQSR manual [External link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR) for more information.


## Higher-order QC

Errors and the like can occur on many levels; the most basic is at the level of the sample (above).
Higher-order effects matter too-- a run with poor phasing impacts all samples; likewise, if you add a boatload of dimers (in a negative, for example), it is argued that that will obliterate the run. Please consider the following:
-  Run-level QC
   -  See:
      -  [External link](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-troubleshooting-list/000001911) on adapter dimers.
      -  [External link](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000003874) short inserts.
   -  We recommend using Illumina's (on-site) metrics on the run. In particular, you may find the following metrics to be important:
      -  Pre/Post phasing and/or Base misincorporation rates
      -  % of clusters >= Q30
      -  % of clusters passing filters
         - Really, the *pair* of these two values, in a scatterplot.
      -  % of clusters occupied. 
-  Offtargets (directory created by Tapir; one per run)
   -  At a high level, this can be something simple and common; you used the wrong sample sheet, or you told Tapir that you expect to see UDI25, but instead you used UDI23. We call these "Offtarget" results.
       -  We recommend using *all* sample indexes when demultiplexing.
       -  See the "Offtargets" results
          -  Should the coverage be too high, this can indicate a problem with the sample sheet.
	  -  We recommend evaluating both the flagstat and samstat.cov files
-  Undetermined (generated by Tapir, one per run)
   -  These results correspond to fastq records that failed to demultiplex. Poor run quality may make these files big.
   -  Less likely, if you play with the demultiplexing sensitivity, that too can cause a high demultiplexing failure rate.
      -  One way that the Undetermined files can be big if you have Ns in your adapter sequences. Consider looking at the "I1" and "I2" fastq.gz files. These files are fastqs for the sample indexes themselves. If (for example) you have >2 Ns in them, you'll probably fail to demultiplex the run. 
	  -  Tapir provides `deML` [github link](https://github.com/grenaud/deML), which demultiplexes by maximum likelihood. It isn't implemented in to Tapir as a workflow, but sophisticated users can consider using it.
	  
For troubleshooting, see [General_Troubleshooting.md](General_Troubleshooting.md)	 
