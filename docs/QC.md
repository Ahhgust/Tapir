# QC

For definitions, file formats and specifics on the software packages used by Tapir, please see the [Glossary](Glossary.md).

## Sample QC

### Alignment statistics
Tapir runs `samtools flagstat` [External link](http://www.htslib.org/doc/samtools-flagstat.html). <br>
Of note, `samtools flagstat` reports both the raw number of reads and the *number* of properly paired reads; both can be informative on sample performance. 
Note that having a low percentage/number of properly paired reads can be symptomatic of many things, 
included in which are bacterial contamination and adapter dimers (adapter dimers can generate fastq records that are blank; blanks are not properly paired).

### Site-level statistics
One of the most important questions of a sample is how much data do I have? Flagstat information does not (readily) provide that; it tells you how many reads you have. 
If your reads are short (after adapters are removed; short reads may be unrelated to the 2x151 bp sequencing you may have conducts), or your templates are short (similar to the above), or if your reads overlap (read overlaps count as 1 measurement, not 2); all affect how much usable data (i.e., coverage) you have. 
Tapir includes custom software to estimate the breadth and depth of coverage, as well as the number/rate of duplicates.

### Mixture status
Step 2 of Tapir (genotyping) only applies to single-source samples (though we're working on it). We use Demixtify to estimate the mixture proportion and the number of contributors (NoC). Although Demixtify only supports 1 or 2 contributor scenarios (really, hypotheses), 
in practice more than 2 contributors are often identified as mixtures as well. <br>

Of note, genomic techniques can be very sensitive; Demixtify can detect mixtures below 1% (see the original publication 10.1016/j.fsigen.2023.102980). <br> If we have a moderate amount of information, to the forensic scientist the term "drop-in" might instead be a better description. <br>
In my humble opinion, values below 1% are unlikely to be true mixtures (perhaps environmental DNAs, and/or bacterial DNAs are to blame).There edge cases that are exceptions (if we have a 1000x genome, a 1% mixture corresponds to a 10x and 990x genome, respectively; in this very non-forensic example, maybe treating the 10x genome as "real" is not so far-fetched).
<br>
We highly recommend developing an interpretation scheme for mixtures. Here is an example, and note that both the scheme and the cutoffs involved should be empirically determined.<br>
Let's imagine three scenarios:
- The mixture hypothesis does not have significant support
  - Treat as single source (if the *base e* LLR is below some threshold)
  - Treat as inconclusive (otherwise)
- The mixture hypothesis has significant support, but the mixture proportion is quite small
  - Treat as single source.
  - If less than 1%, low-level "contamination" may be to blame (i.e., contamination present that is well below what is traditionally detectable by CE).
  - If less than ~5-10%, genotyping may be largely unaffected, but the mixture hypothesis may be realistic. If the sample *should not* be a mixture, try to investigate what happened.
    - Perhaps evaluate sex-linked data (x/y/mt)
	- If you know ground truth, try pulling out the sites that are homozygous and rare. Do you see the common allele too (and more often than the other two alleles)? (if so, then this strong evidence for a mixture).
- The mixture hypothesis has significant support, but the mixture proportion is NOT small
  - Stop. Collaborate. And listen.
  - Really, figure out what happened. The odds the sample is a true mixture are appreciable. Do you have CE results? What do they say? Are other samples on the same run (apparent) mixtures too? Did you merge BAMs?
  
Additional notes:
Demixtify uses a likelihood-ratio *test* (LRT) [Wikipedia article](https://en.wikipedia.org/wiki/Likelihood-ratio_test) to contrast the single source vs mixture hypotheses. <br>
In short, when we use the LRT we treat (twice the) log likelihood ratio (LLR) as being chi-square distributed with 1 degree of freedom. <br>
Procedurally, if the mixture hypothesis is slightly more likely than the single source, the single source hypothesis may yet be correct.
 At a high level, the mixture hypothesis is more complicated than the single source, and you need to correct for the fact that a more complicated model may fit your data better simply because of chance. <br>. For context, if you remember the difference between an adjusted R2 and R2 (in linear regression), the same core issue is addressed here (namely, that a more complicated model will tend to fit data better than a simpler one, even if the simpler model is **correct**). <br>
One way to address this to set your false discovery rate, or alpha. The original publication suggested 0.01 (though again, this should be evaluated empirically); but know that any alpha that you choose will lead to cases where we incorrectly reject the single source hypothesis.<br>
In practice, I worry about this most when I have very little data.

### Additional reports

- Tapir creates a FastQC report. FastQC looks for likely problems in your Fastq records. Note that Tapir only uses aligned/mapped reads to inform Fastqc. Please see the FastQC website [External link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for instructions on how to interpret this report. I would note, that (IMO) some of the "problems" reported by FastQC (e.g., over-represented sequences) matter for *some* genomic applications (eg, de novo assembly), but likely matter very little for characterizing well-described SNPs/SNVs.
- Tapir also creates a BQSR report. As a reminder BQSR (base quality score recalibration) empirically adjusts base quality scores; base qualities are intended to reflect the probability that the basecall is wrong; in practice, they tend to be a bit inflated. For example, the instrument doesn't know how many PCR cycles were done, and we know PCR makes mistakes at times. Likewise, it cannot account for DNA damage. The overall result is that base qualities are often too high (though not always). In practice, we have also found that poor pre-phasing may lead to (highly) biased base quality. In low-pass scenarios, quality really matters.  Many of these issues are fixed by BQSR. Please see GATK's BQSR manual [External link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR) for more information.


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
          -  Should the coverage be too high (as determined empirically; perhaps by considering blanks), this can indicate a problem with the sample sheet.
	  -  We recommend evaluating both the flagstat and samstat.cov files
-  Undetermined (generated by Tapir, one per run)
   -  Fastq records that failed to demultiplex are placed in the `Undetermined` directory. Poor run quality may make these files big (e.g., if none of adapter 2 generated usable sequence data).
   -  Less likely, if you play with the demultiplexing sensitivity, that too can cause a high demultiplexing failure rate.
      -  One way that the Undetermined files can be big if you have Ns in your adapter sequences. Consider looking at the "I1" and "I2" fastq.gz files. These files are fastqs for the sample indexes themselves. If (for example) you have >2 Ns in them, you'll probably fail to demultiplex the run. 
	  -  Tapir provides `deML` [github link](https://github.com/grenaud/deML), which demultiplexes by maximum likelihood. It isn't implemented in to Tapir as a workflow, but sophisticated users can consider using it.
	  
For troubleshooting, see [General_Troubleshooting.md](General_Troubleshooting.md)	 
