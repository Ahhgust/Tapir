# Glossary

-  [**Genomic file formats**](#genomic-file-formats). Tapir uses/creates lots of files. Let's learn about the file formats!
-  [**List of tools**](#list-of-tools) Tapir also uses a lot of tools. What are they? Where do they come from?
-  [**Key terminology**](#key-terminology) Tapir uses a lot of weird words. SNP? SNV? What do these words mean?


## Genomic file formats

- **bam**. Binary alignment map. A compressed (binary) fastq record with alignment information. See [Wikipedia article](https://en.wikipedia.org/wiki/SAM_(file_format))

- **bcl**. Binary base call. Illumina's raw file format. See [External link](https://www.illumina.com/informatics/sequencing-data-analysis/sequence-file-formats.html)

- **fastq**. FASTA ([Wikipedia article](https://en.wikipedia.org/wiki/FASTA_format)) with Quality. 
A MPS read (DNA sequence) with corresponding base qualities. Often compressed with gzip (as .fastq.gz or fq.gz). See [Wikipedia article](https://en.wikipedia.org/wiki/FASTQ_format)

- **sam**. an uncompressed version of the bam file format. BAM files are large; SAM files are far larger. It is poor practice to keep a file in the `.sam` format.

- **vcf**. Variant call format; often (block) compressed as vcf.gz or bcf. A file format commonly used to describe SNVs. See [Wikipedia article](https://en.wikipedia.org/wiki/Variant_Call_Format).

- **23andMe**. A file format used by *23andMe*, apparently when you download your "raw" data (by which they mean genotypes) this is the format. The format is loosely described [here](http://fileformats.archiveteam.org/wiki/23andMe). Note that Tapir does not add RS numbers. Note that RS numbers, despite what they are defined as, are not (entirely) stable identifiers.

## List of tools
### Off the shelf tools:
- **BCFtools**. Tool for working with BCF and VCF file formats; Includes tools for filtering and file format interchange (view) and for estimating genotypes in sequence data (call); BCFtools call (as run) does not consider prior information (neither allele frequency nor LD). Available [here](https://github.com/samtools/bcftools)

- **bcl2fastq** Tool for converting BCL files to Fastq files. Works on older Illumina instruments. This is the preferred tool used by Tapir. Available [here](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/downloads.html)

- **bcl-convert** Tool for converting BCL files to Fastq files. Works on newer/all Illumina instruments. For Tapir, likely requires a custom installation. Available [here](https://www.illumina.com/content/illumina-support/language-master/en/sequencing/sequencing_software/bcl-convert/downloads.html)

- **bwa** Burrows Wheeler Aligner. A tool for aligning Fastq files to the reference genome; Produces sequence alignments (to the reference; the reads are not aligned to each other) in the SAM file format. Available [here](https://github.com/lh3/bwa)

- **bqsr**. Base quality score recalibration. 
A data reprocessing step that recalculates base qualities based on empirical measures 
(matches and mismatches at sites generally invariant in human populations). 
Tapir uses GATK's implementation of BQSR, available [here](https://github.com/broadinstitute/gatk/releases)

- **demixtify** A tool that uses biallelic SNVs to detect mixtures and estimate the mixture fraction. Available [on Github](https://github.com/Ahhgust/Demixtify)

- **fastqc**. Tool for high throughput sequencing which checks the quality control of raw sequence data; 
includes checks on base quality, read length, and extraneous sequences. Available [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- **GLIMPSE**. Tool for estimating genotypes from low-coverage whole genome sequencing; 
GLIMPSE (v2) considers LD (and allele frequency, as a prior) when estimating genotypes; 
GLIMPSE is an imputation method (it is also a genotype refinement method). Available [here](https://github.com/odelaneau/GLIMPSE)

- **samtools** Common tools for working with SAM, BAM, and CRAM file formats; Tools include flagstat (produces summary statistics) and view (file format conversion). Available [here](https://github.com/samtools)

### Custom scripts

- **samStats.py** Custom script that provides information on segregating sites (selected from the GSA), including read depth and breadth (of coverage), duplication rates and template lengths

- **estimateXPloidy.R** A simple tool to estimate the ploidy of the X chromosome. BCFtools and GLIMPSE support haploid and diploid calling models; this tool uses the number of reads on the X vs an autosome (chromosome 7) to estimate the ploidy of the X chromosome (roughly speaking, the biological sex of the sample). In (very) low coverage settings, structural and copy-number variants can affect this tool.

- **bcf223andme.py** A simple tool to convert a VCF file 
(not from a variant caller, but from Tapir; ie, it must include homozygous reference calls) 
and convert it to a "23andme" file. 
Supports filtering on genotype quality (BCFtools), as well as the posterior probability and Bayes factor (Glimpse).
X chromosome genotypes are "diplotyzed" (converted from haploid to diploid, as necessary). No-calls are dropped.

## Key terminology
- **Bayes factor (BF, as provided by Tapir)**. The posterior odds (of the genotype, GP, transformed to odds), relative to the prior odds (naive, estimated using Hardy Weinberg).
Informally, a posterior probability is composed of two parts; the prior 
(which only cares about what we already know; it does not consider the data/evidence) and the likelihood (which doesn't care about what we knew, but instead reflects the evidence provided by the data).
The BF (as implemented by Tapir) attempts to decouple these two pieces of information. 
A BF>1 implies empirical support for a given genotype call. 
Note that most alleles in a population are rare (<<1%). 
Diabolically, in the complete absence of data one can still accurately infer genotypes.
In short you could use the following algorithm: if the allele is rare, call the genotype homozygous major (often homozygous reference); for other SNPs, simply drop them.
The resulting profile will be accurate (with respect to any individual), and *very* matchy (because most everyone has that profile and the genotypes in the profile are accurate). 
However such markers will have BF<1. It should be noted that the BF described is *a* Bayes factor; not *the* Bayes factor.
<br>
It should be noted that by default, we empose a minimum BF of least 1.7 (for higher quality samples; much higher for lower quality samples). So while the above is *true*, these effects are mitigated by using a BF>1.

- **Breadth (of coverage)**. The fraction of the genome that has at least a given number of reads. E.g., Breadth (5×) is the fraction of the genome that has at least 5 reads. See also, Depth, Coverage; estimated using the same sites used to estimate Coverage.

- **Call rate**. For some genotyping panel; what fraction of sites were called (genotyped, not as missing data)?

- **Coverage**. How many times (on average) do we measure sites? Denoted as "×" (as in *times*), for example: 5.2×. 
Coverage is often a misused term; we mean it to be the "redundancy of coverage" of Lander and Waterman [link to original paper](https://doi.org/10.1016/0888-7543(88)90007-9). 
Note that Lander and Waterman's equation to estimate coverage does NOT apply well to shotgun WGS. In practice, we define coverage as the **mean read depth**. We estimate at 10k autosomal sites selected from Illumina's GSA panel. 
If you want to be formal about things, you can think of coverage (C, that thing you want to estimate) and your estimate of coverage $\hat{\Cmath}$, which is provided by Tapir.
Only reads that pass QC are considered. Before your ask, yes, 10k is actually a very large sample size to estimate a mean.

- **Depth**. At a site, how many times was it measured? Only reads that pass QC are considered. If reads are selected at random from a person's genome, per Lander and Waterman, Depth is Poisson distributed with lambda=Coverage.

- **Duplicate**. Of a read; multiple measurements of the same DNA (original template) molecule. Duplicates are often either derived from PCR (PCR Duplicates, most common) or optically (Optical Duplicates, less common). In general, when the same molecule is measured more than once, only a single measure is retained; other measures are "marked" (and ignored).

- **Genotype posterior probability (GP)**. Provided by Glimpse. An estimate of the genotype posterior probability; considers both the genotype likelihood and a LD-informed prior.

- **Genotype quality (GQ)**. Also see **Quality**. Provided by BCFtools.
The probability that the given genotype call is right (assuming a flat prior). Encoded on a **Phred-scale**, estimated as a log-likelihood ratio (most-likely versus second-most).

- **Genotyping**. See also **Variant calling**. We use BCFtools and Glimpse to perform *genotyping*, wherein we provide the tools with a panel (locations and alleles), and the tools are used to estimate genotype at said sites.

- **GRCh37**. A version of the reference genome released in 2009. See [External link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.13/). GRCh37 (NCBI) is functionally equivalent to hg19 (UCSC). Tapir converts all variant calls to hg19; the genome build used by GEDmatch.

- **GRCh38**. A version of the reference genome released in 2013. See [External link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/). GRCh38 (NCBI) is functionally equivalent to hg38 (UCSC). Tapir natively conducts all analyses in hg38.

- **Imputation**. A class of techniques originally developed to infer missing or incorrect data. Modern genetic imputation methods (e.g., GLIMPSE, Beagle) employ genotype refinement. With genotype refinement, uncertain genotype calls (based on their likelihoods, at many linked sites) are converted into (more) certain calls based on both the likelihoods and patterns of LD. Note, all genotype refinement methods are imputation methods, but not all imputation methods are genotype refinement methods

- **Mixture Proportion**. Ratio of the amount of DNA from each contributor included in the total sample that was used for the library preparation process. When estimated using Demixtify [External link](https://doi.org/10.1016/j.fsigen.2023.102980), the mixture proportion is the fraction of reads estimated to come from the minor contributor (2-person mixture scenarios only).

- **Percent Occupied (of Illumina Sequencing)**. Percentage of clusters with signal (regardless of filter on quality).

- **Phasing (of Illumina Sequencing)**. Rate at which single molecules within a cluster become out of sync with each other during the sequencing process; can be notated as a summary value indicating the percentage of clusters for which the sequencing falls behind.

- **Pre-Phasing (of Illumina Sequencing)**. See Phasing. Instead, the percentage of clusters for which the sequencing jumps ahead.

- **Quality**. Quality is commonly used to describe a base call (in a read, how likely is to be called correctly?), an alignment (ostensibly, what is the probability that the read is correctly mapped. Note that this *is not* what mapping quality really is), and a genotype (the likelihood associated with a proposed genotype; often given as a likelihood ratio). Quality is often encoded on a Phred scale [External link](https://en.wikipedia.org/wiki/Phred_quality_score). In most cases, large(r) quality is desired (ie, bigger is better). Notable exceptions include the PL tag in the VCF file format (bigger is... confusingly... worse).

- **Read**. The raw sequence data produced by a sequencing instrument, generated as a result of detection of the sequence of nucleotides in a DNA strand and the translation of that information into digital sequence data.

- **SAV**. Sequencing Analysis Viewer. Illumina software that allows you to view quality metrics generated by the Real-Time Analysis (RTA) software on Illumina sequencers

- **SNP**. Single/Simple/Short nucleotide polymorphism. A SNV with a minor allele frequency >=1%. 
One of the major aims of the 1000 Genomes Project was to assay ~95% of all SNPs. Additionally, the imputation panel used by Tapir includes the 1000 Genomes Project (as well as data from HGDP).

- **SNV**. Single/Simple/Short nucleotide variant. A location in the genome where individuals in a population have different alleles. Most SNVs involve a single nucleotide, often as a substitution. There are far more SNVs than SNPs.

- **Template**. Of shotgun sequencing; the (original) thing you sequenced. Ie, the part of a read that neglects adapters and other bits added to enable sequencing. DNA degradation can shorten template lengths.

- **Variant**. A place in the genome where the/a sample differs from the reference genome. 
Can refer to a location, or a location and an allele.

- **Variant calling**. See also **Genotyping**. 
A technique used to find all of the locations (really, Variants) in the genome where an individual differs from the reference genome. 
Variant calling (of a single individual) does not *explicitly* include homozygous reference calls. 
Tapir primarily uses **Genotyping**.