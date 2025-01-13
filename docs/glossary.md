# Glossary

-  [**Genomic file formats**](#genomic-file-formats). Tapir uses/creates lots of files. Let's learn about the file formats!
-  [**List of tools**](#list-of-tools) Tapir also uses a lot of tools. What are they? Where do they come from?
-  [**Key terminology**](#key-terminology) Tapir uses a lot of weird words. SNP? SNV? What do these words mean?


## Genomic file formats

- **bam**. Binary alignment map. A compressed (binary) fastq record with alignment information. See [Wikipedia article](https://en.wikipedia.org/wiki/SAM_(file_format))

- **bcl**. Binary base call. Illumina's raw file format. See [External link](https://www.illumina.com/informatics/sequencing-data-analysis/sequence-file-formats.html)

- **fastq**. FASTA ([Wikipedia article](https://en.wikipedia.org/wiki/FASTA_format)) with quality. A text file format (often compressed with gzip). See [Wikipedia article](https://en.wikipedia.org/wiki/FASTQ_format)

- **sam**. an uncompressed version of the bam file format. Good people never store data in sam format.

- **vcf**. Variant call format; often compressed as vcf.gz or bcf. A file format commonly used to describe SNVs. See [Wikipedia article](https://en.wikipedia.org/wiki/Variant_Call_Format).


## List of tools
### Off the shelf tools:
- **BCFtools**. Tool for working with BCF and VCF file formats; Includes tools for filtering and file format interchange (view) and for estimating genotypes in sequence data (call); BCFtools call (as run) does not consider prior information (neither allele frequency nor LD). Available [here](https://github.com/samtools/bcftools)

- **bcl2fastq** Tool for converting BCL files to Fastq files. Works on older Illumina instruments. This is the preferred tool used by Tapir. Available [here](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/downloads.html)

- **bcl-convert** Tool for converting BCL files to Fastq files. Works on newer/all Illumina instruments. This likely requires a custom installation. Available [here](https://www.illumina.com/content/illumina-support/language-master/en/sequencing/sequencing_software/bcl-convert/downloads.html)

- **bwa** Burrows Wheeler Aligner. A tool for aligning Fastq files to the reference genome; Produces sequence alignments (to the reference) in the SAM file format. Available [here](https://github.com/lh3/bwa)

- **bqsr**. Base quality score recalibration. Data reprocessing step that recalculates base qualities based on empirical measures (matches and mismatches at sites generally invariant in human populations). Tapir uses GATK's implementatinon of BQSR, available [here](https://github.com/broadinstitute/gatk/releases)

- **demixtify** A tool that uses biallelic SNVs to detect mixtures and estimate the mixture fraction. Available [here](https://github.com/Ahhgust/Demixtify)

- **fastqc**. Tool for high throughput sequencing which checks the quality control of raw sequence data; includes checks on base quality, read length, and extraneous sequences. Available [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- **glimpse2**. Tool for estimating genotypes from low-coverage whole genome sequencing; Glimpse2 considers LD (and allele frequency, as a prior) when estimating genotypes; Glimpse2 is an imputation method (it is also a genotype refinement method). Available [here](https://github.com/odelaneau/GLIMPSE)

- **samtools** Common tools for working with SAM, BAM, and CRAM file formats; Tools include flagstat (produces summary statistics) and view (file format conversion). Available [here](https://github.com/samtools)

### Custom scripts

- **samstats** Custom script that provides information on segregating sites (selected from the GSA), including read depth and breadth (of coverage), duplication rates and template lengths

- **estimateXPloidy** A simple tool to estimate the ploidy of the X chromosome. BCFtools and GLIMPSE support haploid and diploid calling models; this uses the number of reads on the X vs an autosome (chromosome 7) to estimate the ploidy of the X chromosome (roughly speaking, the biological sex of the sample). In (very) low coverage settings, structural and copy-number variants can affect this tool.


## Key terminology
- **Breadth (of coverage)**. The fraction of the genome that has at least a given number of reads. E.g., Breadth (5x) is the fraction of the genome that has at least 5 reads. See also, Depth, Coverage; estimated using the same sites used to estimate Coverage.

- **Call rate**. For some genotyping panel; what fraction of sites were called (genotyped, not as missing data)?

- **Coverage**. The mean read depth; e.g., 5.2$`\times`$. Often a misused term; we mean it to be the "redundancy of coverage" of Lander and Waterman [External link](https://doi.org/10.1016/0888-7543(88)90007-9). Note that Lander and Waterman's equation to estimate coverage does NOT apply well to shotgun WGS. In practice, we take coverage as the mean read depth estimated at 10k autosomal sites selected from Illumina's GSA panel. Only reads that pass QC are considered.

- **Depth**. At a site, how many times was it measured? Only reads that pass QC are considered. If reads are selected at random from a person's genome, per Lander and Waterman, Depth is Poisson distributed with lambda=Coverage.

- **Duplicate**. Of a read; multiple measurements of the same DNA (original template) molecule. Duplicates are often either derived from PCR (PCR Duplicates, most common) or optically (Optical Duplicates, less common). In general, when the same molecule is measured more than once, only a single measure is retained; other measures are "marked" (and ignored).

- **GRCh37**. A version of the reference genome released in 2009. See [External link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.13/). GRCh37 (NCBI) is functionally equivalent to hg19 (UCSC). Tapir converts all variant calls to hg19; the genome build used by GEDmatch.

- **GRCh38**. A version of the reference genome released in 2013. See [External link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/). GRCh38 (NCBI) is functionally equivalent to hg38 (UCSC). Tapir natively conducts all analyses in hg38.

- **Imputation**. A class of techniques originally developed to infer missing or incorrect data. Modern genetic imputation methods (e.g., GLIMPSE, Beagle) employ genotype refinement. With genotype refinement, uncertain genotype calls (based on their likelihoods, at many linked sites) are converted into (more) certain calls based on both the likelihoods and patterns of LD. Note, all genotype refinement methods are imputation methods, but not all imputation methods are genotype refinement methods

- **Mixture Proportion**. Ratio of the amount of DNA from each contributor included in the total sample that was used for the library preparation process. When estimated using Demixtify [External link](https://doi.org/10.1016/j.fsigen.2023.102980), the mixture proportion is the fraction of reads estimated to come from the minor contributor (2-person mixture scenarios only).

- **Quality**. Quality is commonly used to describe a base call (in a read, how likely is to be called correctly?), an alignment (ostensibly, what is the probability that the read is correctly mapped. Note that this *is not* what mapping quality really is), and a genotype (the likelihood associated with a proposed genotype; often given as a likelihood ratio). Quality is often encoded on a Phred scale [External link](https://en.wikipedia.org/wiki/Phred_quality_score). In most cases, large(r) quality is desired (ie, bigger is better). Notable exceptions include the PL tag in the VCF file format (bigger is... worse?).

- **Read**. The raw sequence data produced by a sequencing instrument, generated as a result of detection of the sequence of nucleotides in a DNA strand and the translation of that information into digital sequence data.

- **SNP**. Single/Simple/Short nucleotide polymorphism. A SNV with a minor allele frequency >=1%. One of the major aims of the 1000 Genomes Project was to assay ~95% of all SNPs.

- **SNV**. Single/Simple/Short nucleotide variant. A location in the genome where individuals in a population have different alleles. Most SNVs involve a single nucleotide, often as a substitution. There are far more SNVs than SNPs.

- **Template**. Of shotgun sequencing; the (original) thing you sequenced. Ie, the part of a read that neglects adapters and other bits added to enable sequencing. DNA degradation can shorten template lengths.

- **Variant**. A place in the genome where the/a sample differs from the reference genome. Can refer to a location, or a location and an allele.