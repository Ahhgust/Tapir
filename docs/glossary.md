# Glossary

- **Breadth (of coverage)**. The fraction of the genome that has at least a given number of reads. E.g., Breadth (5x) is the fraction of the genome that has at least 5 reads. See also, Depth, Coverage; estimated using the same sites used to estimate Coverage.

- **Call rate**. For some genotyping panel; what fraction of sites were called (genotyped, not as missing data)?

- **Coverage**. The mean read depth. Often a misused term; we mean it to be the "redundancy of coverage" of [Lander and Waterman](https://doi.org/10.1016/0888-7543(88)90007-9). Note that Lander and Waterman's equation to estimate coverage does NOT apply well to shotgun WGS. In practice, we take coverage as the mean read depth estimated at 10k autosomal sites selected from Illumina's GSA panel. Only reads that pass QC are considered.

- **Depth**. At a site, how many times was it measured? Only reads that pass QC are considered. If reads are selected at random from a person's genome, per Lander and Waterman, Depth is Poisson distributed with lambda=Coverage.

- **Duplicate**. Of a read; multiple measurements of the same DNA (original template) molecule. Duplicates are often either derived from PCR (PCR Duplicates, most common) or optically (Optical Duplicates, less common). In general, when the same molecule is measured more than once, only a single measure is retained.

- **GRCh37**. A version of the reference genome released in 2009. See [link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.13/). GRCh37 (NCBI) is functionally equivalent to hg19 (UCSC). Tapir converts all variant calls to hg19; the genome build used by GEDmatch.

- **GRCh38**. A version of the reference genome released in 2013. See [link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/). GRCh38 (NCBI) is functionally equivalent to hg38 (UCSC). Tapir natively conducts all analyses in hg38.

- **Imputation**. A class of techniques originally developed to infer missing or incorrect data. Modern genetic imputation methods (e.g., GLIMPSE, Beagle) employ genotype refinement. With genotype refinement, uncertain genotype calls (based on their likelihoods, at many linked sites) are converted into (more) certain calls based on both the likelihoods and patterns of LD. Note, all genotype refinement methods are imputation methods, but not all imputation methods are genotype refinement methods

- **Mixture Proportion**. Ratio of the amount of DNA from each contributor included in the total sample that was used for the library preparation process. When estimated using Demixtify [External link](https://doi.org/10.1016/j.fsigen.2023.102980), the mixture proportion is the fraction of reads estimated to come from the minor contributor (2-person mixture scenarios only).

- **Read**. The raw sequence data produced by a sequencing instrument, generated as a result of detection of the sequence of nucleotides in a DNA strand and the translation of that information into digital sequence data.

- **SNP**. Single/Simple/Short nucleotide polymorphism. A SNV with a minor allele frequency >=1%. One of the major aims of the 1000 Genomes Project was to assay ~95% of all SNPs.
- **SNV**. Single/Simple/Short nucleotide variant. A location in the genome where individuals in a population have different alleles. Most SNVs involve a single nucleotide, often as a substitution. There are far more SNVs than SNPs.

- **Template**. Of shotgun sequencing; the (original) thing you sequenced. Ie, the part of a read that neglects adapters and other bits added to enable sequencing. DNA degradation can shorten template lengths.
