# Detailed information


## Command line options
In practice, Tapir simply runs a bunch of UNIX command in a particular order. By default, Tapir uses the following options:
|  Tool     |  Subtool | Additional Parameters | Notes |
| --------  | -------- | --------------------- | ----- |
|  bcl2fastq | -       | --minimum-trimmed-read-length 0 --ignore-missing-bcls --create-fastq-for-index-reads | Enables 0-length FASTQ records |
| bcl-convert | -      | --sample-name-column-enabled true --create-fastq-for-index-reads true --fastq-gzip-compression-level 9 | Set up analogously to bcl2fastq |
| bwa  | mem | -M | -M ensure compatibility with downstream tools (e.g., GATK) |
| sambamba | markdup |  | Functionally equivalent to Picard’s MarkDuplicates tool, but much faster |
| gatk   |  LeftAlignIndels |  | Also serves to merge BAMs|
|        | BaseRecalibrator |  | Masks SNPs/Indels and positions from Woerner et al 2022. This tool is run twice; once to create the recalibration tables, and again to plot the post-BQSR calibration (AnalyzeCovariates) |
|   | ApplyBQSR | |
|   | AnalyzeCovariates | |
|   | LiftoverVCF | --RECOVER_SWAPPED_REF_ALT true -WRITE_ORIGINAL_ALLELES false --WRITE_ORIGINAL_POSITION true --WARN_ON_MISSING_CONTIG true | |
|bcftools | mpileup | -d 512 -q 20 -Q 10 -I -E -a 'FORMAT/DP,FORMAT/AD,FORMAT/SP' | 512 max read depth, 20 min mapping quality, 10 min base quality, additional tags (-a, listed) |
|  | call | -Am -C alleles -P 0. | keep alt alleles, multiallelic calling, limit to the alleles listed, disable the prior. |
| samtools | flagstat | | |
| GLIMPSE2 | phase/ligate | --mapq20 | Chunk sizes are set to twice the default. |
| DeepVariant | | --model_type=WGS | Available in Tapir; makes VCFs/gVCFs only |
| Demixtify  | |  | Estimates the mixture fraction using 586,670 autosomal SNPs from the GSA panel |
| samstats.py* | | | Estimates the mean read depth (an estimate of coverage) using the GSA10K panel |
| bcf223andme.py** | | BCFtools input: -x -q 20 GLIMPSE input: -x -b 1.7 -p 0.99| Prefiltered with bcftools view -e "INFO/OriginalContig!=CHROM" And bcftools norm -m+ |
* Custom software

Note flags affecting threading/multiprocessing/compression are not provided.









