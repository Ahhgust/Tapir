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


