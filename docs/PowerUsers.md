# Detailed information
Tapir can (and should) be tuned to work well on your system.

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
*Custom software

Note flags affecting threading/multiprocessing/compression are not listed above. That information can be gleaned from the config.yaml file, however.

A more complete set of flags are found in the config.yaml files (snakemakes/config_v2_1.yaml.)

## Custom Configurations
Tapir uses a config.yaml file (with some version number).
The config yaml files provide come useful constants, including:
<br>
`/tmpdirprefix: "/tmp"` <br>
Ideally, /tmp (or whatever you set this to) is local fast storage (NVMe is best)
<br>
In addition, each tool supports some level of multithreading; the values in place work reasonably well on our system. They may not work well on yours. For example, you likely don't want the number of threads requested; e.g., `threads: 256`, to exceed what your unix system actually has. Consult your `/proc/cpuinfo` file if you are unsure of how many horses you have.

### Read groups
Genomic tools use *Read Groups* [External link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) to keep track of samples, libraries, and other components of sequencing.
Tapir embeds read-group information when BWA-mem is invoked. Only one parameter, `LB` (the library identifier), can be difficult to glean. 
By default, Tapir uses the `"I7_Index_ID"` in the sample sheet and treats that as the identifier. This parameter can be changed in the `config.yaml` file used by Tapir.
<br>
Why does this matter? It mostly doesn't, however, if you are merging multiple BAMs, PCR duplicates can occur both within and between BAMs; 
specifically, if you take the same library and sequence it twice, the LB code should set to be the *same*. If instead you sequence the same sample
using two different library preps, `LB` (within subjects) should be different. Using the `i7_Index_ID` is safe for the first scenario (same library sequenced twice), but it may or may not be true of the second.
That said, the consequences of calling two different libraries the "same" library are pretty minimal; you might artificially call a few more duplicate reads... but that's a pretty minor issue. 
If you want, however, you can include some column in the sample sheet and dedicate it to `LB`; simply modify your config file accordingly!

## Notes on paths:
the config file assumes one of three kind of paths when referring to some file:
-  no directories
   -  e.g., binary: "samtools"
   -  In this case, the binary is assumed to come from the user environment
   -  In the samtools example, samtools is assumed to be in your PATH.
-  an absolute directory
   -  e.g., binary: "/foo/bar/myfile.py"
   -  the same as above
-  a relative directory
   -  e.g., binary: "bin/sambamba"
   -  tapir assumes that this file is owned by Tapir and the relevant prefix is added (converting it to an absolute path in practice)
   
## Notes on scripts
-  .py
   -  files ending in .py are assumed to be python3; `python3` is prepended to the unix call
-  .pl
   -  assumed to be `perl`
-  .R
   -  assumed to be `Rscript`

## Notes on additional utilities
Tapir has quite a few tools present that aren't actively used; they may be incorporated in the future, so maybe this serves as a "sneak peek!"









