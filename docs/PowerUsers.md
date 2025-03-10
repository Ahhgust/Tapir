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
| DeepVariant | | --model_type=WGS | Available in Tapir; makes VCFs/gVCFs only (no 23andme files) |
| Demixtify  | |  | Estimates the mixture fraction using 586,670 autosomal SNPs from the GSA panel |
| samstats.py* | | | Estimates the mean read depth (an estimate of coverage) using the GSA10K panel |
| bcf223andme.py** | | BCFtools input: -x -q 20 GLIMPSE input: -x -b 1.7 -p 0.99| Prefiltered with bcftools view -e "INFO/OriginalContig!=CHROM" And bcftools norm -m+ |
*Custom software

Note flags affecting threading/multiprocessing/compression are not listed above. That information can be gleaned from the config.yaml file, however.

A more complete set of flags are found in the config.yaml files (defaults to: $TAPIR/snakemakes/config_v_2_standard.yaml.)

## Repeated measurements
Tapir supports repeated measurements on a sample. <br>
Tapir is also incredibly ornery about it. <br>
Repeat after me: If you want merge repeated measurements, <br><br> **They must have the same sample name** <br><br>
e.g.,
It's not <br>
*Sample1_Replicate1* <br>
and <br>
*Sample1_Replicate2* <br>
But it is instead: <br>
*Sample1*
<br>
<br>
Now, if you never ever want to merge the two replicates, the naming above is fine. <br>
<br>
For the record, this isn't some rule that I made up; this rule is a byproduct of how `Read Groups` are used in genomics. See the *Read Groups* section below.
<br>
And for the record, if you wish to merge samples and have mis-named them, the easiest thing to do is just name them correctly (in the sample sheet) and rerun everything.
That is heavy handed, but in fairness, that also makes more work for the computer to do (so really, who cares?).
You can force Tapir to evaluate multiple bams (in the same sample directory) with different/incorrect sample information. Tapir just evaluates: <br>
`*/Bams/Final_Bams/*la.md.bqsr.bam` <br>
so with some creative symbolic linking, you can force it to merge multiple (disparate) files. <br>
Choosing to do this is dangerous and largely untested; on paper, Glimpse should work, while BCFtools (genotyping) would require some changes. 
More importantly, the results will be ( a little bit ) wrong if (and in principle, only if) you sequence the same library multiple times (and fail to specify the right sample identifier).
Said more concretely, you would miss some PCR duplicates...

## Custom Configurations
Tapir supports limited flexibility. Below we describe what things can (and probably should) be configured.
### Configuration files.
Tapir pulls command-line options from a "config" (YAML) file.
A good place to house them is in $TAPIR/configs. We provide a few configs with TAPIR. One is the default 
(`config_v_2_standard.yaml`) and one is for "ultralow" WGS (`config_v_2_ultralow.yaml`). 
Note that the two configs only differ at a single line; the difference (`diff`) is:
```
(tapir) [augustw@hal configs]$ diff config_v_2_standard.yaml config_v_2_ultralow.yaml
186c186
<     glimpse2: "-x -b 1.7 -p 0.99" # (bcf223andme.py) filters for glimpse2; Include the X; Bayes Factor > 1.7; Genotype Posterior > 0.99
---
>     glimpse2: "-x -b 1.7 -p 0.95" # CHANGED HERE (from standard)(bcf223andme.py) filters for glimpse2; Include the X; Bayes Factor > 1.7; Genotype Posterior > 0.95

```
ie, one use a minimum posterior probability of 0.99; ultralow drops that to 0.95 (increasing both the call-rate and the error rate, which you might have to do if you're below the minimum call rate)
(in context, this section of the config describes the filtering performed on `bcf223andme.py`; in words bcf to 23 and me.)

### Local storage
Tapir writes (many) temporary files to local storage devices. In short, if you have an NVMe (or traditional SSD) mounted someplace, update:
<br>
`/tmpdirprefix: "/tmp"` <br>
in the config file(s). 
<br>
Better yet, set `/tmp` to the mount location...
<br>

### Multithreading
In addition, many of the tools that Tapir uses supports some level of multithreading; 
the values in place work reasonably well on our system (a 256-core beast). They may not work well on yours. <br>
For example, you likely don't want the number of threads requested; e.g., `threads: 256`, to exceed what your unix system actually has. 
Consult your cpuinfo file (e.g., `tail -n 40 /proc/cpuinfo`) if you are unsure of how many horses you have.

### Read groups
Genomic tools use *Read Groups* [External link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) to keep track of samples, libraries, and other components of sequencing.
Tapir embeds read-group information when BWA-mem is invoked. Only one parameter, `LB` (the library identifier), can be difficult to glean from a FASTQ record. 
By default, Tapir uses the `"I7_Index_ID"` in the sample sheet and treats that as the identifier. This parameter can be changed in the `config.yaml` file used by Tapir.
<br>
Why does this matter? It mostly doesn't, however, if you are merging multiple BAMs, PCR duplicates can occur both within and between BAMs; 
specifically, if you take the same library and sequence it twice, the LB code should set to be the *same*. If instead you sequence the same sample
using two different library preps, `LB` (within subjects) should be different. Using the `i7_Index_ID` is safe for the first scenario (same library sequenced twice), but it may or may not be true of the second.
That said, the consequences of calling two different libraries the "same" library are pretty minimal; you might artificially call a few more duplicate reads... but that's a pretty minor issue. 
If you want, however, you can include some column in the sample sheet and dedicate it to `LB`; simply modify your config file accordingly!
<br>
Just for utter clarity, Tapir creates (and otherwise assumes) that all BAM files are single-sample (no multisample BAM files are permitted).

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
This pertains to how the `binary` configuration variables are treated by Tapir:
-  .py
   -  files ending in .py are assumed to be python3; `python3` is prepended to the unix call
-  .pl
   -  assumed to be `perl`
-  .R
   -  assumed to be `Rscript`

## Notes on additional utilities
Tapir has quite a few tools present that aren't actively used; they may be incorporated in the future, so maybe this serves as a "sneak peek!"









