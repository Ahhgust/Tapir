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
| Demixtify  | |  | Estimates the mixture proportion and likelihoods for the single-source vs two person mixture hypothesis using 586,670 autosomal SNPs from the GSA panel. |
| samstats.py* | | | Estimates the mean read depth (an estimate of coverage) and related statistics using the GSA10K panel |
| bcf223andme.py** | | BCFtools input: -x -q 20 GLIMPSE input: -x -b 1.7 -p 0.99| Prefiltered with bcftools view -e "INFO/OriginalContig!=CHROM" And bcftools norm -m+ |
*Custom software

Note flags affecting threading/multiprocessing/compression are not listed above. That information can be gleaned from the config.yaml file, however.

A more complete set of flags are found in the config.yaml files (defaults to: $TAPIR/snakemakes/config_v_2_standard.yaml.)

## Repeated measurements (merging runs)
Tapir supports repeated measurements on a sample. <br>
Tapir is also incredibly particular about it. <br>
Repeat after me: If you want merge repeated measurements, <br><br> **They must have the same sample name** <br><br>
e.g.,
In the `SampleSheet`, it's not <br>
*Sample1_Replicate1* <br>
and <br>
*Sample1_Replicate2* <br>
But it is instead: <br>
*Sample1*
<br>
<br>
Now, if you never ever want to merge the two replicates, the first set of naming above is fine. <br>
<br>
For the record, this isn't some rule that I made up; this rule is a byproduct of how `Read Groups` are used in genomics. See the *Read Groups* section below.
<br>
And also for the record, if you wish to merge samples and have mis-named them, the easiest thing to do is just name them correctly (in the sample sheet) and rerun everything.
That is heavy handed, but in fairness, that also makes more work for the computer to do (so really, who cares?).
You can force Tapir to evaluate multiple bams (in the same sample directory) with different/incorrect sample information. Tapir just evaluates: <br>
`*/Bams/Final_Bams/*la.md.bqsr.bam` <br>
so with some creative symbolic linking, you can force it to merge multiple (disparate) files. <br>
Care must be taken, however; the read-group information must be constructed so as to reflect the right library and sample. <br>
Know that GLIMPSE ignores read-group information (and always emits a single genotype in a single-sample VCF file), whereas BCFtools will generate multi-sample VCF (if that's what the read-groups describe); note that multi-sample VCFs will break Tapir.


### Forcibly merging samples within or between runs
Tapir is very conservative in how it merges samples. The sample names need to be exactly the same; by construction, they also need to be run on different flowcells (otherwise we'd have the same name twice, which is a no-no when it comes to how we name files!)
Note that this is not Tapir's fault; this is just how genomics works.
<br>
Sometimes we need to merge runs across samples. For example, let's say we have two hairs: `hair1` and `hair2`. We'd like to think they came from the same person, but maybe they didn't. If neither one gave us enough information to do genotyping, maybe both of their results together is enough!
<br>
Tapir provides *limited* support for doing such things. 
Right now that support is limited to GLIMPSE (though the BCFtools mpileup command can be modified to accommodate this scenario too).
Why is it limited, you say? Part of the issue is that the `Read Groups` information is wrong. That means that Tapir will naively merge the bams, but duplicates between but not within your two hair samples cannot be detected. That's okay in the hair case above because they involve separate library preps, so there are no (real) cases of this happening.
In truth, when we merge BAMs that have different sample ids, we're creating a multi-sample BAM file. So if we do traditional genotyping on that, you'd get a VCF file for two individuals (with BCFtools).
As it stands, GLIMPSE assumes that we're working with a single-sample file, which is convenient because that's what we want to treat the data as.
A more correct solution is to rewrite the read groups to have the "correct" sample identifier (hair1hair2?), but the result is the same either way.
<br>
We provide a python script to help with this process. 
It also supports using BAMs produced by other pipelines, but you do so at your own risk (ie, it's not Tapir's fault if you get a bad result; that's on you).
Let's say we have two BAM files (`file1.bam` and `file2.bam`) and we want to run step2 of Tapir.
Simply:
```
python3  $TAPIR/bin/makeQuasiReplicates.py -s Merged -i file1.bam file2.bam
```
The above will make Tapir-compliant BAM files (in name only, as symbolic links, but in a directory structure consistent with what Tapir expects). 
By default, it makes a sample name "Combined"; in the above we went with the catchy name "Merged"
Afterwards, simply run Step 2. E.g., 
```
cd Merged
snakemake -c 128 -s $TAPIR/snakemakes/bams2genotypes.smk  --until call_glimpse2"
```

this will run Glimpse on the merge/remark-duplicate result. To reiterate (as described above), if you are not careful with your read-groups, the BCFtools genotyping will be wrong; stick with GLIMPSE.


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
(in context, this section of the config describes the filtering performed on `bcf223andme.py`; in words "bcf to 23 and me".)

### Local storage
Tapir writes (many) temporary files to local storage devices. In short, if you have an NVMe (or traditional SSD) mounted someplace, update:
<br>
`/tmpdirprefix: "/tmp"` <br>
in the config file(s). 
<br>
Better yet, set `/tmp` to the mount location of the NVMe drive...
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
specifically, if you take the same library and sequence it twice, the LB code (part of the read group) should set to be the *same*. If instead you sequence the same sample
using two different library preps, `LB` (within subjects) should be different. Using the `i7_Index_ID` is safe for the first scenario (same library sequenced twice), but it may or may not be true of the second.
That said, the consequences of calling two different libraries the "same" library are pretty minimal; you might artificially call a few more duplicate reads, but the probability of that happening is miniscule, and the result is that your coverage might be a touch smaller than what it would be otherwise.
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
   -  e.g., binary: "/absolute/path/of/myfile.py"
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

## DeepVariant
Tapir also comes prepackaged with DeepVariant (`src/DeepVariant/deepvariant_latest.sif`). You can invoke it by:
`snakemake -c 128 -s src/bams2genotypes.smk --until call_deepvariant`
Tapir provides minimal support DeepVariant; singularity must be installed (there are various tutorials on how to do so; none are provided here). Note that *variant calling* (in short, describe all of the places in the genome that look different from the reference sequence) is not the same thing as *genotyping* (in short, given these sites and alleles, what are the corresponding genotype likelihoods?). For most forensic applications, we want to **genotype** samples; not variant call them.








