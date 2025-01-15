![Step1](../images/Stage1.jpg)

Items in teal are for *user-interpretation*. Items in dark blue are programs and tools. The final result (bam) is depicted in green.

# Step 1
## Quick references
-  ![Click here](Glossary.md) for a list of tools, terms and definitions.

## The flow
In words:
-  Raw BCL data (Illumina sequencer) are demultiplex (converted to fastq) by
   -  bcl2fastq (MiSeq, NovaSeq, HiSeq...) **OR**
   -  bcl-convert (NextSeq, ...)
      -  Which is chosen depends on the format of the SampleSheet.
-  In either case, Tapir will demultiple
   -  Samples (per lane) **AND**
   -  Indexes (ie, indexes without an accompanying sample; useful for detecting carryover)
-  All fastqs are aligned using BWA mem.
   - One bam per lane per sample/index
-  A per-sample BAM is created
   -  indels are realigned (gatk)
   -  Input:  1+ bams
   -  Output: 1 bam
-  Duplicates are marked (sambamba) and the following summaries are made:
   -  Flagstat (per-read summary statistics)
   -  Samstats (per-site summary statistics)
   -  Fastqc is run
      -  And its report is created
-  Base quality scores are empirically recalibrated (BQSR)
   -  The results are tabulated/plottted (plot_bqsr)
   -  And the same is tested to see if it is a mixture (Demixtify)
      -  Which generates both a detailed report (.demix) AND
      -  a summary
-  The final bam has been left-aligned around indels (.la), duplicates have been marked (.md), and BQSR has been applied (.bqsr).

## How to
Tapir is a snakemake workflow; snakemake lets us (easily) stitch all of these tools together, as well, snakemake keeps track of program "state". eg, if the power is turned off, snakemake can deduce what it was doing, and restart computation at the appropriate location.
<br>
An Illumina instrument will (by default) write to a directory called "Output". <br>
First, `cd` that directory. On our system, that amounts to:
```
cd /eva/staging/Novaseq/Novaseq/Output/
```
Suppose a sequencing run just finished; let's call it: `250103_A01324_0120_BHWGY2DMXY`

First, attempt a dry-run (-n)
```
snakemake -n  -s /eva/codebase/snakemakes/bcl2bam.smk  -c256 --config Bcldir=250103_A01324_0120_BHWGY2DMXY Outdir=/eva/datums/cooked/NovaSeq001/2025/WGSValidation Experiment=Sensitivity
echo $?
```

Which tells snakemake to evaluate all of the information, and evaluate what needs to be run. Dry-runs are "chatty", so it can be a bit hard to tell if errors or present (or not, depending on the error). The error code is also printed `echo $?`; the last line *should* say 0. If not, you may have mis-specified something.
Parameters (ie, things that start with a "-") of note:
-  \-s
   -  Which snakemake script are we calling?
-  \-c
   -  How many cores/cpus do you want to allocate?
      - Adjust as necessary
   -  Bcldir=
      - Set this to the BCL directory you want to extract
         -  In the case of the NovaSeq, a "CopyComplete.txt" file is expected
   -  Outdir=
      -  Set this to the (base) directory;
   -  Experiment=
         -  This *can* be set in the SampleSheet (Experiment= ...)
         -  Set this to the (parent) directory.
         -  In the above, results are written to:
	 -  `/eva/datums/cooked/NovaSeq001/2025/WGSValidation/Sensitivity/250103_A01324_0120_BHWGY2DMXY`
   - Optional 
      - Samplesheet=
         -  Provide a path (and filename) for a different sample sheet.
	    -  This can be useful if you messed up the sample sheet the first time.
	 -  By default, SampleSheet.csv (in the BCL dir) is used.


The dry-run will also parse the provided sample sheet; there are a lot of ways to mess up a sample sheet! Please see the writeup:
<br>
![here](../examples/sample_sheets/README.md) for a reminder of how what the sample sheet *should* look like.
<br>

If your dry-run is successful, let's run Tapir for real. I recommend:

```
nohup snakemake -s /eva/codebase/snakemakes/bcl2bam.smk  -c256 --config Bcldir=250103_A01324_0120_BHWGY2DMXY Outdir=/eva/datums/cooked/NovaSeq001/2025/WGSValidation Experiment=Sensitivity &
```

Which runs the command in the background (trailing &) and lets you exit/logout of the shell without aborting the call (nohup; no hardware interrupt)

<br>
Note that when Step 1 completes, that would be a good time to evaluate both sample-level QC metrics; see ![here](./QC.md) <br>
As well, run-level metrics are of key importance.

### Files made

Step1 of Tapir takes a single run, and creates demultiplexed results; one for each sample. Additionally, "Offtargets" data are created; these correspond to fastqs records that are supported by your library prep method, but you *shouldn't* see any results for. Because demultiplexing is imperfect (really, we observe reads with error), some (tiny amount) of data may be present, anyways.

/eva/datums/cooked/NovaSeq001/2025/WGSValidation

## TMI
Snakemake workflows are best implemented without arguments (ie, things that change how the programs are run). Tapir does not entirely adhere to this philosphy. Specifically, in Step 1, a traditional Snakemake workflow would write to the BCL directory (to preserve state). That is a bad idea; the original data should be treated as a WORM (write once, read many times); that way it can be easily backed up as soon as it comes off the instrument (and MD5s and whatnot will match).

