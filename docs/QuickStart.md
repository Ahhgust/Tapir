![Tapir](../images/Tapir.png)

# Tapir
## Quick start
## Step 0:Install Tapir
After Tapir has been [installed](Install.md), you (or, if you're careful about things, any/select other users on the unix machine) can take data from any Illumina instrument and convert these data into a pre-processed BAM file (step 1); a BAM is created for each sample in the sample sheet (optionally, for each index as well). If the sample passes QC, step2 in Tapir calls genotypes and produces a GEDmatch-compatible file.
<br>
As a sanity test, if you:
```
echo $TAPIR
```
it should give you the path of where TAPIR is installed. Anything less (e.g., an empty line), means you haven't finished installing Tapir.
## Step 1: Sample sheets
- Pick a sample sheet format
  - bcl-convert or
  - bcl2fastq
    -- which format you choose will tell Tapir which program to run.
- Add your run information to the samplesheet.
  - *This is the most common source of error!*
  - Do
    - Set the `Experiment Name`
	- Set the sample names
	  - numbers of letters (or -) **only**
	  - sample identifiers CANNOT start with UD
	    - (this is how bioinformatic negatives are marked)
	  - do make sure you have an identifier for each *potential* UDI/UDP that you use.
	  
## Dry run


By example, load Tapir
```
mamba activate tapir
```
(replace `mamba` with `conda` if you opted for a conda installation) <br>
Try a dry run:

```
snakemake -n  -s $TAPIR/snakemakes/bcl2bam.smk  -c16 --config Bcldir=/eva/staging/Novaseq/Novaseq/Output/250103_A01324_0120_BHWGY2DMXY Samplesheet=./MySampleSheet.csv Experiment=WriteDataHere --configfile $TAPIR/configs/config_v_2_low_mem.yaml
```

In words: 
- do a dry run (-n)
  - This means, just evaluate the information and tell me what things will be run.
- using this snakemake file (bcl2bam)
- Use 16 cores
  - If you have more than 16 available/free, set this to a higher value!
- Override the default configuration with:
  - Extract this BCL directory ( 250103_A01324_0120_BHWGY2DMXY )
  - You *want* to override the BCL directory.
- The following overrides are *optional*
  - Write the fastqs to `WriteDataHere/`
     - doing this will OVERWRITE the default `Experiment Name` used in the sample sheet
  -  You can also add (after the --config)
     - The path to a *different* sample sheet using `Samplesheet=ModifiedSampleSheet.csv`
     - (otherwise, the sample sheet in the BCL directory is used.
  -  And you can use an *entirely different* configuration file using `--configfile`	 
<br>
This will make a bunch of output, but it should tell you that bcl-convert of bcl2fastq will be run. Which downstream programs are run will depend on the output of these tools.

If you are happy with the output, try running Tapir for real:

## Live run
Take whatever variant of the above command you used, and just remove the `-n`. Running in the background is probably helpful too:
For example (modify the below to match the flags/additions you made from the above)
```
nohup snakemake -s -s $TAPIR/snakemakes/bcl2bam.smk  -c16 --config Bcldir=/eva/staging/Novaseq/Novaseq/Output/250103_A01324_0120_BHWGY2DMXY   &> tapirOutput.oe &
```	
which removes the -n (not a dry run, but the real thing), and tells unix to write all of the output to `tapirOutput.oe`, to run this command in the background (final `&`), and even if I log out, let this program run (`nohup`).
<br>
If Tapir finishes quickly, check `tapirOutput.oe`. Odds are you specified some directory that you don't have permission to read/write to.
Other than that, Tapir should run (perhaps for the next day or 3), depending on what it is you're doing.

## Step2
QC your run; Do you FASTQs make sense? Are there a lot of reads in a UDI read?
It's very hard for your samplesheet to be right; most people do something wrong. If the results look suspect, you can delete the directories and run the workflow again.
If you're happy with the results, do step 2. Step 2 is run once per sample.
so `cd` into `Sample_Data` directory. 