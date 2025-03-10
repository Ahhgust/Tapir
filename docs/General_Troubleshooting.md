# General troubleshooting
Tapir didn't do what you wanted it to do. Let's see what we can do about that.

## SampleSheet problems
Sample sheets are by *far* the easiest place to make a mistake.
<br>
A dry-run (`-n` when you run snakemake) will (amongst other things) parse the sample sheet. If you made a mistake, now is your time to fix it.
Common mistakes:
- bad sample names
  -  Illumina forces us to use numbers, letters and - in sample names. **that's it**
     - so **no** periods (0.5ng)
	 - underscores (foo_bar) 
	 - or spaces (sample 5)
  -  Sample names also cannot start with UD
     - This prefix is reserved for indexes.
- bad file formats
  -  Tapir (and all tools downstream) need the sample sheet to be in plain text.
  -  I recommend using [Notepad++](https://notepad-plus-plus.org/) (Windows users)
  -  If you are exporting data from Excel, you want to Export as CSV MS-DOS (**NO UTF-8**)
     -  You can convert utf-8 to plain text using:
	    - ```iconv -f utf-8 -t ascii//TRANSLIT BadSamplesheet.csv > GoodSamplesheet.csv```

<br>
In addition, remember you can specify a modified SampleSheet when you run `bcl2bam`
add (`--config Samplesheet=mySampleSheet.csv`) to your command line.

## Restarting a run and/or incomplete runs
If a run gets killed (did you accidentally close Virtual Box?), you can have snakemake pick up where it left off with:
```
snakemake --rerun-incomplete --keep-going (...)
```
where `(...)` is the rest of your command.
```keep-going``` will tolerate program failures. e.g., BCFtools call can (on occassion) just crash. This is rare, and usually occurs when you have very (very) little data.
Using ```keep-going``` will let the rest of the workflow complete, whereas the default is for snakemake to stop working when any component of it crashes.

## Logs
Whenever a different tool is run, Tapir keeps the output of that tool in a log.
<br>
Some tools are really chatty, however, so don't confuse a big log with a bad result. <br>
If you think something broke (did you need to say ```--keep-going```?), you should investigate why.
<br>
In general, it pays to look at the beginning of the log (head -n 20 foo.log) and the end of the log (tail -n 20 foo.log).
Errors (most often) are reported in these places.
### Run-level logs
bcl-convert and bcl2fastq are run-level tools. A file called `Fastqs/bcl2fastq.outerr` (or bcl-convert.outerr) has all of the chatter
### Per sample (per run) logs
Each sample should have a bam file
and `VCFs`, `Uploads` and a `Logs` directory. They should also have 1 or more Run directories (ie, directories with BAMs in them)
<br>
The log names should be self-explanatory. (this is not meant in a condescending way; we want the user to be able to look at a file, and have a solid idea what the file *is*, based on its name alone)
<br>E.g.:
```
S002-4cm-1-05x.la.md.bqsr.merged.log
```
has the:
- sample name (S002-4cm-1-05x)
  - (it is up to the user to make a good sample name; i.e., one that describes what the sample *is*)
- and the bam file was:
  -  left aligned (la)
  -  duplicates were marked (md)
  -  bqsr was run (bqsr)
-  and the log corresponds to the merging (also, left-align indels) step in Tapir.

## Poor demultiplexing performance
bcl-convert and bcl2fastq use heuristics to demultiplex your samples; by default, the sequences must match (tolerating up to 1 mismatch). <br>
We also provide a binary for [deML](https://github.com/grenaud/deML), which demultiplexes by maximum likelihood. 
If your run had very poor results for the sample indexes (_I1_ and the _I2_ fastq files), then your demultiplexing will be poor.
`deML` can rescue some of these reads (and in fact, it's the reason we extract the indexes in fastq format). Adding deML support is something we hope to support in the future.
However, you (the user) can manually run it....
