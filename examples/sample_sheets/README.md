## SampleSheet Templates
Provided are two sample sheets; each sample sheet is specific to the software used to convert BCLs to FASTQs <br>
[SampleSheet.bcl2fastq.csv](SampleSheet.bcl2fastq.csv)
<br>
and
<br>
[SampleSheet.bcl-convert.csv](SampleSheet.bcl-convert.csv)
<br>
bcl2fastq [External Link](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) is (in the process of being) deprecated; it is also the bcl conversion tool used in the developemental validation of whole genome sequence (WGS) at the Center for Human Identification.<br>
bcl2fastq can be used on the following instruments: MiSeq, HiSeq, iSeq, NovaSeq, and *some* NextSeq instruments
<br><br>
bcl-convert [External Link](https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html) is intended as a replacement for bcl2fastq.<br>
bcl-convert can be used on the following instruments: MiSeq, HiSeq 4000, iSeq, NovaSeq 6000, and NextSeq 500 instruments

### Rules on sample names

Illumina disallows the following sample names: <br>
all, undetermined <br>
<br>
In addition, Tapir treats sample names that start with UD (eg, UDI, UDP) as adapter names (without a sample); eg, these "sample names" are dedicated to adapter sequences that are the result of carryover.
<br>
**In other words, a REAL sample name cannot start with UD**


### Additional Notes
- Both sample sheets are set up to extract multiple samples. At the CHI we extract out all adapter sequences, including those not intended/expected, as those data are important for detecting carryover and/or contamination.
- bcl-convert places many extraction parameters into the sample sheet itself; bcl2fastq instead specifies these arguments on the command line. By default, we attempt to extract out as many fastq records as we can (even if they correspond to very short fragments). These records may or may not persist in the downstream files (e.g., BWA considers a minimum read length).
  - bcl-convert is run with the following flags: --sample-name-column-enabled true --create-fastq-for-index-reads true --fastq-gzip-compression-level 9
    - note that mask-short-adapter-reads and  --minimum-trimmed-read-length flag-equivalents in bcl-convert are specified in the sample sheet.
  -- bcl2fastq is run with the folliwing flags: --mask-short-adapter-reads 0   --minimum-trimmed-read-length 0 --ignore-missing-bcls --fastq-compression-level 9

### Note on binaries
Tapir provides compiled versions of both bcl2fastq and bcl-convert. bcl2fastq is statically compiled and short port to any x86 (64-bit) linux distribution. bcl-convert is dynamically linked. I've tried to fix that, but I ~~can't~~ haven't found a way around yet. Use the `ldd` command to verify that the necessary libraries are found; alternatively, you can download bcl-convert from Illumina and install it yourself.

