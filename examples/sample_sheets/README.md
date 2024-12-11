## SampleSheet Templates
Provided are two sample sheets; each sample sheet is specific to the software used to convert BCLs to FASTQs <br>
[SampleSheet.bcl2fastq.csv](SampleSheet.bcl2fastq.csv)
<br>
and
<br>
[SampleSheet.bcl-convert.csv](SampleSheet.bcl-convert.csv)
<br>
[bcl2fastq](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) (External Link) is (in the process of being) deprecated; it is also the bcl conversion tool used in the developemental validation of whole genome sequence (WGS) at the Center for Human Identification.<br>
bcl2fastq can be used on the following instruments: MiSeq, HiSeq, iSeq, NovaSeq, and *some* NextSeq instruments
<br><br>
[bcl-convert](https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html) (External Link) is intended as a replacement for bcl2fastq.<br>
bcl-convert can be used on the following instruments: MiSeq, HiSeq 4000, iSeq, NovaSeq 6000, and NextSeq 500 instruments

<br><br>
Tapir digests the sample sheet, and uses the bcl converter that matches the format of the sample sheet. IE, you tell Tapir which extraction method you want by selecting the appropriate sample sheet format. If you wish to deviate from the templates provided, please note the flags/options (at the bottom of this document).


### Rules on sample names

Illumina disallows the following sample names (presumably case insensitive): <br>
all, undetermined <br>
<br>
In addition, Tapir treats sample names that start with UD (eg, UDI, UDP) as adapter names (without a corresponding sample); eg, these are UDIs that are used in the lab (generally), but are not intended for this run.
<br>
**In other words, a REAL sample name cannot start with UD**


### Additional Notes
- Both sample sheets are set up to extract multiple samples. At the CHI we extract out all adapter sequences, including those not intended/expected, as those data are important for detecting carryover and/or contamination.
- bcl-convert places many extraction parameters into the sample sheet itself; bcl2fastq instead specifies these arguments on the command line. By default, we attempt to extract out as many fastq records as we can (even if they correspond to very short fragments). These records may or may not persist in the downstream files (e.g., BWA considers a minimum read length).
  - bcl-convert is run with the following flags: --sample-name-column-enabled true --create-fastq-for-index-reads true --fastq-gzip-compression-level 9
    - note that mask-short-adapter-reads and  --minimum-trimmed-read-length flag-equivalents in bcl-convert are specified in the sample sheet.
  -- bcl2fastq is run with the folliwing flags: --mask-short-adapter-reads 0   --minimum-trimmed-read-length 0 --ignore-missing-bcls --fastq-compression-level 9

### Note on binaries
Tapir provides compiled versions of both bcl2fastq and bcl-convert. bcl2fastq is statically compiled and should port to any x86 (64-bit) linux distribution. bcl-convert is dynamically linked. I've tried to fix that, but I ~~can't~~ haven't found a way around yet. Use the `ldd` command to verify that the necessary libraries are found (e.g., `ldd bin/bcl-convert` and check for the keyword "not found"; if something is "not found" you need to fix that, either by updating your LD_LIBRARY_PATH environmental variable OR by placing symlinks, OR by updating/installing various packages.); alternatively, you can download bcl-convert from Illumina and install it yourself.

