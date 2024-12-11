## SampleSheet Templates
Provided are two sample sheets; each sample sheet is specific to the software used to convert BCLs to FASTQs <br>
[SampleSheet.bcl2fastq.csv](SampleSheet.bcl2fastq.csv)
<br>
and
<br>
[SampleSheet.bcl-convert.csv](SampleSheet.bcl-convert.csv)
<br>
bcl2fastq [External Link](https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) is (in the process of being) deprecated; it is also the bcl conversion tool used in the developemental validation of whole genome sequence (WGS) at the Center for Human Identification.
bcl2fastq can be used on the following instruments: MiSeq, HiSeq, iSeq, NovaSeq, and *some* NextSeq instruments
<br><br>
bcl-convert [External Link](https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html) is intended as a replacement for bcl2fastq. 
bcl-convert can be used on the following instruments: MiSeq, HiSeq 4000, iSeq, NovaSeq 6000, and NextSeq 500 instruments

### Additional Notes
- Both sample sheets are set up to extract multiple samples. At the CHI we extract out all adapter sequences, including those not intended/expected, as those data are important for detecting carryover and/or contamination.
- bcl-convert places many extraction parameters into the sample sheet itself; bcl2fastq instead specifies these arguments on the command line. By default, we attempt to extract out as many fastq records as we can (even if they correspond to very short fragments). These records may or may not persist in the downstream files (e.g., BWA considers a minimum read length).


