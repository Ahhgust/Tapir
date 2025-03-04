

### Where the data came from

Files were drawn from the Gambian Genome Variation Project [GGVP](https://www.internationalgenome.org/gambian-genome-variation-project/) <br>
A single individual (female) from the Wolof population was selected (SC_GMWOF5428715) and downsampled to a 1x genome.
The cram file was exported to a pair of fastq files using:

```
samtools sort -u -m 50G -n SC_GMWOF5428715.alt_bwamem_GRCh38DH.20151208.WOLOFF.gambian_lowcov.1.0x.cram | samtools view -s 1.5 | samtools fastq -1 GMWOF5428715.r1.fastq.gz -2 GMWOF5428715.r2.fastq.gz --reference /eva/edatums/reference_materials/reference_genomes/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

Just to make life a bit easier (For me), I went ahead and reduced the coverage by a factor of 2 (~1x -> 0.5x; samtools view -s 1.5 does this). <br><br>

Tapir requires that the reads come from a (relatively modern) illumina instrument. See the Wikipedia article [here](https://en.wikipedia.org/wiki/FASTQ_format); specifically, the Illumina sequence identifiers.
The GGVP data contain non-standard fastq identifiers. To make the Fastq records compatible with Illumina, let's go ahead and re-make the fastq ids:

```
zcat GMWOF5428715.r2.fastq.gz | perl -e 'srand(1); while (<>) { if ($. % 4 == 1) { chomp; $o='@A01324:76:HFGHFDRX3:1:' . int(rand(1000)) . ":" . int(rand(10000)) . ":" . int(rand(10000)); print $o , "\n";} else { print $_; }}' | gzip -9 > GMWOF5428715_S89_L001_R2_001.fastq.gz
zcat GMWOF5428715.r1.fastq.gz | perl -e 'srand(1); while (<>) { if ($. % 4 == 1) { chomp; $o='@A01324:76:HFGHFDRX3:1:' . int(rand(1000)) . ":" . int(rand(10000)) . ":" . int(rand(10000)); print $o , "\n";} else { print $_; }}' | gzip -9 > GMWOF5428715_S89_L001_R1_001.fastq.gz
```

Of note, "making up" the IDs may impact the rate of optical duplicates. As this is an example, I am less worried about that; if this were real data, I would be more concerned.
Also note, the format of the FASTQ file; this is the format that Illumina uses. It is: SampleName_SampleID_Lane_R1/2_001.fastq.gz



