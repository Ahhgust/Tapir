# Changelog

### Beta changes (post-publication)

### Alpha changes (pre-publication)
4/2025
-  Bugfix with GATK's BQSR
	-  In some verybizarre cases, BQSR produces a corrupt BAM file.
		-  In the one and only case, the BAM file had two unmapped read pairs (one of each pair was an empty fastq record)
	-  When BQSR fails, Tapir now produces an empty BAM file (based on the header to the pre-bqsr file)
	   -  No BQSR plots are produced in this case, however.
	-  The final rule (repro) no longer requires the BQSR *plots*
	   -  Though it does require BQSR bams 
	   -  This affects "snakemake --keep-going", which does not work as expected (IMO, a bug in snakemake)

	
3/2025
-  Support for multiple config files added.
	-  at most 1 config can be used at any 1 time.
	-  config files and command line arguments are copied (bcl2bam and bams2genotypes)
	-  for the latter, config files copied per genotyper.
-  Major refactor of directory structures.
   -  all unix routines are unchanged.
   -  refactor also retools fastq input.
      - bugfix 3/11 and 3/24
-  Tape backup routines added (see: tape_deck/)
-  .la.md bam files are set as temporary files
   -  use --notemp (snakemake directive) to keep all temporary files instead
	
	
2/2025
-  Added support for FASTQ input (as well as BCL)
-  Dexixtify now run iff multiple bam files are merged.
-  bcl2bam detects non-ascii text samplesheets (and throws an error)
	
11/2024
-  Added X chromosome support
   -  Estimates X ploidy (1 or 2) by maximum likelihoood
   -  Uses this information to genotype (haploid or diploid)
-  BCFtools call
   -  turned off the "prior" (-P 0.0)
   -  leads substantial improvements in lower-coverage genomes
   -  performance is still (much) worse than GLIMPSE (when below ~5x)


	
