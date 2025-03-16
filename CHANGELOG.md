# Changelog


### Beta
3/2025
-  Support for multiple config files added.
   -  at most 1 config can be used at any 1 time.
-  Major refactor of directory structures.
   -  all unix routines are unchanged.
   -  refactor also retools fastq input.
      - bugfix 3/11
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


	
