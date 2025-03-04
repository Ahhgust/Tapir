# Changelog


### Beta
3/2025
-  Support for multiple config files added.

2/2025
-  Added support for FASTQ input (as well as BCL)
-  Mixtify now run iff multiple bam files are merged.
-  bcl2bam detects non-ascii text samplesheets (and throws an error)
	
11/2024
-  Added X chromosome support
   -  Estimates X ploidy (1 or 2) by maximum likelihoood
   -  Uses this information to genotype (haploid or diploid)
-  BCFtools call
   - turned off the "prior" (-P 0.0)
   - leads substantial improvements in lower-coverage genomes


	
