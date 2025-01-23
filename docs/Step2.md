![Step2](../images/Stage2.jpg)

Items in teal are for *user-interpretation*. Items in dark blue are programs and tools. The final result (bam) is depicted in green.

# Step 2
## Quick references
-  [Click here](Glossary.md) for a list of tools, terms and definitions.

## Before you start
In words:
- In **Step 1**, we took our raw BCLs 
   -  n>=1 samples
- and made *n* + 2 (+ Offtargets and Undetermined) directories
   -  with processed bam files
      -  .la.md.bqsr.bam
	     - left-align (indel realignment)
		 - mark-duplicates (pcr/optical duplicates marked)
		 - bqsr (base-quality scores have been empirically recalibrated)
-  Internally, we found that:
   -  at least 0.10× Coverage is needed
      - for Glimpse (imputation)
   -  and 5-10× Coverage is needed
      -  for BCFtools (maximum likelihood; no imputation)
-   Likewise, mixture fractions <10% (idealing <1%) are highly recommended
   
In **Step 2** you select *one* of the *n* samples		 
-  You've completed **Step 1**
   - Checked the [QC](QC.md)
     - of each BAM
   -  And it passes *your labs* quality metrics

## The flow
Within a biological sample (with *m* libraries)
-  Composed of *m* measurments (bams)
   -  bams are merged (GATK)
      -  per-bam read-group information ensures a correct merger. [Click here](PowerUsers.md) for an explanation.
   -  duplicates are re-marked
      -  in case the same library was sequenced (at least) twice.
   -  Tapir optimizes the common case (m=1)
      - symbolic links when only 1 bam file is available; no merging, no PCR duplicate remarking.
-  Autosomal genotyping is performed using
   -  BCFtools OR
   -  GLIMPSE
      -  Which tool depends on the user
-  The ploidy of the X chromosome (X vs XX hypothesis) is estimated
   -  Ploidy-aware genotyping is performed on the X
   -  And the results are merged with the autosomal data
-  The resulting VCF file 
   - hg38 reference genome
-  is converted to hg19 
   - LiftoverVCF
-  and a "23andme" file is made.
   -  uses Tapir's `bcf223andme.py` script
   -  applies threshold embedded in the `config.yaml` file
      - By default: 
	     - BCFtools: Genotype quality >20
		 - Glimpse: Genotype posterior probability >0.99 AND a Bayes Factor > 1.7