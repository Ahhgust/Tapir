# Resources
Tapir **uses** quite a few genomic resources. Here is a brief accounting of them. <br>
Specifically, here's an accounting of the files in the `resources/` directory in Tapir. <br>
See [File.md](Files.md) for an accounting of the files *made* by Tapir.
<br>
Unless stated otherwise, all coordinates are relative to the hg38 reference genome. <br>
Woerner et al. 2022 refers to: [External Link](https://doi.org/10.1016/j.fsigen.2022.102785) <br>
Woerner et al. 2025 refers to: (The tapir manuscript; currently under review. TODO: update when published) <br>
Koenig et al. refers to: [External Link](https://doi.org/10.1101/gr.278378.123) <br>
Poznik et al. [External Link](https://doi.org/10.1126/science.1237619) <br>

-  7_matching_X.txt
    - See `X_nonpar.txt`
    -  plain text file containing: "chr7:2781480-155701382"
	-  Matching genomic coordinates to the NON-pseudoautosomal region on the X, but on chromosome 7.
	-  File used to estimate the ploidy of the X chromosome.
-  atlas.excludesites.*
	-  A "blacklist" of sites; composed of the mask of Woerner et al 2022. AND snvs from Koenig et al.
	-  See: `Readme.txt`
    -  various extensions include compressed files (.gz), and "zero-based coordinate" beds (atlas needs these, though that's not bed format!)
    - the "atlas.include\*" files are the set complement of the exclude files. 
-  bcftools_call/
	-  SNVs (sites + alleles) used for genotyping with bcftools; tsvs, in `bcftools call` format.
    -  See Woerner et al. 2025 for extended details
	-  In brief, HGDP+1000 Genomes samples from Koenig et al. , relatives removed, sites with a PP (phasing probability) annotation removed.
	-  For the X see the chrX/ documentation
-  chrX/
    -  Information necessary to genotype the X chromosome.
	-  The glimpse2 directory also includes information on making the "chunks" in GLIMPSE.
    -  see `chrX/xsnppanel.readme`
-  gsa_24v3-0_A2.GRCh38.rand10k.bed
	-  10,000 autosomal sites selected from Illumina's Global Screening Array (24v3-0_A2)
	-  masking out positions difficult to assay by WGS (see Woerner et al 2022)
-  GSA-24v3-0_A2.hg38.gnomadannos.autos.sites2include.justafs.*
    -  586,624 autosomal sites from Illumina's Global Screening Array (24v3-0_A2)
	-  masking out positions difficult to assay by WGS (see see Woerner et al 2022)
	-  including allele frequency annotations (AF tags from `bcftools annotate`) from gnomAD (v3.1.2)
-  hg19/
   -  The hg19 reference genome as downloaded from UCSC.
   -  Includes various indexes (BWA, and a .dict file)
      -  See README.txt
   -  The GATK_Resources/ data are not used.
   -  chrAll.standardChroms files (canonical chromosomes; chr1-22, X,Y,MT) are are used by Tapir.
-  hg38/
   -  The hg38 reference genome.
   -  The same version used by the *1000 Genomes Project* [External ftp link](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)
   -  Includes various indexes (BWA, and a .dict file)   
-  hg38ToHg19.over.chain
   -  A "liftover" file, used to convert between hg38 to hg19. Developed by UCSC.
   -  See [External Link](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/)
-   hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp*
   -  Autosomal SNV panel used by Tapir, with various variations.
   -  `hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.1per.noindels.hg19.filt.vcf.gz`
      -  SNPs (1%MAF), no indels, hg19 coordinates.
	  -  `filt`ered to remove sites where the chromosomes (hg38 vs hg19) disagree.
   -  `hgdp1kgp_autos.filtered.SNV_INDEL.phased.shapeit5.unrelateds.nopp.merged.2kbslop.bed`
      -  All autosomal markers merged into discrete non-overlapping regions at least 2kb in size.
	     -  bedtools slop 2kb; and merged (set union)
		 -  See: notes.slop
	  -  Not used as of yet (may be helpful for Atlas)
-  ibis.sexaveraged\*gmap
    -  Not currently used by Tapir.
    -  Compatible with IBIS. [External Link](https://github.com/williamslab/ibis)
	-  hg38 (eg,  chr1) and GRCH38 (eg, 1) genetic maps.
	-  X chromosome rates are population rates; not the rates in females.
	   -  i.e., 2/3 of the female rate.
-  plinkgsa\/
    -  Not currently used by Tapir.
	-  Genotypes from Koenig et al., GSA Sites, Plink format.
-  Poznik.YChromCallable.hg38.bed
    -  Not currently used by Tapir.
    -  Genomic positions from Poznik et al. 
	   -  Sites argued to be suitable for variant calling on the Y chromosome.
	      - Their "callability" mask.
	   -  May be used soon for X/Y ploidy estimation (sex/mixture assessments)
	   -  Genomic coordinates were converted to hg38 using liftOver (command line tool; defaults)
-  snppanel_autos_x.vcf.gz
   -  True "SNPS" used by TAPIR. (SNVS with MAF>1%)
   -  hg19 reference genome.
   -  Includes sites where the chromosomes do not match *hg19 vs hg38*
-  split/
	-  See: `Readme.txt` 
	-  The reference panel used by GLIMPSE2; binary format
-  X_nonpar.txt
    -  See `7_matching_X.txt`
    -  plain text file containing: "chrX:2781480-155701382"
	-  The canonical coordinates of the non-pseudoautosomal region of chromosome X