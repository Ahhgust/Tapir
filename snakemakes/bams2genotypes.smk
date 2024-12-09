import glob, os, sys, re, time, gzip
import GenomixHelper

configfile: os.path.join(workflow.current_basedir, "config_v2_1.yaml")

GX= GenomixHelper.GenomixHelper(workflow.current_basedir, config)

### Let's grab some resources/binaries that we'll need a lot
#
# All GX.get methods return an absolute path
# paths are 

# reference genome
hg38=GX.getResources("Reference")

# bcftools binary:
bcftools=GX.getParam("bcftools", "binary")

sampleName=config["Samplename"]
# non-empty suffix; can either be 'la.md.bqsr'
# (left-align, mark duplicates, base quality score recalibration; aka, bqsr)
# OR
# la.md
# (Not using BQSR doesn't make a lot of sense though)
suffix=config["Suffix"]

# be default, it assumes you're in the directory of some sample within some experiment
# in which case, the samplename is the name of the directory
if sampleName=="":
    sampleName=os.path.basename( os.getcwd() )

# left-aligned, mark duplicates, then base quality score recalibration...
suffix="la.md.bqsr"

  
bams=glob.glob(f"*/Bams/*{suffix}.bam")
if len(bams)==0:
    print("No bam files detected?!", "Are you sure you're in right directory?", sep="\n", file=sys.stderr)
    exit(1)
#TODO: refine the BAM list using tags/decorators
#TODO: Check for {outdir}/{rundir}/RunAnalysisComplete.txt"
    

def bamGenerator(wildcards):
    """
    Convenience function; lists all of the final merged/markdup bams needed
    """
    return expand("{samplename}.{suffix}.merged.md.bam", samplename=sampleName, suffix=suffix)



def listBams(wildcards):
    """
    Lists the bams (in directory subtree)
    that are to be used for merging
    TODO: mark/signal BAMs that "fail"
    """
    return bams


module bcl2bam:
    snakefile: 
        "bcl2bam.smk"
    config: config
#TODO:
# split off the summary stats routines into their own snakemake     
# TODO: 
# make summary statistics of the merged bam.
#use rule make_bam_samstats from bcl2bam as bcl2bam_make_bam_samstats


print("Processing: ", sampleName, suffix)

# current_basedir returns the path of this snakemake script
# bin/ and resources/ are in sibling directories; refer to them relative to ROOT
ROOT=os.path.abspath(  os.path.join(workflow.current_basedir, ".."))

print("Path: ", ROOT)


wildcard_constraints:
    suffix=suffix

rule call_bcftools:
    input:
        bamGenerator,
        expand("Uploads/{caller}/{samplename}.{suffix}.{caller}.hg19.23andme.tsv.gz", samplename=sampleName, suffix=suffix, caller=['bcftools'])

rule call_glimpse2:
    input:
        bamGenerator,
        expand("Uploads/{caller}/{samplename}.{suffix}.{caller}.hg19.23andme.tsv.gz", samplename=sampleName, suffix=suffix, caller=['glimpse2'])

rule call_glimpse2_and_bcftools:
    input:
        bamGenerator,
        expand("Uploads/{caller}/{samplename}.{suffix}.{caller}.hg19.23andme.tsv.gz", samplename=sampleName, suffix=suffix, caller=['glimpse2', 'bcftools'])

rule call_deepvariant:
    input:
        bamGenerator,
        expand("VCFs/deepvariant/{samplename}.{suffix}.deepvariant.vcf.gz", samplename=sampleName, suffix=suffix)
        
rule call_all:
    input:
        bamGenerator,
        expand("VCFs/deepvariant/{samplename}.{suffix}.deepvariant.vcf.gz", samplename=sampleName, suffix=suffix),
        expand("Uploads/{caller}/{samplename}.{suffix}.{caller}.hg19.23andme.tsv.gz", samplename=sampleName, suffix=suffix, caller=['glimpse2', 'bcftools'])

# This takes in 1+ bams, all of the same sample (but possibly different libraries)
# and creates a merged bam file of them (note; no index)
# optimized for the common-case; this creates symlinks iff there is 1 bam to "merge"        
rule merge_bams:
    input:
        listBams
    output:
        temp("{samplename}.{suffix}.merged.bam")
    params:
        binary=GX.getBinary("gatk"),
    log:
        "Logs/{samplename}.{suffix}.merged.log"        
    run:
        inputs_all= " ".join([" I=" + repr(i) for i in input])
        if len(input)>1:
            shell( # note: sambamba generates .bais by default
            """
            {params.binary} --java-options {config[gatkParams][javaoptions]} MergeSamFiles CREATE_INDEX=false O={output} {inputs_all} &> {log}
            """)
        else: # common case; a trivial "merge" of a single bam
            shell("ln -s {input} {output} &> {log}")


# takes in a merged file (merge_bams) and marks duplicates
# note that this happens w/in a sample, within a library (LB tag in the bam file)
# optimized for the common case (the merged file is a symlink; aka, it comes from a single sample/library/run, 
# in which case the file has already been marked.
rule remark_duplicates:
    input:
        "{samplename}.{suffix}.merged.bam"
    output:
        bam="{samplename}.{suffix}.merged.md.bam",
        bai="{samplename}.{suffix}.merged.md.bam.bai"
    threads: 
        config["sambambaParams"]["threads"]
    params:
        tmpprefix=config["tmpdirprefix"],
        binary=GX.getBinary("sambamba") #config["sambambaParams"]["binary"]
    log:
        "Logs/{samplename}.{suffix}.markdup.log"   
    run: # todo: check bais..
        if os.path.islink( str({input}.pop()) ): # only true if no data were merged; as such, duplicates have already been marked
            shell("ln -sL `readlink {input}` {output.bam} &> {log} && ln -sL `readlink {input}`.bai {output.bai} &>> {log}")
        else:
            shell("{params.binary} markdup  --tmpdir={params.tmpprefix} -l 9 -t {threads} {input} {output.bam} &> {log}")

# calls genotypes using bcftools.
# only weakly parallelized ; supports per-chromosome level parallelism
rule bcftools_genotype_a:
    input:
        bam="{samplename}.{suffix}.merged.md.bam"
    output:
        "VCFs/bcftools/{samplename}.{suffix}.autos.bcftools.vcf.gz"
    log:
        "Logs/{samplename}.{suffix}.autos.bcftools.log"
    threads: 
        config["bcftoolsParams"]["callthreads"]
    params:
        ref=hg38,
        binary=bcftools,
        panel=GX.getParam("bcftools", "panel"),
        mpileopt=config["bcftoolsParams"]["mpileup_options"]
    shell:
        """
        if ! outdir=`mktemp -d --tmpdir={config[tmpdirprefix]} XXXXXXbcftools_call`; then
            echo "Failed to make tempdir!"
            exit 1
        fi
        for chrom in {config[chroms]}; do
          echo "{params.binary} mpileup {params.mpileopt} --fasta-ref {params.ref} -r chr$chrom {input.bam} |\
            {params.binary} call -Am -P 0. -C alleles -T {params.panel}$chrom.tsv.gz -Ob1 -o $outdir/chr$chrom.bcf"
        done | parallel -j {threads} 2>> {log}
        
        {params.binary} concat -Oz9 -o {output} --threads {threads} $outdir/chr?.bcf $outdir/chr??.bcf 2>> {log}
        {params.binary} index --threads {threads} {output} 2>> {log}
        rm -rf $outdir
        """
        #"""
        #({params.binary} mpileup {params.mpileopt} --fasta-ref {params.ref} -T {params.panel} {input.bam} |\
        #{params.binary} call -Am -P 0. -C alleles -T {params.panel} -Ov9 --threads {threads} -o {output}) 2> {log} &&\
        #{params.binary} index --threads {threads} {output} 2>> {log}
        #"""


rule combine_genotypes_xa:
    input:
        autos="VCFs/{caller}/{samplename}.{suffix}.autos.{caller}.vcf.gz",
        xchrom="VCFs/{caller}/{samplename}.{suffix}.chrX.{caller}.vcf.gz"
    output:
        "VCFs/{caller}/{samplename}.{suffix}.{caller}.vcf.gz"
    log:
        "Logs/{samplename}.{suffix}.{caller}.log"
    threads: 
        config["bcftoolsParams"]["threads"]
    params:
        binary=bcftools
    shell:
        """
        {params.binary} concat -Ov9 -o {output} {input.autos} {input.xchrom} 2> {log}
        {params.binary} index --threads {threads} {output} 2>> {log}
        """

rule glimpse2:
    input: 
        bam="{samplename}.{suffix}.merged.md.bam",
    output:
        "VCFs/glimpse2/{samplename}.{suffix}.autos.glimpse2.vcf.gz"
    log:
        phase="Logs/{samplename}.{suffix}.autos.glimpsephase.log",
        ligate="Logs/{samplename}.{suffix}.autos.glimpseligate.log"
    threads:
        config["glimpse2Params"]["threads"]
    params:
        phasebinary= GX.getParam("glimpse2Params","phase_binary"),
        phaseopt=config["glimpse2Params"]["phase_options"],
        ligatebinary= GX.getParam("glimpse2Params","ligate_binary"),
        panelprefix=GX.getParam("glimpse2Params", "reference_database_prefix"),
        bcftoolsbinary= bcftools        
    shell: 
        """
        if ! outdir=`mktemp -d --tmpdir={config[tmpdirprefix]} XXXXXXglimpse2`; then
            echo "Failed to make tempdir!"
            exit 1
        fi
        
        for reference in {params.panelprefix}.chr*bin; do
            bn=`basename $reference`
            chrom=`echo $bn | cut -f2 -d '_'`
            start=`echo $bn | cut -f3 -d '_'`
            stop=`echo $bn | cut -f4 -d '_'`
            echo "{params.phasebinary} {params.phaseopt} --bam-file {input.bam} --reference $reference --output $outdir/$chrom.$start.$stop.bcf "
        done | parallel -j {threads} &> {log.phase}
        mkdir -p $outdir/ligate
        # note that curly braces are a total pain. hence the 1..22 (explicitly) below
        for chrom in {config[chroms]}; do
            echo "{params.ligatebinary} --input <(ls -1v $outdir/chr$chrom.*bcf) --output $outdir/ligate/chr$chrom.bcf"
        done | parallel -j {threads} &> {log.ligate}
        # and concatenate to make an autosomal vcf
        {params.bcftoolsbinary} concat -Oz9 --threads {threads} -o {output} $outdir/ligate/chr?.bcf $outdir/ligate/chr??.bcf &>> {log.ligate}
        {params.bcftoolsbinary} index --threads {config[bcftoolsParams][threads]} {output} &>> {log.ligate}
        rm -rf $outdir  
        """

rule glimpse2_x:
    input: 
        bam="{samplename}.{suffix}.merged.md.bam",
        xploidy="VCFs/{samplename}.{suffix}.xploidy.txt"
    output:
        vcf="VCFs/glimpse2/{samplename}.{suffix}.chrX.glimpse2.vcf.gz",
        samplefile="VCFs/glimpse2/{samplename}.{suffix}.chrX.samplesfile.tsv"
    log:
        phase="Logs/{samplename}.{suffix}.chrX.glimpsephase.log",
        ligate="Logs/{samplename}.{suffix}.chrX.glimpseligate.log"
    threads:
        config["glimpse2Params"]["threads"]
    params:
        phasebinary= GX.getParam("glimpse2Params","phase_binary"),
        phaseopt=config["glimpse2Params"]["phase_options"],
        ligatebinary= GX.getParam("glimpse2Params","ligate_binary"),
        panelprefix=GX.getParam("glimpse2Params", "reference_database_prefix_chrX"),
        bcftoolsbinary= bcftools        
    shell: 
        """
        if ! outdir=`mktemp -d --tmpdir={config[tmpdirprefix]} XXXXXXglimpse2chr`; then
            echo "Failed to make tempdir!"
            exit 1
        fi
        
        samplename=`basename {input.bam} .bam`
        ploidy=`cut -f4 {input.xploidy}`
        echo "$samplename $ploidy" | tr ' ' '\t' > {output.samplefile}
        
        for reference in {params.panelprefix}.chrX*bin; do
            bn=`basename $reference`
            chrom=`echo $bn | cut -f2 -d '_'`
            start=`echo $bn | cut -f3 -d '_'`
            stop=`echo $bn | cut -f4 -d '_'`
            echo "{params.phasebinary} {params.phaseopt} --samples-file {output.samplefile} --bam-file {input.bam} --reference $reference --output $outdir/$chrom.$start.$stop.bcf "
        done | parallel -j {threads} &> {log.phase}
        
        # note that GLIMPSE_ligate fails on haploid samples. We use bcftools concat --ligate instead
        {params.bcftoolsbinary} concat --ligate -o {output.vcf} -Oz9 `ls -1v $outdir/chrX*bcf` &> {log.ligate}
        #{params.ligatebinary} --input <(ls -1v $outdir/chrX*bcf) --output {output.vcf}
        {params.bcftoolsbinary} index --threads {config[bcftoolsParams][threads]} {output.vcf} &>> {log.ligate}
        rm -rf $outdir  
        """


rule get_x_ploidy:
    input: 
        "{samplename}.{suffix}.merged.md.bam",
    output:
        "VCFs/{samplename}.{suffix}.xploidy.txt"
    log:
        "Logs/{samplename}.{suffix}.xploidy.txt"
    params:
        samtoolsbinary=GX.getParam("samtools", "binary"),
        aregion=GX.getParam("sexchromosomeParams", "aregion"),
        xregion=GX.getParam("sexchromosomeParams", "xregion"),
        binary=GX.getParam("sexchromosomeParams", "xestimator"),
        flags=GX.getParam("sexchromosomeParams", "samflags") # removes duplicates, read 2 (if it exists), secondary maps and fails QC
    shell: 
        """
        xreg=`cat {params.xregion}`
        areg=`cat {params.aregion}`
        
        #xreads=`{params.samtoolsbinary} view -q 20 -c -F {params.flags} {input} $xreg`
        #areads=`{params.samtoolsbinary} view -q 20 -c -F {params.flags} {input} $areg`
        #{params.binary} $areads $xreads > {output} 2> {log}
        
        {params.binary} `{params.samtoolsbinary} view -q 20 -c -F {params.flags} {input} $areg` `{params.samtoolsbinary} view -q 20 -c -F {params.flags} {input} $xreg` > {output} 2> {log}
        """    


# calls genotypes using bcftools.
# only weakly parallelized (but it's fast enough, so let's not bother)
rule bcftools_genotype_x:
    input:
        bam="{samplename}.{suffix}.merged.md.bam",
        ploidy="VCFs/{samplename}.{suffix}.xploidy.txt"
    log:
        "Logs/{samplename}.{suffix}.chrX.bcftools.log"
    output:
        "VCFs/bcftools/{samplename}.{suffix}.chrX.bcftools.vcf.gz"
    threads: 
        config["bcftoolsParams"]["threads"]
    params:
        ref=hg38,
        binary=bcftools,
        panel=GX.getParam("bcftools", "xpanel"),
        mpileopt=config["bcftoolsParams"]["mpileup_options"]
    shell: # needs testing. I may have to specify the sex, then provide a ploidy file.
        """
        ploidy=`cut -f4 {input.ploidy}`
        ({params.binary} mpileup {params.mpileopt} --fasta-ref {params.ref} -r chrX {input.bam} |\
        {params.binary} call --ploidy $ploidy -Am -P 0. -C alleles -T {params.panel} -Ov9 --threads {threads} -o {output}) 2> {log} &&\
        {params.binary} index --threads {config[bcftoolsParams][threads]} {output} 2>> {log}
        """

# This rule is NOT implemented; future versions may debug and implement "matchiness" filtering
# (what this does), but for now using Bayes factors appears to do well enough.
rule generate_sharing:
    input:
        "VCFs/glimpse2/{samplename}.{suffix}.glimpse2.vcf.gz"    
    output:
        seg="VCFs/glimpse2/{samplename}.{suffix}.ibisseg.tsv.gz",
        coef="VCFs/glimpse2/{samplename}.{suffix}.ibiscoef.tsv.gz",
    log:
        phase="Logs/{samplename}.{suffix}.matchiness.log"
    threads:
        config["matchinessParams"]["ibis_threads"]
    params:       
        ibis_binary= GX.getParam("matchinessParams","ibis_binary"),
        ibis_script= GX.getParam("matchinessParams", "ibis_script"),
        plink_binary= GX.getParam("matchinessParams","plink_binary"), # plink version 1.9 (NOT 2.x)!
        make_filters= GX.getParam("matchinessParams","filterer"),
        panel=GX.getParam("matchinessParams", "panel"),
        gmap=GX.getParam("matchinessParams", "gmap"),
        gsa_vcf=GX.getParam("matchinessParams", "gsavcf"),
        bcftoolsbinary= bcftools  
    shell: 
        """
        set +o pipefail
        if ! outdir=`mktemp -d --tmpdir={config[tmpdirprefix]} XXXXXXibis`; then
            echo "Failed to make tempdir!"
            exit 1
        fi
        # apply a barrage of filters; restrict to the GSA (autosomes)
        {params.bcftoolsbinary} view -V indels,mnps -R {params.gsa_vcf} {input} | {params.make_filters} | {params.bcftoolsbinary} norm -m+ | {params.bcftoolsbinary} view -M 2 -Oz1 -o ${outdir}/mfilt.vcf.gz
               
        # convert to plink file format
        {params.plink_binary} --vcf ${outdir}/mfilt.vcf.gz --set-missing-var-ids '@:#[b38]\$1,\$2' --make-bed --out $outdir/plink
            
        # TODO FIX (remove {} and/or replace AWK with something else)
        awk '{print $1,$1":"$4,$3,$4,$5,$6}' $outdir/plink.bim> $outdir/foo && mv $oudir/foo $outdir/plink.bim      
        {params.plink_binary} --bfile $outdir/plink --extract {params.panel}.bim --make-bed --out $outdir/plink_overlappedGSA

        # we have multiallelic variants. Let's deal with those.
        # generates merging failures. exit code 3 (from multiallelics)
        {params.plink_binary} --bfile $outdir/plink --bmerge {params.panel} --make-bed --out $outdir/merged_plink
        # remove merging failures
        {params.plink_binary} --bfile $outdir/plink_overlappedGSA --exclude $outdir/merged_plink.missnp --make-bed --out $outdir/plink_nomiss_overlappedGSA
        # from both files...
        {params.plink_binary} --bfile {params.panel} --exclude $outdir/merged_plink.missnp --make-bed --out $outdir/1kg
 
        {params.plink_binary} --bfile $outdir/1kg --bmerge $outdir/plink_nomiss_overlappedGSA --make-bed --out $outdir/
        
        # add the cM positions 
        {params.ibis_script} -cm $outdir/merged.bim {params.gmap} > $outdir/foo && mv foo $outdir/merged.bim
        
        {params.ibis_binary} -bfile $outdir/merged -printCoef -t {threads} -o $outdir/ibis
        
        cat $outdir/ibis.coef | gzip -9 > {output.coef}&
        cat $outdir/ibis.seg | gzip -9 > {output.seg}
        wait
        
        rm -rf $outdir
        """
        
# Note these are *variant* calls; consult the gvcf file if you wish
# to recover sites that are homozygous reference. 
rule deepvar:
    input:
        "{samplename}.{suffix}.merged.md.bam"
    output:
        vcf = "VCFs/deepvariant/{samplename}.{suffix}.deepvariant.vcf.gz",
        gvcf = "VCFs/deepvariant/{samplename}.{suffix}.deepvariant.g.vcf.gz"
    log:
        "Logs/{samplename}.{suffix}.deepvariant.log"
    threads: config["deepvarParams"]["dvthreads"]
    params:
        singularityimage=GX.getParam("deepvar", "image"),
        binary=GX.getParam("deepvar", "binary"), # note that this is a path WITHIN the singularity image
        path=os.path.abspath( os.path.join(GX.ROOT, "..")),
        ref=hg38,
        tmppath=config["tmpdirprefix"]
    run:
        shell(
            """
	    if ! outdir=`mktemp -d --tmpdir={params.tmppath} XXXXXXdv`; then
               echo "Failed to make tempdir for deepvariant!"
               exit 1
            fi
	    
            singularity run -B /usr/lib/locale,{params.path},{params.tmppath} {params.singularityimage} \
            {params.binary} --model_type=WGS \
            --ref={params.ref} \
            --reads={input} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards={threads} --intermediate_results_dir $outdir &> {log}           
            rm -rf $outdir
            """
        )


    """
    Our workflow is in hg38; gedmatch is in hg19...
    this fixes that problem :)
    IE, this runs GATK/Picard's liftovervcf tool
    by default, snps associated strand flips (A in hg38, T in hg19) are retained
    """
rule liftover2hg19:
    input:
        "VCFs/{caller}/{samplename}.{suffix}.{caller}.vcf.gz"
    output:
        vcf="VCFs/{caller}/{samplename}.{suffix}.{caller}.hg19.vcf.gz",
        fail="VCFs/{caller}/{samplename}.{suffix}.{caller}.liftover_fail.vcf.gz"
    params:
        binary=GX.getBinary("gatk"),
        options=config["gatkParams"]["liftoverOptions"],
        hg19=GX.getResources("Hg19"),
        chain=GX.getParam("gatk", "liftoverChain")
    log:
        "Logs/{samplename}.{suffix}.{caller}.liftover.log"
    shell:
        "{params.binary} --java-options {config[gatkParams][javaoptions]} LiftoverVcf {params.options} -C {params.chain} -I {input} -O {output.vcf} -R {params.hg19} --REJECT {output.fail} 2> {log}"

def getVcfFilters(wildcards):
    return config["bcf223andmeParams"][wildcards.caller]
        

rule hg19_to_23andme:
    input:
        "VCFs/{caller}/{samplename}.{suffix}.{caller}.hg19.vcf.gz"    
    output:
        "Uploads/{caller}/{samplename}.{suffix}.{caller}.hg19.23andme.tsv.gz"
    log:
        "Logs/{samplename}.{suffix}.{caller}.23andme.log"
    params:
        bcftools=bcftools,
        binary=GX.getBinary("bcf223andme"),
        filt=getVcfFilters
    shell: # in words; remove no-calls (-U), and places where the original and lifted chromosome do not match are excluded (-e ...)
           # merge (normalize) to reinstate multiallelic sites
           # filter the genotypes (depending on the caller; GQ20 equivalent by default)
           # and create a "23andme" file (using bcf223andme.py)
        """
        ({params.bcftools} view -Ou -U -e "INFO/OriginalContig!=CHROM" {input} | {params.bcftools} norm -Ov -m+  | \
        {params.binary} {params.filt} -o | gzip -9 > {output} )  2> {log}
        """



