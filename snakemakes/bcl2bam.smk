import glob, os, sys, re, time, gzip
configfile: "/eva/codebase/snakemakes/config_v2.yaml"

def sanitizeString(s):
    """
    strips out non-alphanumeric symbols from a string s (input), and returns it    
    """
    return "".join(c for c in s if c.isalnum() or c == "_")



SAMPLES=dict()
OFFTARGETS=dict() # eg, UDIs that you extract but are not associated with a sample
LIBRARIES=dict()
EXPERIMENT=sanitizeString( config["Experiment"] ) # updated in parseSampleSheet
# this is the PARENT directory of where the output will be written
# OUTDIR/EXPERIMENT is taken as the root directory
OUTDIR=config["Outdir"]
SAMPLESHEET=config["Samplesheet"]

BCLPATH= os.getcwd()

if config["Bcldir"] != ".":
    if not os.path.isdir(config["Bcldir"]):
        print("You specified a directory that doesn't exist: ", config["Bcldir"], file=sys.stderr)
        exit(1)
    BCLPATH=config["Bcldir"]
    if OUTDIR=="..":
        OUTDIR=os.path.join(BCLPATH, OUTDIR)

# the name of the BCL run directory.
# eg, 240325_A01324_0104_BH2FFYDSXC
BCLDIR= os.path.basename(BCLPATH )

# if the sample sheet is not found
if not os.path.isfile(SAMPLESHEET):   
    # check and see if it's in the BCL directory
    SAMPLESHEET=os.path.join(config["Bcldir"], SAMPLESHEET)
    if not os.path.isfile(SAMPLESHEET):   
        print(config["Samplesheet"] , " does not exist...", file=sys.stderr)
        exit(1)

def parseSamplesheet(sheet):
    """
    Takes in a sample sheet (SampleSheet.csv is most common)
    and does sanity testing on it. 
    further, it deduces the names of the libraries and the name of the experiment
    (globals)
    """
    getData=False
    getHeader=False
    idCol = -1
    nameCol=-1
    libCol=-1
    errs=0
    duptest= set()
    global EXPERIMENT
    global LIBRARIES
    
    with open(sheet) as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith("[Data]"):
                getHeader=True
            elif not getHeader and EXPERIMENT=="":
                if line.startswith("Experiment Name"):
                    # this is a CSV element. grab the first non-empty element.
                    s = line.split(",")[1:]
                    for elem in s:
                        if elem != "":
                            EXPERIMENT=sanitizeString(elem)
                            break
                
            elif getHeader:
                getData=True
                getHeader=False
                sp = line.split(",")
                i=0
                for colid in sp:
                    if colid== "Sample_ID":
                        idCol=i
                    elif colid=="Sample_Name":
                        nameCol=i
                    
                    if colid==config["SamplesheetLibraryColumn"]:
                        libCol=i
                    i += 1
                if idCol< 0 and nameCol < 0:
                    print("Your sample sheet is corrupt: failed to extract IDs/Names from: ", sheet, sep="\n", file=sys.stderr)
                    exit(1)
                if libCol<0:
                    print("Failed to extract the library from", sheet, config["SamplesheetLibraryColumn"], sep="\n", file=sys.stderr)
                    exit(1)
                    
            elif getData:
                sp = line.split(",")
                sid=sname=""
                if idCol >= 0:
                    sid = sp[idCol]
                
                if nameCol>=0:
                    sname=sp[nameCol]
                
                # if no sample name is given, the sample name is the sid
                name = sid
                if sname != "":
                    name = sname
                
                if name == "":
                    errs += 1
                    print("Blank sample name detected..", line, sep="\n", file=sys.stderr)
                    
                LIBRARIES[name]=sp[libCol]
                
                if name in duptest:
                    print("Duplicate sample names detected; that's a problem", name, line, sep="\n", file=sys.stderr)
                    errs+=1
                elif name.find("_") >= 0:
                    print("Underscores are not allowed in sample identifiers. Maybe use a - instead?", name, file=sys.stderr)
                    errs+=1
                duptest.add(name)
    if errs:
        print("Your sample sheet has problems. Please fix them!", file=sys.stderr)
        exit(1)
    
    if EXPERIMENT == "":
        print("Failed to parse the Experiment Name (config[Experiment]) from the sample sheet", file=sys.stderr)
        exit(1)
        


def parseFastqs(basedir):
    """
    Input: a directory name with *.fastq.gz files AS MADE BY BCL2FASTQ2
    Output: none; updates the SAMPLES and OFFTARGETS data structures.
    """
    files=glob.glob(os.path.join(basedir, "*fastq.gz")) # note: later on we only consider read1 (_R1_). If there's a read2, the BWA will find it.
    SAMPLENAME=[os.path.basename(x) for x in files]
    global SAMPLES
    global OFFTARGETS
    validFiles = []
    # makes a pairing; each sample name corresponds to 1+ lanes of data (novaseq S4 it's 4 lanes)
    # the sample name and lane information is determined by the fastq files. 
    for x in SAMPLENAME:
        if re.match('.*_?.*_L\d+.*',x) is not None:
            if x.find("_R2_") > -1: # we ignore read 2 (and test for its presence later on)
                continue
            
            # file format is:
            # UDI-002_S2_L002_R1_001.fastq.gz
            # samplename_index_lane_readnum_001.fastq.gz
            # note that if the sample name starts with UDI-
            # it's assumed to be off-target
            s = x.split("_")
            if len(s) < 4:
                print("Critical error: Your fastqs do not conform to the BCL2FASTQ standard", s, sep="\n", file=sys.stderr)
                print("Please rename your files to look like: samplename_index_lane.*001.fastq.gz", file=sys.stderr)
                exit(1)
                
            samp= s[0]
            
            if samp == "Offtargets":
                print("That's just diabolical. Choose a better sample name than 'Offtargets' please", file=sys.stderr)
                exit(1)
            
            d = SAMPLES
            if samp.startswith( config["Offtarget_index_prefix"] ):
                d = OFFTARGETS
            
            if samp not in d:
                d[samp] = dict()
                d[samp]['sampleindex']=""
                d[samp]['lanes']=set()
                d[samp]['dirname']=""
                
            
            inner = d[samp]
            # adding: sampleindex,lane
            inner['lanes'].add(s[2])
            inner['sampleindex']= s[1]
            
            if d == OFFTARGETS:
                inner['dirname'] = "Offtargets"
            else:
                inner['dirname'] = samp
            
            validFiles.append( os.path.join(basedir, x))
            
    return validFiles

def aggregate_fastqs(wildcards):
    """
    Helper function; used for the checkpoint
    it gathers the fastqs that were extracted, and sets the appropriate fields (associating the samplename w/ the set of fastqs)
    """


    # .output refers to: "Fastqs/Stats/Stats.json" ; we just want the Fastqs directory
    fqdir = os.path.dirname( os.path.dirname(str(checkpoints.extract_fastqs.get(**wildcards).output)) )
    
    validFiles = parseFastqs(fqdir) # update the SAMPLES and OFFTARGETs dictionaries.
    
    return validFiles
    


def glob_bams(wildcards):
    """
    helper function; used to generate the files (with relative path) of all bam files (offtargets too); bams
    are per-lane.
    """
    # it's unclear why I have to update the fastqs information (again)
    # but this is necessary to populate the SAMPLES and OFFTARGETS dictionaries.
    if len(SAMPLES) + len(OFFTARGETS)==0:
        parseFastqs( os.path.join( wildcards.outdir, wildcards.rundir, "Fastqs"))
    
    if wildcards.samplename in SAMPLES:
        inner = SAMPLES[wildcards.samplename]
    elif wildcards.samplename in OFFTARGETS:
        inner = OFFTARGETS[wildcards.samplename]
    else:
        print("glob_bams error; unexpected sample_name", wildcards.samplename, file=sys.stderr)
        print("Samples are: ", SAMPLES)
        print("Wildcards are: ", wildcards)
        exit(1)

    files = expand("{outdir}/{dirname}/{rundir}/Bams/{samplename}_{sampleindex}_{lane}.bam", samplename=wildcards.samplename, \
        sampleindex= inner['sampleindex'], lane=inner['lanes'], dirname=inner['dirname'],outdir=wildcards.outdir, rundir=wildcards.rundir, allow_missing=True)
    
  
    return files

def gather_all_reports(wildcards):
    """
    Note this is run *after* the checkpoint
    this emits all files that should be created (at the level of the sample)
    """

    fqdir = os.path.join(wildcards.outdir, wildcards.rundir, "Fastqs")   
    
    parseFastqs(fqdir)
    
    
    files = []
    
    for samp in SAMPLES.keys():
        files.extend( \
            expand("{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam", \
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname=samp, samplename=samp)) # markdup-bams
        files.extend( \
            expand("{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.{suffix}", \
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname=samp, samplename=samp,\
            suffix=["demix.summary", "r1_fastqc.zip", "flagstat", "samstats.cov", "bqsr_summary.pdf"]))
            #mixure analysis, fastqc, samtools flagstat, samstats (coverage estimate), and bqsr-post hoc summary

            
    for samp in OFFTARGETS.keys():
        files.extend( \
            expand("{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam", \
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname="Offtargets", samplename=samp))# markdup-bams
        files.extend( \
            expand("{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.flagstat",
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname="Offtargets", samplename=samp)) # flagstat
        files.extend( \
            expand("{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats.cov",
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname="Offtargets", samplename=samp)) # samstats 
    print(files)        
    return files
    
        
parseSamplesheet(SAMPLESHEET) # determines the name of the Experiment and the library prep; performs a sanity check on the samples


rule bcl:
    input:
        expand("{outdir}/{rundir}/SampleSheet.csv", outdir=os.path.join(OUTDIR,EXPERIMENT), rundir=BCLDIR), # extract fastqs (successfully)...
        expand("{outdir}/{rundir}/RunAnalysisComplete.txt", outdir=os.path.join(OUTDIR,EXPERIMENT), rundir=BCLDIR)
        
checkpoint extract_fastqs:
    input:
        os.path.join(BCLPATH, "CopyComplete.txt")
    params:
        rl=config["bcl2fastqParams"]["params"],
        samplesheet=SAMPLESHEET,
        bclpath=BCLPATH,
        binary=config["bcl2fastqParams"]["binary"]
    threads:
        config["bcl2fastqParams"]["threads"]
    log:
        "{outdir}/{rundir}/Fastqs/bcl2fastq.outerr"
    output: 
        "{outdir}/{rundir}/Fastqs/Stats/Stats.json"
    shell:
        """
        bn=`dirname {output}`  
        bn=`dirname $bn` # removes Stats/Stats.json
        {params.binary} {params.rl} -R {params.bclpath} -w {threads} -p {threads} -o $bn --sample-sheet {params.samplesheet} --fastq-compression-level 9 2> {log}
        """

rule gather_fastqs:
    input:
        aggregate_fastqs # updates the SAMPLES and OFFTARGETS dictionaries.
    output:
        samplesheet="{outdir}/{rundir}/SampleSheet.csv",
        version="{outdir}/{rundir}/version.txt"
    params:
        samplesheet=SAMPLESHEET,
        version=config["Version"]
    shell:
        """
        cp {params.samplesheet} {output.samplesheet} && echo "bcl2bam version " {params.version} > {output.version}
        """

def getLibrary(wildcards):
    """
    Helper function (for bwa_mem_map)
    returns the library identifier for a particular sample
    (or just the sample name if no library is found. e.g., the "Undetermined" fastqs/bams)
    """
    global LIBRARIES
    sn=str(wildcards.samplename)
    
    if sn in LIBRARIES:
        return LIBRARIES[sn]
    return sn

rule bwa_mem_map:
    input:
        fq="{outdir}/{rundir}/Fastqs/{samplename}_{sampleindex}_{lane}_R1_001.fastq.gz",
        samplesheet="{outdir}/{rundir}/SampleSheet.csv"
    output:
        bam=temp("{outdir}/{dirname}/{rundir}/Bams/{samplename}_{sampleindex}_{lane}.bam") # note; no indexing required.
    wildcard_constraints:
        lane="L\d+"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}_{sampleindex}_{lane}_bwa_mem.log"   
    threads: config["bwamemParams"]["threads"]
    resources:
        mem_mb=config["samtoolsParams"]["mem_mb_sort"]
    params:
        samtools_binary=config["samtoolsParams"]["binary"],
        bwa_binary=config["bwamemParams"]["binary"],
        library=getLibrary 
    run:
        read1=input[0]
        read2=input[0].replace('_R1_','_R2_')
        
        # limited support for single-end reads
        if not os.path.isfile(read2):
            read2=""
            print("No read2 detected; assuming single-end sequencing", file=sys.stderr)
            
        # the set +o pipefail is necessary for the zcat | head (which generates a fail)
        shell(
			"""
            set +o pipefail
			header=$(zcat {input.fq} | head -1)
			seqID=$(echo $header | cut -d':' -f 1 | cut -d' ' -f 1 | cut -c 2-)
			id=$(echo $header | cut -d':' -f 2-4 --output-delimiter='_')
			sm=$(echo {wildcards.samplename} | cut -d'_' -f 1)
			pu=$(echo $header | cut -d':' -f 3,4 --output-delimiter='_')
			rg=$(echo "@RG\tID:$seqID"_"$id\tSM:$sm\tLB:{params.library}\tPL:ILLUMINA")
			({params.bwa_binary} mem  -R $(echo "@RG\\tID:$sm"_"$seqID"_"$id\\tSM:$sm\\tLB:{params.library}\\tPL:ILLUMINA\\tPU:$pu") \
			-t {threads} \
			-M {config[refs][Reference]} \
			{read1} {read2} | \
			{params.samtools_binary} sort -@ {config[samtoolsParams][samtoolsThreads]}  -m {config[samtoolsParams][sortMem]} - > {output.bam} ) 2> {log}
			"""
        )

       

# alpha version of the pipeline had two separate steps; merge, then left align.
# leftalign takes multiple input files, soo, let's combine those operations.
#{outdir}/{dirname}/{rundir}/Bams/{samplename}_{sampleindex}_{lane}.bam"
rule merge_and_leftalign_bams:
    input:
        glob_bams
    output:
        bam=temp("{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.bam"),
        bai=temp("{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.bai")
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.merge_leftalign.log"
    params:
        binary=config["gatkParams"]["binary"]
    run:
        inputs_all= " ".join(["-I " + repr(i) for i in input])
        shell("{params.binary} LeftAlignIndels --verbosity ERROR --QUIET true -R {config[refs][Reference]} {inputs_all} -O {output.bam} &> {log}")   


rule markdup_bams:
    input:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.bam"
    output:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    threads: config["sambambaParams"]["threads"]
    params:
        tmpprefix=config["tmpdirprefix"],
        binary=config["sambambaParams"]["binary"]
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.markdup.log"   
    shell:
        "{params.binary} markdup  --tmpdir={params.tmpprefix} -l 9 -t {threads} {input} {output} 2> {log}"
 

rule make_bam_flagstats:
    input:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.flagstat"
    threads: config["samtoolsParams"]["samtoolsThreads"]
    params:
        binary=config["samtoolsParams"]["binary"]
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.flagstat.log"
    shell:
        "{params.binary} flagstat -@ {threads}  {input} > {output} 2> {log}"

rule make_bam_samstats:
    input:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.la.md.samstats.log"
    params:
        binary=config["samstatsParams"]["binary"]
    threads: 1
    shell:
        "{params.binary} {input} {config[samstatsParams][panel]}   > {output} 2> {log}"

rule make_bam_cov:
    input:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats"
    output:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats.cov"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.la.md.samstatscov.log"
    threads: 1
    params:
        binary=config["samstatsParams"]["summarizer"] 
    shell:
        "{params.binary} {input}  > {output} 2> {log}"

rule bam_estimate_mix:
    input:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.demix"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.la.md.mf.log"
    params:
        binary=config["demixtifyParams"]["binary"] 
    threads: 
        config["demixtifyParams"]["threads"]
    shell:
        "{params.binary} -t {threads} -b {input} -v  {config[demixtifyParams][panel]} > {output} 2> {log}"
    
rule mix_summary:
    input:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.demix"
    output:
        "{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.demix.summary"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.la.md.mf.summary.log"
    params:
        binary=config["demixtifyParams"]["summarizer"]
    shell:
        "{params.binary} {input} > {output} 2> {log}"

# note: as written, this cannot be applied to the off-target bams. 
rule fastqc:
    input:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        r2="{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.r2_fastqc.html",
        o2="{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.r2_fastqc.zip",
        r1="{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.r1_fastqc.html",
        o1="{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.r1_fastqc.zip"
    params:
        samtools_binary=config["samtoolsParams"]["binary"],
        fastqc_binary=config["fastqcParams"]["binary"],
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.fastqc.log"
    shell: # note: this cannot be applied to the UDI- data.
        """
        outdir=`dirname {output.r1}`
        {params.samtools_binary} fastq -0 /dev/null -2 /dev/null -1 /dev/stdout -F 2308 {input} | {params.fastqc_binary} --quiet --delete -o $outdir /dev/stdin 2> {log}
        mv $outdir/stdin_fastqc.html {output.r1}
        mv $outdir/stdin_fastqc.zip {output.o1}
        {params.samtools_binary} fastq -0 /dev/null -1 /dev/null -2 /dev/stdout -F 2308 {input} | {params.fastqc_binary} --quiet --delete -o $outdir /dev/stdin 2>> {log}
        mv $outdir/stdin_fastqc.html {output.r2}
        mv $outdir/stdin_fastqc.zip {output.o2}
        """
        # fastqc is a bit annoying to work with.
        # if you pipe it information, it'll write the file names as "stdin"
        # the above is done (in series); note that the wildcards.samplename removes the possibility of namespace collisions.
        # some would have you run fastqc on each individual (fastq) file. that is TMI
        # considering sample level (R1 and R2) make a lot more sense
        # it makes sense to only consider mapped reads, not supplemental alignments, and primary alignments (2308);  

rule gatk_learn_bqsr_table:
    input:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.recal.csv"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.bqsr.recal.log"
    params:
        binary=config["gatkParams"]["binary"],
        args=config["gatkParams"]["bqsrargs"],
        mask=config["atlasParams"]["bqsrmask"]
    shell:
        "{params.binary} --java-options {config[gatkParams][javaoptions]} BaseRecalibrator {params.args} -R {config[refs][Reference]} -O {output} -I {input}  --known-sites {params.mask}.gz &> {log}"

rule gatk_apply_bqsr:
    input:
        bam="{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bam",
        recal_file="{outdir}/{dirname}/{rundir}/Bams/{samplename}.recal.csv"
    output:
        "{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bqsr.bam"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.bqsr.apply.log"
    params:
        binary=config["gatkParams"]["binary"]
    shell:# the ln -s command makes a foo.bam.bai index file ; by default, GATK makes their index files as foo.bai instead (which breaks some downstream tools)
        """
        {config[gatkParams][binary]} --java-options '{config[gatkParams][javaoptions]} -Dsamjdk.compression_level=9' ApplyBQSR --QUIET true --create-output-bam-index true -R {config[refs][Reference]} -O {output} --bqsr-recal-file {input.recal_file} -I {input.bam} &> {log}
        # trims off the last letter in the bam (.bam -> bai), and sets that as the source of the symlink
        ln -s `echo {output} | rev | cut -b2- | rev`i {output}.bai 
        """

# note: requires R packages: gplots and gsalib (legacy code; may be difficult to install)
rule plot_gatk_bqsr:
    input:
        bam="{outdir}/{dirname}/{rundir}/Bams/{samplename}.la.md.bqsr.bam",
        recal_file="{outdir}/{dirname}/{rundir}/Bams/{samplename}.recal.csv"
    output:
        report="{outdir}/{dirname}/{rundir}/Reports/{samplename}.la.md.bqsr_summary.pdf",
        posthoc_file="{outdir}/{dirname}/{rundir}/Bams/{samplename}.recalafter.csv"
    log:
        "{outdir}/{dirname}/{rundir}/Logs/{samplename}.bqsr.plots.log"            
    params:
        binary=config["gatkParams"]["binary"],
        args=config["gatkParams"]["bqsrargs"],
        mask=config["atlasParams"]["bqsrmask"]
    shell:
        """
        {params.binary} --java-options {config[gatkParams][javaoptions]} BaseRecalibrator {params.args} -R {config[refs][Reference]} -O {output.posthoc_file} -I {input.bam} --known-sites {params.mask}.gz &> {log}
        {params.binary} --java-options {config[gatkParams][javaoptions]} AnalyzeCovariates --QUIET --before {input.recal_file} --after {output.posthoc_file} --plots {output.report} &>> {log}
        """    
  


rule complete:
    input:
        aggregate_fastqs,
        gather_all_reports
    output:
        "{outdir}/{rundir}/RunAnalysisComplete.txt"
    shell:
        """
        echo input
        touch {output}
        """

        