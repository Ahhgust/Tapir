import glob, os, sys, re, time, gzip
import GenomixHelper

VERSION=0.30
# List of changes
# Added support for the MiSeq (by changing the CopyComplete.txt dependency when extracting BCLs)

ROOT=os.path.abspath(  os.path.join(workflow.current_basedir, ".."))

# defaults
configfile: os.path.join(ROOT, "configs", "config_v_2_standard.yaml")
configpath=os.path.join(ROOT, "configs", "config_v_2_standard.yaml")
i=0


# They python variable `configpath` is equivalent to the snakemake variable `configfile`
# (though this script cannot access `configfile` directly)
if '--configfile' in sys.argv:
    i = sys.argv.index('--configfile')

if i:
    if i < len(sys.argv)-1:
        configpath=sys.argv[i+1]
    else: # I do not think it's possible to enter this Else statement, but just in case.
        print("--configfile argument is misspecified (?)", file=sys.stderr)
        exit(1) # 
    

# GX provides an alternative way to get information from the config file.
# In practice, it creates the illusion of full paths (when they are needed),
# by stitching together that path based on the directory of this script.
# it also prepends python3/Rscript where needed
GX= GenomixHelper.GenomixHelper(workflow.current_basedir, config)
# you can still use the config dictionary directly, but only do so
# when want you want is independent of the path

# we use this a lot; let's just pull it out.
# reference genome
hg38=GX.getResources("Reference")

# bcftools binary:
bcftools=GX.getParam("bcftools", "binary")

# samtools binary
samtools=GX.getParam("samtools", "binary")


def sanitizeString(s):
    """
    strips out non-alphanumeric symbols from a string s (input), and returns it    
    """
    return "".join(c for c in s if c.isalnum() or c == "_" or c == '-' or c == '/')

def removeTrailingSlash(s):
    """
    Converts a foo/ or foo/bar/ to foo or foo/bar, respectively
    """
    if len(s) and s[-1] == "/":
        return s[:-1]
    return s    
    
SAMPLES=dict()
OFFTARGETS=dict() # eg, UDIs that you extract but are not associated with a sample
LIBRARIES=dict()
SAMPLESHEET=config["Samplesheet"]
if "SampleSheet" in config:
    SAMPLESHEET=config["SampleSheet"]

EXPERIMENT=""
    
# TODO: ADD CHECK FOR FASTQ INPUT
if not os.path.isdir(config["Bcldir"]) and "Fastqdir" not in config:
    print("You specified a directory that doesn't exist: ", config["Bcldir"], file=sys.stderr)
    exit(1)
    
# FYSA, snakemake does not like trailing slashes;
# even if they're syntactically preferred (as a unix convention)
BCLDIR=removeTrailingSlash( config["Bcldir"] )
RUN=os.path.basename(BCLDIR)

# this is a lie, but a necessary one...
# in short, we make Fastq input look like BCL input (but after the Fastqs have been extracted from BCLs)
if "Fastqdir" in config:
    if not os.path.isdir(config["Fastqdir"]):
        print("You're trying to use FASTQ files (shame on you). And you failed to specify a (valid/real) *directory*. Please fix the directory issue...", config["Fastqdir"], file=sys.stderr)
        exit(1)       
    BCLDIR=removeTrailingSlash( config["Fastqdir"] )
    RUN=os.path.basename(BCLDIR)   

# if the sample sheet is not found
if not os.path.isfile(SAMPLESHEET):   
    # check and see if it's in the BCL directory
    SAMPLESHEET=os.path.join(config["Bcldir"], SAMPLESHEET)
    if not os.path.isfile(SAMPLESHEET):   
        print(config["Samplesheet"] , " does not exist...", SAMPLESHEET, file=sys.stderr)
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
    ret = "bcl2fastq" # default to bcl2fastq
    with open(sheet) as fh:
        for line in fh:
            line=line.rstrip()
            if not all(ord(c) < 128 for c in line):
                print("Your Samplesheet is not an ASCII-text file. That's not okay! If you're exporting from Excel, export as 'CSV (MS-DOS)')", file=sys.stderr)
                exit(1)
            if len(line) < 1: # tolerate extra whitespace (typically at the end of the file)
                continue
                
            if line.startswith("[Data]"):
                getHeader=True
            elif not getHeader:
                if EXPERIMENT=="":
                    if line.startswith("Experiment Name"):
                        # this is a CSV element. grab the first non-empty element.
                        s = line.split(",")[1:]
                        for elem in s:
                            if elem != "":
                                EXPERIMENT=removeTrailingSlash( sanitizeString(elem) )
                                break
                if line.startswith("AdapterRead1,"): # one of the handful of syntactic differences between a bcl-convert sample sheet and a bcl2fastq sample sheet
                    ret="bcl-convert"
            elif getHeader and not getData:
                getData=True
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
                elif name.find(" ") >= 0:
                    print("Spaces are not allowed in sample identifiers. Maybe use a - instead?", name, file=sys.stderr)
                    errs+=1
                elif name.find(".") >= 0:
                    print("Periods are not allowed in sample identifiers. Maybe use a - instead?", name, file=sys.stderr)
                    errs+=1
                duptest.add(name)
    if errs:
        print("Your sample sheet has problems. Please fix them!", file=sys.stderr)
        exit(1)
    
    if EXPERIMENT == "":
        print("Failed to parse the Experiment Name (config[Experiment]) from the sample sheet", file=sys.stderr)
        exit(1)
    print( " ".join(duptest) , " are the samples in the sample sheet")    
    return ret

def isEmptyGzip(f):
    """
    returns True/False iff a gzipped file is "empty" (may have nonzero file size, but no contents)
    """
    if os.path.getsize(f) == 0:
        return True
    h= gzip.open(f, "rb")
    empty=True
    for line in h:
        empty=False
        break
    h.close()
    return empty

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
    i=-1
    for x in SAMPLENAME:
        i += 1
        if re.match('.*_?.*_L\d+.*',x) is not None:
            if x.find("_R2_") > -1: # we ignore read 2 (and test for its presence later on)
                continue
            # bcl2fastq does not emit empty fastqs; bcl-convert does. This effectively skips the empties
            if isEmptyGzip(files[i]): # was x, but we need the path... TODO: double check w/ bcl-convert
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
    #print("ParseFastqs: ", basedir, "\n", validFiles)             
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
        parseFastqs( os.path.join( wildcards.outdir, "Run_Fastqs", wildcards.rundir, "Fastqs"))
    
    if wildcards.samplename in SAMPLES:
        inner = SAMPLES[wildcards.samplename]
    elif wildcards.samplename in OFFTARGETS:
        inner = OFFTARGETS[wildcards.samplename]
    else:
        print("glob_bams error; unexpected sample_name", wildcards.samplename, file=sys.stderr)
        print("Samples are: ", SAMPLES)
        print("Wildcards are: ", wildcards)
        exit(1)

    files = expand("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}_{sampleindex}_{lane}.bam", samplename=wildcards.samplename, \
        sampleindex= inner['sampleindex'], lane=inner['lanes'], dirname=inner['dirname'],outdir=wildcards.outdir, rundir=wildcards.rundir, allow_missing=True)
    
  
    return files

def gather_all_reports(wildcards):
    """
    Note this is run *after* the checkpoint
    this emits all files that should be created (at the level of the sample)
    """
    global RUN
    global EXPERIMENT
    print("WILDCARDS", wildcards)
    
    if "outdir" not in wildcards:
        wildcards.outdir=EXPERIMENT 
    if "rundir" not in wildcards:
        wildcards.rundir=RUN

    fqdir = os.path.join(wildcards.outdir, "Run_Fastqs", wildcards.rundir, "Fastqs") 
    parseFastqs(fqdir)

    files = []
    
    for samp in SAMPLES.keys():
        files.extend( \
            expand("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam", \
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname=samp, samplename=samp)) # markdup-bams
        files.extend( \
            expand("{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.{suffix}", \
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname=samp, samplename=samp,\
            suffix=["bqsr.demix.summary", "r1_fastqc.zip", "flagstat", "samstats.cov", "bqsr_summary.pdf"]))
            #mixure analysis, fastqc, samtools flagstat, samstats (coverage estimate), and bqsr-post hoc summary

            
    for samp in OFFTARGETS.keys():
        files.extend( \
            expand("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam", \
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname="Offtargets", samplename=samp))# markdup-bams
        files.extend( \
            expand("{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.flagstat",
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname="Offtargets", samplename=samp)) # flagstat
        files.extend( \
            expand("{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats.cov",
            outdir=wildcards.outdir, rundir=wildcards.rundir, dirname="Offtargets", samplename=samp)) # samstats 
    print(files)  

    return files
    
def bcl_complete(wildcards):
    """
    helper function; 
    The Novaseq makes a "CopyComplete.txt" file when, well, the copying is complete
    while the MiSeq makes no such file. RTAComplete seems like equivalent, though
    """
    global BCLDIR
    file=os.path.join(BCLDIR, "CopyComplete.txt")
    if os.path.isfile(file):
        return file
    return os.path.join(BCLDIR, "RTAComplete.txt")

def touchy(filename):
    """
    unix touch (but not really)
    makes an empty file iff it doesn't exist.
    does NOT update the timestamp (unlike unix's touch)
    """
    if not os.path.isfile(filename):
        han=open(filename, "w")
        han.close()


# EXPERIMENT set here (from sample sheet)
bclConverter=parseSamplesheet(SAMPLESHEET) # determines the name of the Experiment and the library prep; performs a sanity check on the samples

# check to see if the exeperiment (name) has been overloaded in the command line
exp=removeTrailingSlash ( sanitizeString( config["Experiment"] )  )
# command line trumps the sample sheet.
if exp != "":
    EXPERIMENT= exp

# Fastq input?
if "Fastqdir" in config:
    print("Fastq input detected...")
    fastqs=glob.glob(os.path.join(BCLDIR, "*fastq.gz"))
    for fq in fastqs:
        sp= os.path.basename(fq).split("_")
        if len(sp) != 5:
            print("Fastq files MUST be made by bcl2fastq/bcl-convert (or made to look like it). The expetect format is: {samplename}_{sampleid}_{lane}_{R?}_001.fastq.gz", "not", fq, file=sys.stderr, sep="\n")
            exit(1)
            
    # a little unseemly; write into the Fastq directory.
    touchy( os.path.join(BCLDIR, "CopyComplete.txt") )
    d=os.path.join(EXPERIMENT, "Run_Fastqs", RUN, "Fastqs", "Stats")
    
    if os.system(f"mkdir -p {d}"):
        print(f"Failed to make directory: {d}")
        exit(1)
    
    # make the Fastq directory (as if Tapir had made it from the BCL files...)
    #touchy( os.path.join(d, "Stats.json") ) # made by rule bcl:
    d=os.path.dirname(d) # trim off the Stats
    for fq in fastqs:
        apath=os.path.abspath(fq)
        if os.system(f"ln -sf {apath} {d}"):
            print("Failed to make symlink", apath, d, sep="\n", file=sys.stderr)
            exit(1)
    
    # last but not least, put the sample sheet in the correct location
    d=os.path.dirname(d) # trim off the Fastqs/
    d=os.path.join(d, "SampleSheet.csv")
    os.system(f"cp -L {SAMPLESHEET} {d}")
    bclConverter='none'
    
    
print("Converter: ", bclConverter)
print("Experiment: ", EXPERIMENT)
print("Bcldir: ", BCLDIR)
print("Samplesheet: ", SAMPLESHEET)


rule bcl:
    input:
        expand("{outdir}/Run_Fastqs/{rundir}/SampleSheet.csv", outdir=EXPERIMENT, rundir=RUN), # extract fastqs (successfully)...
        expand("{outdir}/Run_Fastqs/{rundir}/RunAnalysisComplete.txt", outdir=EXPERIMENT, rundir=RUN)
    # set to either bcl2fastq or bcl-convert, depending on the format of the sample sheet
if bclConverter=="bcl2fastq":
    checkpoint extract_fastqs:
        input:
            bcl_complete
        params:
            rl=GX.getParam("bcl2fastq", "params"),
            samplesheet=SAMPLESHEET,
            bclpath=BCLDIR,
            binary=GX.getBinary("bcl2fastq"),
        threads:
            config["bcl2fastqParams"]["threads"]
        log:
            "{outdir}/Run_Fastqs/{rundir}/Fastqs/bcl2fastq.outerr"
        output: 
            "{outdir}/Run_Fastqs/{rundir}/Fastqs/Stats/Stats.json"
        shell:            
            """
            bn=`dirname {output}`  
            bn=`dirname $bn` # removes Stats/Stats.json
            {params.binary} {params.rl} -R {params.bclpath} -w 0 -r {threads} -p {threads} -o $bn --sample-sheet {params.samplesheet} 2> {log}
            """
 
elif bclConverter=="bcl-convert":
    checkpoint extract_fastqs:
        input:
            bcl_complete
        params:
            samplesheet=SAMPLESHEET,
            bclpath=BCLDIR,
            obinary=GX.getBinary("bclconvert"),
            oparams=GX.getParam("bclconvert", "params")
        threads:
            config["bclconvertParams"]["threads"]
        log:
            "{outdir}/Run_Fastqs/{rundir}/Fastqs/bcl-convert.outerr"
        output: 
            "{outdir}/Run_Fastqs/{rundir}/Fastqs/Logs/FastqComplete.txt"
        shell:
            """
            bn=`dirname {output}`  
            bn=`dirname $bn` # removes Logs/FastqComplete.txt
            {params.obinary} {params.oparams} --bcl-input-directory {params.bclpath} --output-directory $bn --bcl-num-conversion-threads {threads} --sample-sheet {params.samplesheet} 2> {log}
            """
elif bclConverter=='none': # in this case, the fastqs are already made...
    checkpoint extract_fastqs:
        input:
            bcl_complete
        params:
            rl=GX.getParam("bcl2fastq", "params"),
            samplesheet=SAMPLESHEET,
            bclpath=BCLDIR,
            binary=GX.getBinary("bcl2fastq"),
        threads:
            1
        log:
            "{outdir}/Run_Fastqs/{rundir}/Fastqs/bcl2fastq.outerr"
        output: 
            "{outdir}/Run_Fastqs/{rundir}/Fastqs/Stats/Stats.json"
        shell:            
            """
            touch {output} 2> {log} 
            """
else:
    print("Should never happen; failed to deduce the right bcl to fastq converter!?", file=sys.stderr)
    exit(1)
    
rule gather_fastqs:
    input:
        aggregate_fastqs # updates the SAMPLES and OFFTARGETS dictionaries.
    output:
        samplesheet="{outdir}/Run_Fastqs/{rundir}/SampleSheet.csv",
        version="{outdir}/Run_Fastqs/{rundir}/version.txt"
    params:
        samplesheet=SAMPLESHEET,
        version=config["Version"]
    shell:
        """
        cp -L {params.samplesheet} {output.samplesheet} && echo "bcl2bam version " {params.version} > {output.version}
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
        fq="{outdir}/Run_Fastqs/{rundir}/Fastqs/{samplename}_{sampleindex}_{lane}_R1_001.fastq.gz",
        samplesheet="{outdir}/Run_Fastqs/{rundir}/SampleSheet.csv" # builds dependency to gather_fastqs
    output:
        bam=temp("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}_{sampleindex}_{lane}.bam") # note; no indexing required.
    wildcard_constraints:
        lane="L\d+"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}_{sampleindex}_{lane}_bwa_mem.log"   
    threads: config["bwamemParams"]["threads"]
    resources:
        mem_mb=config["samtoolsParams"]["mem_mb_sort"]
    params:
        samtools_binary=samtools,
        bwa_binary=GX.getBinary("bwamem"),
        ref=hg38,
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
			-M {params.ref} \
			{read1} {read2} | \
			{params.samtools_binary} sort -@ {config[samtoolsParams][samtoolsThreads]}  -m {config[samtoolsParams][sortMem]} - > {output.bam} ) 2> {log}
			"""
        )

       

# alpha version of the pipeline had two separate steps; merge, then left align.
# leftalign takes multiple input files, soo, let's combine those operations.
#{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}_{sampleindex}_{lane}.bam"
rule merge_and_leftalign_bams:
    input:
        glob_bams
    output:
        bam=temp("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.bam"),
        bai=temp("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.bai")
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.merge_leftalign.log"
    params:
        binary=GX.getBinary("gatk"),
        ref=hg38
    run:
        inputs_all= " ".join(["-I " + repr(i) for i in input])
        shell("{params.binary} --java-options {config[gatkParams][javaoptions]} LeftAlignIndels --verbosity ERROR --QUIET true -R {params.ref} {inputs_all} -O {output.bam} &> {log}")   



rule markdup_bams:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.bam"
    output: temp("{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam")
    threads: config["sambambaParams"]["threads"]
    params:
        tmpprefix=config["tmpdirprefix"],
        binary=GX.getBinary("sambamba")
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.markdup.log"   
    shell:
        "{params.binary} markdup  --tmpdir={params.tmpprefix} -l 9 -t {threads} {input} {output} 2> {log}"
 

rule make_bam_flagstats:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.flagstat"
    threads: 
        config["samtoolsParams"]["samtoolsThreads"]
    params:
        binary=samtools
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.flagstat.log"
    shell:
        "{params.binary} flagstat -@ {threads}  {input} > {output} 2> {log}"

rule make_bam_samstats:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.la.md.samstats.log"
    params:
        binary=GX.getBinary("samstats"),
        panel=GX.getParam("samstats", "panel")
    shell:
        "{params.binary} {input} {params.panel} > {output} 2> {log}"

rule make_bam_cov:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.samstats.cov"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.la.md.samstatscov.log"
    params:
        binary=GX.getSummarizer("samstats")
    shell:
        "{params.binary} {input}  > {output} 2> {log}"

rule bam_estimate_mix:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bqsr.bam"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.bqsr.demix"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.la.md.mf.bqsr.log"
    params:
        binary=GX.getBinary("demixtify"),
        panel=GX.getParam("demixtify", "panel")
    threads: 
        config["demixtifyParams"]["threads"]
    shell:
        "{params.binary} -t {threads} -b {input} -v  {params.panel} > {output} 2> {log}"
    
rule mix_summary:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.bqsr.demix"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.bqsr.demix.summary"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.la.md.mf.bqsr.summary.log"
    params:
        binary=GX.getSummarizer("demixtify")
    shell:
        "{params.binary} {input} > {output} 2> {log}"

# note: as written, this cannot be applied to the off-target bams. 
rule fastqc:
    input:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        r2="{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.r2_fastqc.html",
        o2="{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.r2_fastqc.zip",
        r1="{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.r1_fastqc.html",
        o1="{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.r1_fastqc.zip"
    params:
        samtools_binary=samtools,
        fastqc_binary=GX.getBinary("fastqc")
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.fastqc.log"
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
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.recal.csv"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.bqsr.recal.log"
    params:
        binary=GX.getBinary("gatk"),
        args=GX.getParam("gatk", "bqsrargs"),
        mask=GX.getParam("gatk", "bqsrmask"),
        ref=hg38
    shell:
        "{params.binary} --java-options {config[gatkParams][javaoptions]} BaseRecalibrator {params.args} -R {params.ref} -O {output} -I {input}  --known-sites {params.mask} &> {log}"

rule gatk_apply_bqsr:
    input:
        bam="{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bam",
        recal_file="{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.recal.csv"
    output:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bqsr.bam"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.bqsr.apply.log"
    params:
        binary=GX.getBinary("gatk"),
        ref=hg38,
        compression_level=9
    shell:# the ln -s command makes a foo.bam.bai index file ; by default, GATK makes their index files as foo.bai instead (which breaks some downstream tools)
        """
        {params.binary} --java-options '{config[gatkParams][javaoptions]} -Dsamjdk.compression_level={params.compression_level}' ApplyBQSR --QUIET true --create-output-bam-index true -R {params.ref} -O {output} --bqsr-recal-file {input.recal_file} -I {input.bam} &> {log}
        # makes a .bam.bai version of the .bai file, too. 
        cp `echo {output} | rev | cut -b2- | rev`i {output}.bai 
        """

# note: requires R packages: gplots and gsalib (legacy code; may be difficult to install)
rule plot_gatk_bqsr:
    input:
        bam="{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.la.md.bqsr.bam",
        recal_file="{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.recal.csv"
    output:
        report="{outdir}/Sample_Data/{dirname}/{rundir}/Reports/{samplename}.la.md.bqsr_summary.pdf",
        posthoc_file="{outdir}/Sample_Data/{dirname}/{rundir}/Bams/{samplename}.recalafter.csv"
    log:
        "{outdir}/Sample_Data/{dirname}/{rundir}/Logs/{samplename}.bqsr.plots.log"            
    params:
        binary=GX.getBinary("gatk"),
        args=GX.getParam("gatk", "bqsrargs"),
        mask=GX.getParam("gatk", "bqsrmask"),
        ref=hg38
    shell:
        """
        {params.binary} --java-options {config[gatkParams][javaoptions]} BaseRecalibrator {params.args} -R {params.ref} -O {output.posthoc_file} -I {input.bam} --known-sites {params.mask} &> {log}
        {params.binary} --java-options {config[gatkParams][javaoptions]} AnalyzeCovariates --QUIET --before {input.recal_file} --after {output.posthoc_file} --plots {output.report} &>> {log}
        """    
  

# for reproducibility
# record the config file used, the exact command line arguments, and the version number (this script)
rule repro:
    input:
        aggregate_fastqs,
        gather_all_reports
    output:
        marker="{outdir}/Run_Fastqs/{rundir}/RunAnalysisComplete.txt",
        args="{outdir}/Run_Fastqs/{rundir}/args.txt",
        cfile="{outdir}/Run_Fastqs/{rundir}/config.csv"
    params:
        configFile=configpath,
        versionNumber=VERSION,
        arguments=" ".join(sys.argv)
    shell:
        """
        echo {params.arguments} > {output.args}
        cp -L {params.configFile} {output.cfile}
        echo "Version: " {params.versionNumber} > {output.marker}
        """

        