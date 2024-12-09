import glob, os, sys, re, time, gzip
configfile: "/eva/codebase/snakemakes/config_v2.yaml"

if not os.path.isdir("Fastqs"):
    print("I need a Fastqs/ directory to work...", file=sys.stderr)
    exit(1)


SAMPLES=dict()
OFFTARGETS=dict() # eg, UDIs that you extract but are not associated with a sample

def parseFastqs(basedir):
    """
    Input: a directory name with *.fastq.gz files AS MADE BY BCL2FASTQ2
    Output: a list of the valid filenames (with paths); this is often ignored.
    # Side effects: updates the SAMPLES and OFFTARGETS data structures.
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

parseFastqs("Fastqs")


def get_sample_bams(wildcards):
    # in this case, wildcards is empty; just populate what BAM files we'd like to see
    # note each bam is composed of 1+ Lanes worth of data.
    ret = expand("{samplename}/Bams/{samplename}.la.bam", samplename=SAMPLES.keys())
    ret.extend( expand("Offtargets/Bams/{samplename}.la.bam", samplename=OFFTARGETS.keys()) )
    print("get sample bams: ", ret)
    return ret    

def glob_bams(wildcards):
    """
    helper function; used to generate the files (with relative path) of all bam files (offtargets too); bams
    are per-lane.
    """
    if wildcards.samplename in SAMPLES:
        inner = SAMPLES[wildcards.samplename]
    elif wildcards.samplename in OFFTARGETS:
        inner = OFFTARGETS[wildcards.samplename]
    else:
        print("glob_bams error; unexpected sample_name", wildcards.samplename, file=sys.stderr)
        exit(1)
        
    files = expand("{dirname}/Bams/{samplename}_{sampleindex}_{lane}.bam",samplename=wildcards.samplename, \
        sampleindex= inner['sampleindex'], lane=inner['lanes'], dirname=inner['dirname'])
    print("glob bam, files: ", files)
    return files

def generate_files(wildcards):
        
        print("SAMPLES: ", SAMPLES)
        
        # files for on-target samples
        # per-lane bam files
        on= []
        for s in SAMPLES.keys():
            inner = SAMPLES[s]           
            #on.extend( expand("{dirname}/Bams/{samplename}_{sampleindex}_{lane}.bam", zip, \
            #dirname=s, samplename=s, sampleindex=inner["sampleindex"], lane=inner["lanes"]))

        # per sample, markdup, left align bams
        on.extend(expand("{dirname}/Bams/{samplename}.la.md.bam", zip, dirname=SAMPLES.keys(), samplename=list(SAMPLES.keys())))

        # and off-target samples
        off = []
        for s in OFFTARGETS.keys():
            inner = OFFTARGETS[s]           
            #off.extend( expand("{dirname}/Bams/{samplename}_{sampleindex}_{lane}.bam", zip, \
            #dirname="Offtargets", samplename=s, sampleindex=inner["sampleindex"], lane=inner["lanes"]))
            
        off.extend( expand("{dirname}/Bams/{samplename}.la.md.bam", dirname="Offtargets", samplename=list(OFFTARGETS.keys())) )
        
        on.extend(off)
        print("\n".join(on))
        return on
        
rule all:
    input:
        generate_files

        
rule bwa_mem_map:
    input:
        "Fastqs/{samplename}_{sampleindex}_{lane}_R1_001.fastq.gz",
        "SampleSheet.csv" # marker from bcl2fastq
    output:
        bam=temp("{dirname}/Bams/{samplename}_{sampleindex}_{lane}.bam") # note; no indexing required.
    wildcard_constraints:
        lane="L[0-9]+",
        sampleindex="S[0-9]+",
        #samplename="[A-Za-z0-9]+",
        #dirname="[A-Za-z0-9]+",
    log:
        "{dirname}/Logs/{samplename}_{sampleindex}_{lane}_bwa_mem.log"
    threads: config["bwamemParams"]["threads"]
    resources:
        mem_mb=config["samtoolsParams"]["mem_mb_sort"]
    run: # TODO: Extract Library Information!
        read1=input[0]
        read2=input[0].replace('_R1_','_R2_')
        
        if not os.path.isfile(read2):
            read2=""
            print("No read2 detected; assuming single-end sequencing", file=sys.stderr)
            
        # the set +o pipefail is necessary for the zcat | head (which generates a fail)
        shell(
			"""
			set +o pipefail
			header=$(zcat {input} | head -1)
			seqID=$(echo $header | cut -d':' -f 1 | cut -d' ' -f 1 | cut -c 2-)
			id=$(echo $header | cut -d':' -f 2-4 --output-delimiter='_')
			sm=$(echo {wildcards.samplename} | cut -d'_' -f 1)
			pu=$(echo $header | cut -d':' -f 3,4 --output-delimiter='_')
			rg=$(echo "@RG\tID:$seqID"_"$id\tSM:$sm\tLB:$sm\tPL:ILLUMINA")
			({config[bwamemParams][binary]} mem  -R $(echo "@RG\\tID:$sm"_"$seqID"_"$id\\tSM:$sm\\tLB:$sm\\tPL:ILLUMINA\\tPU:$pu") \
			-t {threads} \
			-M {config[refs][Reference]} \
			{read1} {read2} | \
			{config[samtoolsParams][binary]} sort -@ {config[samtoolsParams][samtoolsThreads]}  -m {config[samtoolsParams][sortMem]} - > {output.bam} ) 2> {log}
			"""
        )


# alpha version of the pipeline had two separate steps; merge, then left align.
# leftalign takes multiple input files, soo, let's combine those operations.
rule merge_and_leftalign_bams:
    input:
        glob_bams
    output:
        bam=temp("{dirname}/Bams/{samplename}.la.bam"),
        bai=temp("{dirname}/Bams/{samplename}.la.bai")
    log:
        "{dirname}/Logs/{samplename}.merge_leftalign.log"
    run:
        inputs_all= " ".join(["-I " + repr(i) for i in input])
        shell("{config[gatkParams][binary]} LeftAlignIndels --verbosity ERROR --QUIET true -R {config[refs][Reference]} {inputs_all} -O {output.bam} &> {log}")   

rule markdup_bams:
    input:
        "{dirname}/Bams/{samplename}.la.bam"
    output:
        "{dirname}/Bams/{samplename}.la.md.bam"
    threads: config["sambambaParams"]["threads"]
    log:
        "{dirname}/Logs/{samplename}.markdup.log"   
    shell:
        "{config[sambambaParams][binary]} markdup  --tmpdir=/tmp -l 9 -t {threads} {input} {output} 2> {log}"

