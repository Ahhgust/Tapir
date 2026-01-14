# Virtual Tapirs!

Tapir is built for a Unix environment. That doesn't mean you can't run it on Windows, however, you just need to run Unix (on Windows).

## Setting Up Virtual Box

Virtual box is a free, open source virtual machine. To unix Tapir on your Windows (or Linux) machine first:

### Minimum requirements
We recommend having
-  at least 16gb of RAM (an ABSOLUTE minimum)
-  200Gb of hard drive space
   -  and some form of solid state drive (NVME is best; traditional SSDs are fine too).
   
I was able to successfully run Tapir on my mid-grade laptop (where the minima above came from). <br>
In general, the cores you throw at some workflow, the more memory you use. If Tapir runs out of memory, use less cores.
I also made a "low memory" config (`$TAPIR/configs/config_v_2_low_mem.yaml`). This was suitable for my laptop.
e.g.,
try:
```
snakemake -c 3 (...) --configfile $TAPIR/configs/config_v_2_low_mem.yaml
```
where (...) specifies what workflow you want to run , and `-c 3` requests three CPUs.
It goes without saying that while Tapir *can* run on cheapo hardware, if you use it a fair amount, it would be wise to invest in something nicer. In general, you want more CPUs (as many as possible), more memory (~0.5Tb), and fast IO (NVMEs).
For reference, on our ~3.5 year old servers at the CHI, a miseq run can be Tapir'd in about an hour.


### Install Virtual Box
This will let you use Virtual Machines on your computer..

-  Download the Virtual Box installer
   -  [External Link](https://www.virtualbox.org/wiki/Downloads)
-  Click the appropriate OS (probably Windows or Mac)
   -  Download and Install it
      -  sign the EULA and agree to everything
      -  stick to the defaults

### Configure your Windows 11 (or Unix?) machine

Tools such as `bcl-convert` and `GLIMPSE` use special machine instructions (AVX); these instructions may not be enabled by default. AVX instructions also means that Tapir won't work on a Mac (unless you're rocking an old Intel Mac)><br>
`Control Panel` -> `Programs` -> `Turn Windows features on or off` (Under Programs and Features) <br><br>
And make sure that the `Virtual Machine Platform` and `Windows Hypervisor Platform` are turned off.


### Download the Tapir image
I would let this run overnight; the image is big, and because the powers that be make sharing large files exceptionally hard, I'm hosting this file on one of my own dime...
Click [here](https://www.dropbox.com/scl/fi/6rcdh0wttdctm894h9zx1/TapirOnUbuntu.ova?rlkey=5v29kk1ks1ieo8ow9tbwa2mhs&st=04p52yti&dl=1) (Note: dropbox link) to begin the download.

### Checking the Download
Windows says that the download was successful. Was it? When working with very large files, the answer can be "no" more often than you would think. <br>
In Windows, type:
```
certutil -hashfile "TapirOnUbuntu.ova" MD5
```
To compute a md5sum; compare it to the one in `MD5s.txt` (e.g., `cat $TAPIR/MD5s.txt`) <br>
(and note that the VM image needs to be in the same directory for the above command to work).
If the two hashsums don't agree, something went wrong with the download and you should do it again.


### Import the Virtual box image

-  Import the VM
   - File->Import Appliance
-  Configure the VM
   -  Adjust
      - The number of CPUs
      - The amount of RAM 
		- 12Gb is a *minimum* (note this is the amount of memory dedicated to your VM. The host OS will need some too).
		- The more, the better.
      - And the disk-space
		- Aim for at least 200Gb
   -  To match your machine.

### Notes
On my machine, I get errors (below) when I import the image. I honestly don't know why, but if you just open (green plus sign) the VM, it works just fine.
Could not create the imported medium (...) (VERR_VD_VMDK_INVALID_FORMAT).

### Logging in

When you use Tapir, log in as: `tapir_user`.
```
User:
tapir_user
Password:
4tapir!
```

`tapir_user` has administrative privileges and `Tapir` is installed locally to their account. 
`vboxuser` retains the default password:  `changeme`.  Feel free to, well, change it.


### Sanity test
To make sure everything went "okay" with the Tapir VM installation, let's first try a sanity check/text, just to make sure nothing is crazy. <br>
When logged in as tapir_user, try:
```
~/src/Tapir/bin/bcl-convert -h
```

If you get a nice help-message, that's a good thing; if it says "Illegal instruction" (and a fair amount of other text), 
that means you need to adjust how your host OS interacts with the Virtual Machine. See **Configure your Windows 11 machine** above.


### Using tapir

When you log in as `tapir_user` the Tapir mamba environment is preloaded. 
Tapir is installed to `/home/tapir_user/src/Tapir` (i.e., that is the value of `$TAPIR`).
And remember, if you're using a lower-end machine, use the low-memory config file (tuned as appropriate).

## Practical Example

We provide a virtual machine (VM) of Tapir. In brief, a VM let's you run a computer (operating system, OS) as software.
The OS in question is a 64-bit linux distribution with Tapir installed on it. For this to be useful, you want data on your machine (the "host" machine)
to be readable by the VM (`In` directory) and vice versa (`Out` directory). <br>
Note that Tapir creates/uses symbolic links (symlinks). Symlinks are supported on Windows, but they take some doing. 
In practical terms, this means you cannot have Tapir *write* to `Out` directly unless Windows has been configured appropriately. <br>
A lazier solution is to just copy the requisite files to `Out` instead (which is what we'll be doing in this walkthrough).

### Install Guest Additions:
Install Guest Additions. 
<br>
Try: (from the running VM):
-  Devices
-  "Insert Guest Additions CD Image"

That (for me) was enough to enable Guest Additions (which lets you share folders).
See additional instructions here: 
[External Link](https://docs.oracle.com/en/virtualization/virtualbox/7.2/user/guestadditions.html#guestadd-install)
<br>
and read about folder sharing (somewhat incomplete instructions from Oracle are found here) [External Link](https://docs.oracle.com/en/virtualization/virtualbox/7.2/user/guestadditions.html#guestadd-intro)


### Folder sharing

Make a directory on your computer called `Shared` <br>
<br><br>

In Oracle VirtualBox Manager, click:

-  Select the VM (TapirOnUbuntu); top left
-  Click on Settings (orange gear, also top left)
-  Click on Shared Folders
-  And Global Folders
   -  And the Add Folder icon

You will be asked to select a Folder Path. Navigate to the `Shared` directory you made earlier. <br>

Keep the `Folder Name` as `Shared` <br>
the Mount Point as `/media`  <br>
And put checkmarks beside: Auto-mount,  Make Global (which also selects Make Machine-permanent) <br>
Select `OK` (and OK again to exit out of the Setting menu) <br>
And restart your VM (if it is on, go to File->Close->Send the shutdown signal)
<br>

On your host machine (likely, the Windows machine), navigate to the `Shared` directory (In File Explorer)
You'll find a directory here called `vboxuser`. 
Inside *that* directory, make an `In` (bcl/fastq data go here) directory and an `Out` directory (Tapir produces data here). 

### Example

Let's try processing the FASTQ files from the `examples/` directory in Tapir. 
Please place the entire `fastq_example` directory in the `In` directory on the Host machine.<br>
Boot up the VM, and log in as `tapir_user`. The `In/` directory (in your home) is the same as the `In/` directory on the Host machine.

 From Guest (Ubuntu), type `ls -1 ~In/Fastqs`. It should return:

```
GMWOF5428715_S89_L001_R2_001.fastq.gz
MadeupSamplesheet_ForExample.csv
GMWOF5428715_S89_L001_R1_001.fastq.gz
```

If it's blank, or you have some kind of error. <br>
<br>
Note that Tapir (much) prefers BCLs, not FASTQs..., but FASTQs make much easier examples to work with.
To perform Step1 of Tapir, do:
```
nohup snakemake --keep-going  -c 4 -s $TAPIR/snakemakes/bcl2bam.smk --config Fastqdir=In/fastq_example \
	Samplesheet=In/fastq_example/MadeupSamplesheet_ForExample.csv \
	Experiment=FauxTest  \
	--configfile $TAPIR/configs/config_v_2_low_mem.yaml &
```
Note in the above, if your system supports making symlinks, you can `cd` to the `Out\` directory (or give it's absolute path under `Experiment` <br>
This will let all files that Tapir makes to be visible to the host.
<br>

Which in words will: 
-  Run bcl2bam 
   -  Using 4 cores/cpus (`-c 4`; ~20G of memory recommended)
      -  Adjust +/- as needed for your system
   -  If any program fails, `--keep-going`
-  Read from the In/ directory (host machine)
-  Write to the FauxTest directory (VM)
-  Using the "low memory" config file from Tapir.

Note that if this runs out of memory, Tapir simply dies. A hallmark is an empty (or nonexistent) bam file.
If this happens, try reducing -c  (`-c 3` perhaps?)
<br>
To run Step2, cd into the single sample's data directory:
```
cd FauxTest/Sample_Data/GMWOF5428715
```
And check to see that the right BAMs have been made:

```
ls fastq_example/Bams/
```
Which should have a handful of bams; one should be called `GMWOF5428715.la.md.bqsr.bam`
<br>
In words gives the sampleid (GMWOF...), which has been left-aligned `la`, duplicates have been marked `md`, and base quality score recalibration has been conducted `bqsr`. 
To run Step2 of Tapir, do:
```
snakemake -c 4 -s $TAPIR/snakemakes/bams2genotypes.smk  --configfile $TAPIR/configs/config_v_2_low_mem.yaml  --until call_glimpse2
```

Which runs the GLIMPSE2 pipeline on this particular sample.
Once this completes (in maybe an hour?), you can export the final reports files to the `Out` directory

```
cp -rf Final_Reports ~/Out
```
The `Final_Reports` directory contains various summaries of your data.

As well as the GEDmatch uploads

```
cp -rf Uploads ~/Out
```
This contains the genotype calls, compressed by gzip. Note you likely want to upload the `snps` file (less markers, but a higher proportion of informative markers)













