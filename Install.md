# Tapir installation overview

Installation involves three steps, which are as follows:

-   [**Download:**](#download) Tapir is first downloaded (github and zenodo); it uses many big files, I strongly recommend using a fast local disk (SSD/NVMe is recommended)
-  [**Configuration:**](#configure-host-os) Tapir also requires (somewhat minimal) user configuration. This involves making an isolated environment (Conda), installing a few niche softwares (parallel and fastqc), as well as modifying your bashrc.
  - Tapir can be configured to run using a conda environment
  - Tapir can also be run natively
    -  A handful of "extra" applications may not work (namely, Fastqc)
    -  And a handful more are assumed to already be installed (samtools, bcftools, bwa).
-  **Customization** <br>Tapir also supports performance tweaks; this includes references to local (fast) storage devices (NVME recommended; likely different than found in the `Install`.) as well as load-balancing (increasing/decreasing the number of threads used by various routines)

## Download

Tapir uses many large files. Because I'm ~~cheap~~ a fan of reproducible science, the largest of these files are stored on zenodo.

### Installation location

Genomic files are big; many genomic computations are IO bound (ie, the speed of the disk dictates the speed of the program). Tapir leverages fast (local) storage whereever possible. We recommend two *fast* local devices. The first (faster/larger) I would mount directly as /tmp, and the second (possibly slower) is wwhere I'd put Tapir. On (one of) our systems Tapir is installed on /mnt/ref/Tapir/. First, change directories to `/mnt/ref/` (adjusing the path as necessary to your system).

```
git clone https://github.com/Ahhgust/Tapir.git
cd Tapir
```
<br><br>

Download the relevant files from Zenodo:
```
# Imputation panel; this is big
wget https://zenodo.org/records/14171544/files/tapir_imputation.tar?download=1 && tar -xf tapir_imputation.tar && rm tapir_imputation.tar
# General genomic resources; including the reference genome and various liftOver files

# And last, source-code for relevant packages. Includes deepvariant (sif) and GATK (jars)

```
And place/mv ```Tapir``` to some centeral location and set the permissions
E.G., if you have a nice fast solid state drive mounted to /mnt/disk0/, try

```
cd ..
mv -i Tapir /mnt/disk0/
chmod -R 775 /mnt/disk0
```

### Configure Host OS

Tools like GATK and Sambamba make a ridiculous number of files. Likewise, it is nice to keep your command history for a good chunk of time and to handle concurrency. <br><br>
Modify your .bashrc to contain:

```
HISTCONTROL=ignoreboth

#append to the history file, don't overwrite it
shopt -s histappend

#for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=10000000
HISTFILESIZE=20000000
#added to keep the command history alive when there are multiple logins
#see: https://askubuntu.com/questions/80371/bash-history-handling-with-multiple-terminals/80882#80882
export PROMPT_COMMAND='history -a; history -r'

alias 'll=ls -tlr'

#increase the number of open files and the stack size
#the former is *needed* by sambamba/gatk
ulimit -n 65768
ulimit -s unlimited
```
(e.g., by editing ~/.bashrc, e.g., ```gedit ~/.bashrc``` and add the above to it)

<br><br>
If you're using some other shell, you'll need to translate these commands as necessary. Odds are, if you are using another shell (on purpose), then this won't be an issue for you.


### Conda install (recommended)

Conda is by far the easiest way to go about configuring Tapir.
*Skip this step if Conda is already installed*
```
cd ~/
mkdir -p src
cd src
curl -O https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh
bash ./Anaconda3-2024.06-1-Linux-x86_64.sh
```
and read (ha) and accept the licensing agreement; choose defaults for everything EXCEPT <br>
set conda to auto-load
```
You can undo this by running `conda init --reverse $SHELL`? [yes|no]
[no] >>> yes
```

Now, let's set up your environment:
```
conda env create --name tapir --file PATH_TO_TAPIR/snakemakes/tapir.yaml
```
where PATH_TO_TAPIR is the, well, the path to this program (`/mnt/disk0/Tapir` in the example).


### Native install

As of this writing, I still have yet to get Fastqc working outside of Conda (in my testing, the java versions conflict w/ GATK). I honestly see little value in Fastqc, but perhaps I'm wrong. <br>
Tapir is largely self-contained. It is assumed that the following are already in your environment:

- R
  - tidyverse (scripts by Tapir)
  - gplots (the following are required by GATK)
  - gsalib
  - reshape2
- Java
  - version 8 (ie, SDK 1.8)
- samtools
  - tested on version 1.9 or higher; older versions likely work too
- bcftools
  - version 1.15 or higher (tested on 1.15.1-20-g68fa3be)
- bwa
  - version 0.71 or higher (tested on 0.7.17-r1188)
- Fastqc
y- Pysam
- GNU parallel
- singularity (optional)
  - version 3.8.7-1.el7
  - only necessary for DeepVariant (which we don't, in practice, use)

<br><br>
Note, it is quite likely that any recent version of bwa, samtools and bcftools will work, and these are all very common tools.

