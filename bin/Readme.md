## Portable binaries/files used by Tapir

Herein are the binaries (x86_64) and scripts (minimal dependencies, perl and python3) used by Tapir.
I have done my best to ensure that the relevant binaries are statically compiled; they should be portable to any linux system that is sane
(sorry, 32-bit linux is not sane for genomics applications). There are three noteworthy exceptions:

### bcl-convert
Note, bcl-convert has not been officially added to Tapir.
bcl-convert has dynamically linked libraries. They are minimal (Try `ldd bin/bcl-convert` to list them), and are likely to be portable
but the onus is on the user to ensure portability (eg, install the relevant dependencies and/or make links to the relevant paths, or reinstall
this program from Illumina.

### atlas
atlas has numerous dependencies. You can try updating the various libraries, but a better answer is to compile this from scratch. Or see:
https://github.com/yassineS/ATLAS-Docker

### fastqc
fastqc must be installed outside of Tapir. Tapir will run in the absence of fastqc, but it is worth mentioning that here. IMO, fastqc is fairly resource intensive
and I have yet to see value in it, however maybe the use-case just hasn't presented itself yet...


