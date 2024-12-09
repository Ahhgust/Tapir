## Snakes on Tapirs!

This directory includes snakemake routines, a python script (GenomixHelper.py) to make quality of life improvements to snakemake, as well as a config file.

<br>
For those interested in modifying Tapir...

### Config files
Notes on paths for binaries/executables:
- Relative paths (bin/deML) in the config file are interpreted as relative paths, relative to the parent directory of this file. If you wish to add a new script/binary, add it to bin/ and use this path-style
- Path-less binaries are inherited from the environment. This includes fastqc (niche, cumbersome to install) and samtools (very generic).
- Absolute paths (/home/augustw/bin/) are treated as-is.

