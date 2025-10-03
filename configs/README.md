# Custom configuration
Many tools by Tapir take command line arguments; some of the arguments are specified in the configuration file.
Two configuration files are provided; the default (`config_v_2_standard.yaml`) and one that is more permissive in its genotyping (`config_v_2_ultralow.yaml`)
## Config format
*Do not modify config files in place!*
You'll want to make a new config file instead. It's not terribly easy to read config files. I have (loosely) standardized the following:
-  *binary*
   - Refers to the executable file/script that is invoked
-  *threads*
   - Refers to the degree of multithreading.
-  Notes on paths
   -  Any file that is refered to in plain text is from the environment (Mamba or native, likely)
      -  eg, `bcftools` is assumed to be in your PATH
   -  Any file that has a directory in front of it (`bin/samStats.py`) is found local to Tapir 
      -  eg, `bin/samStats.py` is equivalent to saying `$TAPIR/bin/samStats.py`   
-  Notes on Scripts
   - .R, .py and .pl extensions are run with `Rscript`, `python3` and `perl` prepended, respectively
      -  (where each script interpreter, eg `python3` is from your environment)


Note, many more tools are supported in the config that are (currently) used by Tapir.
(so if you see something that refers to an unsupported function, that's fine).

## Modifying the config
Make a copy of a config file and change it as you see fit. <br>
Please reach out to me (august.woerner AT unthsc.edu) for help. <br>
Note: no config is optimal, though I do try and make them sensical. Use at your own risk <br>
I've included several modified config files
- `config_v_2_ultralow.yaml`
   - Final GLIMPSE genotypes:  posterior probability from 0.99 to 0.50; Bayes Factor 100 <br>
- `config_v_2_low_mem.yaml`
   - Parameters adjusted to run on a machine with 16Gb of RAM
- `config_v_2_standard_degraded.yaml`
   - Sets alignment parameters from: doi.org/10.1002/ece3.7056 (-k 19 -r 2.5) (affects bcl2bams)
   - Performs hard filtering (removes readS): samtools view -e 'sclen/qlen<0.20 && qlen > 39' (no more than 20% soft clipping; at least 40 bases aligned)


## Using a modified config
Example, dry run (-n):
```
snakemake -c 40 -n -s $TAPIR/snakemakes/bams2genotypes.smk --until call_glimpse2_and_bcftools --configfile $TAPIR/configs/config_v2_ultralow.yaml
```

