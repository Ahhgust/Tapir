# Custom configuration
Many tools by Tapir take command line arguments; some of the arguments are specified in the configuration file.
Two configuration files are provided; the default (`config_v_2_standard.yaml`) and one that is more permissive in its genotyping (`config_v_2_ultralow.yaml`)
## Config format
*Do not modify config files in place!*
You'll want to make a new config file instead. It's not terribly easy to read th config files. I have (loosely) standardized the following:
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
I've included one modified config file (`config_v_2_ultralow.yaml`), with a change in the second-to-last line. In it, we dropped the posterior probability from 0.99 to 0.95. <br>
This lets us emit more sites (increasing call rate) at the expense of accuracy.

## Using a modified config
Example, dry run (-n):
```
snakemake -c 40 -n -s $TAPIR/snakemakes/bams2genotypes.smk --until call_glimpse2_and_bcftools --configfile $TAPIR/configs/config_v2_ultralow.yaml
```

