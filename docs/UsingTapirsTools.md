# Using Tapir's components
Tapir provides a fair number of genomic tools, all prepackaged for (relative) ease of use.
Further, some of the tools have yet to be integrated into Tapir. Examples include:
-  Ibis
-  Plink (v1.9)
-  IGV (`igv.sh`)

We put any program that you can execute in the `bin` directory. To see all of the programs, do:

```
ls $TAPIR/bin
``` 

As you can see, we provide some nicities, like igv (`igv.sh`), that may be useful for advanced diagnostics.
It goes without saying, you can run any of these programs that you like.
For example, the below will ask GATK to print out all of the subtools that it provides.

```
$TAPIR/bin/gatk --list
```

For those that are truly keen on doing this routinely (and would like to NOT write `$TAPIR` a bunch of times, you can always add `$TAPIR/bin` to your PATH. 
Here is an [External link](https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path) that provides generic instructions for adding a directory to your PATH.
Doing so will let you just type:
```
gatk --list
```
for example.


