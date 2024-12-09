# Must be run from a DIRECTORY. Specifically the one created by the NovaSeq for a particular run
# it is assumed that this is run as ROOT by cron

nCores=128

if ! [ -z "$1" ]
then
    nCores=$1
    if ! [[ "$nCores" =~ ^[0-9]+$ ]]
    then
        echo "Optional argument $1 needs to be an integer (number of cores)."
        exit 1
    fi
fi

if [ `whoami` != "jarvis" ]
then
    echo `whoami`
    echo "This script needs to be run as jarvis..."
    exit 1
fi

if [ ! -f "CopyComplete.txt" ]; then
    echo "Copying not complete"
    exit 2
fi

if [ -d "Fastqs" ]; then
    echo "There is already a Fastq directory... nothing left to do."
   exit 3
fi

# count how many fasta files need to get written. Optimize concurrency therein"
# taken as number of lines after [Data] flag in the sample sheet (which also includes the column headers, so i starts as -1)
nOut=`cat *csv | perl -e '$i=-1; $j=0; while (<>){ if (m/^\[Data\]/) { $j=1; } elsif ($j) { $i++; } } print "$i\n";'`

# beware: we want $PWD to stay as a literal string, but we want interpolation (nCores nOut) elsewhere
# strong assumption: there is one *.csv (one sample sheet)
c='/usr/local/bin/bcl2fastq -R $PWD --minimum-trimmed-read-length 10 -o Fastqs -p '
c="$c $nCores --fastq-compression-level 9 --writing-threads $nOut"
c="$c"' --sample-sheet <(/eva/codebase/bin/BCLToFastqHokeyPokey.py < *.csv)'
c="$c"' &> Fastqs/bcl.outerr && /eva/codebase/bin/BCLToFastqPokeyHokey.py ./Fastqs'

#'/usr/local/bin/bcl2fastq -R $PWD -o Fastqs -p 60 --sample-sheet <(/usr/local/bin/BCLToFastqHokeyPokey.py < *.csv) --fastq-compression-level 9 --writing-threads 8  &> Fastqs/bcl.outerr && /usr/local/bin/BCLToFastqPokeyHokey.py ./Fastqs
# echo $c
# UNTESTED! (in the script at least this is untested)
mkdir -p Fastqs && chown "jarvis:jarvis" Fastqs && bash -c "$c"
