#!/usr/bin/env python3
# Written by August Woerner
#
# This is syntactic sugar, allowing jobs to be more easily run on a HPC environment (slurm)
#
# It writes slurm submit scripts for you.
# bonus points, it can set things up to run in parallel within a node (-p) and it supports job arrays (-a)
# (either or both).
# It will either takes command line arguments
# or it will take standard input
# for the latter, each line of input can be put into one of N job arrays (-a N)
# and within a job array index, all tasks can be run in series (-a N) or in parallel (-a N -p)
#
# TODO;
# Things left undone: sub-node allocation funniness.
#
# Usage:
#


import argparse
import os
import sys


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog=sys.argv[0], usage="\n\ncat aBunchOfUnixCommands | %(prog)s [options] | sbatch \n\nOR\n\n%(prog)s [options] 'some bash command' | sbatch")
    parser.add_argument('-m', '--module_load', dest='M', help="different modules to load. Takes multiple arguments (module load M)", nargs='+', default=[])
    parser.add_argument('-e', '--exclude', dest='E', help="Excludes nodes (via --exclude); -e g005 to exclude node 5", default="")
    parser.add_argument('-n', '--nice', dest='N', help="Sets the niceness; 0 means run it right away (max priority), 2147483645 means run it eventually (min priority)", default=0, type=int)
    parser.add_argument('-N', '--NICE', dest='NN', help="Sets niceness to it's maximum (2147483645; min priority)", default=0, type=int)
    parser.add_argument('-c', '--conda_activate', dest='C', default="", help="Activates the named conda environment")
    parser.add_argument('-l', '--logfile_name', dest='L', default="test_%j.log", help="The name of the logfile created")
    parser.add_argument('-p', '--parallel', dest='P', default=False, action="store_true", help="Only works with pipes (1st usage style). Runs tasks in parallel WITHIN a node. Note all jobs will run at the same time with no load balancing!")
    parser.add_argument('-a', '--array', dest='A', default=1, type=int, help="Only works with pipes (1st usage style). Makes a job array of size A; jobs run in parallel BETWEEN nodes. Note: -p AND -a can be used together")
    


    (results, xtra) = parser.parse_known_args(sys.argv[1:])


    if len(xtra) and xtra[0] == '--': # argparse bug. If you use -- (end of argument string), it gets included in the remaining bits
        xtra= xtra[1:]
    
    print("#!/bin/bash")
    print("#SBATCH --output=" , results.L, sep="")
    if results.E:
        print("#SBATCH --exclude=" , results.E, sep="")

    if results.N != 0:
        print("#SBATCH --nice=" , results.N, sep="")
    elif results.NN != 0:
        print("#SBATCH --nice=2147483645" , sep="")

        
    if results.A>1:
        if xtra:
            print("\n\nThat's an error. I only can make job arrays if you pass me standard input. \n\n", xtra, sep="\n", file=sys.stderr)
            exit(1)
            
        print("#SBATCH --array=0-", results.A-1, sep="")
    
    if len(results.M):
        for mod in results.M:
            print("module load " , mod, sep="")
    

    print()

    if len(results.C):
        print("conda activate" , results.C, sep=" ")
        print()
        
    
    print("cd $SLURM_SUBMIT_DIR")
    print()
    print()

    if xtra:
        print( " ".join(xtra) )

    elif results.A>1: # job arrays. This isn't the prettiest solution. However, it is the solution that works on a stream (not knowing it's length)
        i=0
        for line in sys.stdin:

            print("if [ ", i%results.A , " -eq  ${SLURM_ARRAY_TASK_ID} ]; then")
            
            if results.P:
                print(line.rstrip() , "&")
            else:
                print(line)

            print("fi")
            
            i += 1


        if results.P:
            print("wait")


        if i / results.A > 511:
            print("That's not ideal Each process has (at least) 512 concurrent subprocesses. That's too much!", file=sys.stderr)
            
            
    elif results.P:        
        i=0
        for line in sys.stdin:
            print(line.rstrip() , "&")
            i += 1


        print("wait")
        if i > 511:
            print("Not a good idea, yo! " , i , " commands is way to many. Keep it less than 512!", file=sys.stderr)
            
    else:
        for line in sys.stdin:
            print(line, end="")

        
    print()
    # close the door (I honestly don't know if this matters, but it seems like a good idea
    if len(results.C):
        print("conda deactivate" , results.C, sep=" ")

    
