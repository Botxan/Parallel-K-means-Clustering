#!/bin/bash

# Paths to databases
dbgen=~/ARC/Genetics/dbgen.dat
dbdise=~/ARC/Genetics/dbdise.dat

if [[ $1 == "s" ]];
then # serial execution
    clear
    if [[ $# -eq 1 ]]; # with all the elements
    then
        echo "[*] Running serial program with all the elements [*]"
        echo "~/genetics/serial/gengroups_s $dbgen $dbdise"
        ~/genetics/serial/gengroups_s $dbgen $dbdise
    elif [[ $# -eq 2 ]] && [[ $2 -eq 1000 ]]; # for 1000 elements
    then
        echo "[*] Running serial program with 1000 elements [*]"
        echo "~/genetics/serial/gengroups_s $dbgen $dbdise 1000"
        ~/genetics/serial/gengroups_s $dbgen $dbdise 1000 2>&1 |& tee >(tail -14 >~/genetics/serial/memusage_s.out)
    else
        echo "Use: run s [1000]"
    fi
elif [[ $1 == "p" ]] && [[ $2 =~ ^(1|2|4|8|16|32|64|128)$ ]]; # parallel execution
then
    # Run script with source run.sh in order to make the variable change persistent
    export OMP_NUM_THREADS=$2
    clear
    if [[ $# -eq 2 ]]; # with all the elements
    then
        echo "[*] Running parallel program with $2 threads and all the elements [*]"
        echo "~/genetics/parallel/gengroups_p $dbgen $dbdise"
        ~/genetics/parallel/gengroups_p $dbgen $dbdise |& tee >(tail -14 >~/genetics/serial/memusage_s.out)
    elif [[ $# -eq 3 ]] && [[ $3 -eq 1000 ]]; # with all the elements
    then
        echo "[*] Running parallel program with $2 threads and 1000 elements [*]"
        echo "~/genetics/parallel/gengroups_p $dbgen $dbdise 1000"
        ~/genetics/parallel/gengroups_p $dbgen $dbdise 1000 |& tee >(tail -14 >~/genetics/serial/memusage_s.out)
    else
        echo "Use: run p numthreads (1|2|4|8|16|32|64|128) [1000]"
    fi
else
    echo "Use: run (s (for serial) [1000] | p (for parallel) numthreads (1|2|4|8|16|32|64|128) [1000])"
fi