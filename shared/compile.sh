#!/bin/bash

# Check compile mode is correct
if [[ $# -eq 1 ]];
then
    clear
    if [[ $1 == "s" ]];
    then
        echo "[*] Compiling serial program [*]"
        echo "gcc -O2 -lm -o ~/genetics/serial/gengroups_s ~/genetics/serial/gengroups_s.c ~/genetics/serial/fungg_s.c"
        gcc -O2 -o ~/genetics/serial/gengroups_s ~/genetics/serial/gengroups_s.c ~/genetics/serial/fungg_s.c -lm
    elif [[ $1 == "p" ]];
    then
        echo "[*] Compiling parallel program [*]"
        echo "gcc -O2 -fopenmp -lm -o ~/genetics/parallel/gengroups_p ~/genetics/parallel/gengroups_p.c ~/genetics/parallel/fungg_p.c"
        gcc -O2 -fopenmp -lm -o ~/genetics/parallel/gengroups_p ~/genetics/parallel/gengroups_p.c ~/genetics/parallel/fungg_p.c
    else
        echo "Invalid compile mode $1"
    fi
else
    echo "Use: compile.sh (s (for serial) | p (for parallel))"
fi