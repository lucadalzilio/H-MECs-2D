#!/bin/sh

bkill -J hmec_cpp

rm -rf EVO*
rm -rf lsf*
rm -rf *txt

printf "0" >> file.txt

export OMP_NUM_THREADS=1

module load gcc/6.3.0 eigen/3.3.9

g++ -std=gnu++11 -O3 h_mec_rsf_v2.cpp -o h_mec_cpp

bsub -W "24:00" -R 'rusage[mem=10000]' -J hmec_cpp -o output.txt ./h_mec_cpp

