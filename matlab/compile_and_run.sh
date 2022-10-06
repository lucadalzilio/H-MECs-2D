#!/bin/sh

bkill -J hmec_t1

rm -rf *mat
rm -rf EVO*
rm -rf lsf*
rm -rf *txt

printf "0" >> file.txt

export OMP_NUM_THREADS=1

module load matlab/R2017b

bsub -W "24:00" -R 'rusage[mem=10000]' -J hmec_t1 matlab -nodisplay -nojvm -singleCompThread -r h_mec_RSF

export OMP_NUM_THREADS=1
bsub -W "24:00" -R 'rusage[mem=10000]' -J hmec_t1 -w "ended(hmec_t1)" matlab -nodisplay -nojvm -singleCompThread -r h_mec_RSF
bsub -W "24:00" -R 'rusage[mem=10000]' -J hmec_t1 -w "ended(hmec_t1)" matlab -nodisplay -nojvm -singleCompThread -r h_mec_RSF
bsub -W "24:00" -R 'rusage[mem=10000]' -J hmec_t1 -w "ended(hmec_t1)" matlab -nodisplay -nojvm -singleCompThread -r h_mec_RSF
bsub -W "24:00" -R 'rusage[mem=10000]' -J hmec_t1 -w "ended(hmec_t1)" matlab -nodisplay -nojvm -singleCompThread -r h_mec_RSF

rm -rf lsf*


