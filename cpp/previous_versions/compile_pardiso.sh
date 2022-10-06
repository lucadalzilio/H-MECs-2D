#!/bin/sh

CPU_TYPE="XeonGold_6150"
MAIN_FILE="h_mec_rsf_v4.cpp"

bkill -J hmec_cpp

rm -rf EVO*
rm -rf lsf*
rm -rf *txt

printf "0" >> file.txt

export OMP_NUM_THREADS=1

echo "Load Intel Compiler and Eigen Library..."
module load intel/19.1.0 eigen/3.3.9

echo "Compile Program..."
icpc -std=c++17 -O3 -I/usr/include/eigen3 -DEIGEN_USE_MKL_ALL -DMKL_LP64 -I/intel/mkl/include $MAIN_FILE -L/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -o h_mec_cpp_pardiso

echo "Submit Job..."
if [ -z "$CPU_TYPE" ]
then
	echo "No CPU Type specified..."
	bsub -W "24:00" -R "rusage[mem=10000]" -J hmec_cpp -o output.txt ./h_mec_cpp_pardiso
else
	echo "Only using $CPU_TYPE"
	bsub -W "24:00" -R "rusage[mem=10000]" -R "select[model=$CPU_TYPE]" -J hmec_cpp -o output.txt ./h_mec_cpp_pardiso
fi
