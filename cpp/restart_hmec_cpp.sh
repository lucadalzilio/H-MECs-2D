#!/bin/sh

CPU_TYPE="XeonGold_6150"
declare -i NUM_THREADS=4
declare -i START_TIMESTEP=90

echo
echo "=============================================================================================="
echo "Restarting ${MAIN_FILE} using $NUM_THREADS cores of type ${CPU_TYPE} and starting from timestep $START_TIMESTEP"
echo
echo "Terminate running job..."
echo
bkill -J restart_hmec

echo
echo "Remove old lsf files"
rm -rf lsf*
rm output.txt
rm file.txt

printf "$START_TIMESTEP" >> file.txt

echo
echo "Loading Intel Compiler, Eigen Library, HDF5 library and libszip module..."
module load intel/19.1.0 eigen/3.3.9 hdf5/1.10.7 libszip/2.1.1

echo
echo "Set max threads to $NUM_THREADS"
export OMP_NUM_THREADS=${NUM_THREADS}

echo
echo "Submitting Job..."
echo

if [ -z "$CPU_TYPE" ]
then
	echo "No CPU Type specified..."
	bsub -W "4:00" -R "rusage[mem=10000]" -J restart_hmec -o output.txt -n $NUM_THREADS ./cpp_pardiso.out
else
	echo "Only using cpu type: $CPU_TYPE"
	bsub -W "4:00" -R "rusage[mem=10000]" -R "select[model=$CPU_TYPE]" -J restart_hmec -o output.txt -n $NUM_THREADS ./cpp_pardiso.out
fi

echo "=============================================================================================="
