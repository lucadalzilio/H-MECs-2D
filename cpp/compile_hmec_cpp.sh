#!/bin/sh

CPU_TYPE="XeonGold_6150"
MAIN_FILE="h_mec_rsf_v9.cpp"
LINK_FILES="hdf5.cpp"
declare -i NUM_THREADS=4
declare -i START_TIMESTEP=0

echo
echo "=============================================================================================="
echo "Submitting ${MAIN_FILE} using $NUM_THREADS cores of type ${CPU_TYPE} and starting from timestep $START_TIMESTEP"
echo
echo "Terminate running job..."
echo
bkill -J compile_hmec
bkill -J hmec_cpp

echo
echo "Removing old files..."
rm -rf EVO*
rm -rf lsf*
rm -rf *txt
if [ $START_TIMESTEP -eq 0 ]
then
	echo "Starting new simulation -> removing old .h5 files..."
	rm -rf *.h5
fi

printf "$START_TIMESTEP" >> file.txt

echo
echo "Set max threads to 1 for compilation"
export OMP_NUM_THREADS=1

echo
echo "Loading Intel Compiler, Eigen Library, HDF5 library and libszip module..."
module load intel/19.1.0 eigen/3.3.9 hdf5/1.10.7 libszip/2.1.1

echo "Compiling Program..."
echo

if [ -z "$CPU_TYPE" ]
then
	bsub -W "4:00" -R "rusage[mem=50000]" -J compile_hmec icpc -std=c++17 -O3 -I/usr/include/eigen3 -I/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/include -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -ftree-vectorize -march=corei7-avx -mavx -DEIGEN_USE_MKL_ALL -DMKL_LP64 -I/intel/mkl/include -qopenmp -qoverride-limits $MAIN_FILE $LINK_FILES -L/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_hl_cpp.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_cpp.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_hl.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5.a -lhdf5 -lsz -lz -lrt -Wl,-rpath -Wl,/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib -o cpp_pardiso.out	
else
	bsub -W "4:00" -R "rusage[mem=50000]" -R "select[model=$CPU_TYPE]" -J compile_hmec icpc -std=c++17 -O3 -I/usr/include/eigen3 -I/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/include -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE -ftree-vectorize -march=corei7-avx -mavx -DEIGEN_USE_MKL_ALL -DMKL_LP64 -I/intel/mkl/include -qopenmp -qoverride-limits $MAIN_FILE $LINK_FILES -L/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_hl_cpp.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_cpp.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5_hl.a /cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib/libhdf5.a -lhdf5 -lsz -lz -lrt -Wl,-rpath -Wl,/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib -o cpp_pardiso.out
fi

echo
echo "Set max threads to $NUM_THREADS"
export OMP_NUM_THREADS=${NUM_THREADS}

echo
echo "Submitting Job..."
echo

if [ -z "$CPU_TYPE" ]
then
	echo "No CPU Type specified..."
	bsub -W "4:00" -R "rusage[mem=10000]" -J hmec_cpp -w "ended(compile_hmec)" -o output.txt -n $NUM_THREADS ./cpp_pardiso.out
else
	echo "Only using cpu type: $CPU_TYPE"
	bsub -W "4:00" -R "rusage[mem=10000]" -R "select[model=$CPU_TYPE]" -J hmec_cpp -w "ended(compile_hmec)" -o output.txt -n $NUM_THREADS ./cpp_pardiso.out
fi

echo "=============================================================================================="
