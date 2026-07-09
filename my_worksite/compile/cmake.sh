#!/bin/bash
#SBATCH --job-name=cmake
#SBATCH --cpus-per-task=48
#SBATCH --time=06:30:00
#SBATCH -o /data/mfulghieri/Nyx/outputs/compilation/CMake.out
#SBATCH -e /data/mfulghieri/Nyx/outputs/compilation/CMake.err


# Overwrite the path: from the $PATH command, exlude the anaconda prefixes at the left and copy only the remaining
export PATH="/opt/share/sw/gcc-11.3.0/slurm-22.05.7/bin:/opt/share/libs/gcc-11.3.0/ucx-1.13.1/bin:/opt/share/libs/gcc-11.3.0/pmix-4.2.2/bin:/opt/share/libs/gcc-11.3.0/hwloc-2.9.0/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Unset the conda variables form the path
unset CXX
unset GXX
unset CC
unset CPP_FOR_BUILD
unset READELF
unset CMAKE_PREFIX_PATH
unset LD_LIBRARY_PATH
unset LIBRARY_PATH
unset CXXFLAGS
unset CFLAGS
unset LDFLAGS
unset CPPFLAGS


# Load the hpc modules
module purge  
module load gcc-11.3.0/ompi-4.1.4_nccl
module load gcc-11.3.0/hdf5-1.14.1
module load gcc-11.3.0/ompi-4.1.4_nccl
module load cmake-3.22.1


cd /data/mfulghieri/Nyx/build

cmake -DNyx_MPI=no -DCMAKE_BUILD_TYPE=Release ..



