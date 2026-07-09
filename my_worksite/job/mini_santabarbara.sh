#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=64G
#SBATCH --time=12:00:00 
#SBATCH -J Mini_SB
#SBATCH -p long
#SBATCH -o /data/mfulghieri/Nyx/outputs/mini_sb/logs/job_%j.out
#SBATCH -e /data/mfulghieri/Nyx/outputs/mini_sb/logs/job_%j.err
#SBATCH --mail-type=END 
#SBATCH --mail-user=m.fulghieri@campus.unimib.it

# Useful paths
BASE_DIR="/data/mfulghieri/Nyx/build/Exec/MiniSB/"
OUT_DIR="/data/mfulghieri/Nyx/outputs/mini_sb/msb_${SLURM_JOB_ID}" 

# Make (after having used CMake)
cd /data/mfulghieri/Nyx/build/Exec/MiniSB
make -j 48

# Output dir
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# Symbolic link for the executable and the inputs in the output dir
ln -sf "${BASE_DIR}/nyx_MiniSB" .
ln -sf "${BASE_DIR}/inputs.32" .
ln -sf "${BASE_DIR}/ic_sb_32.ascii" .


# Run simulation: [./exe] [input] [overwrite input params from CL]
./nyx_MiniSB ./inputs.32 amr.max_step=5