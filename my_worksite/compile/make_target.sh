#!/bin/bash
#SBATCH --job-name=make
#SBATCH --cpus-per-task=48
#SBATCH --time=06:30:00
#SBATCH -o /data/mfulghieri/Nyx/outputs/compilation/make.out
#SBATCH -e /data/mfulghieri/Nyx/outputs/compilation/make.err

cd /data/mfulghieri/Nyx/build/Exec/MiniSB

make -j 48