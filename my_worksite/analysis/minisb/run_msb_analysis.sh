#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:30:00
#SBATCH -o /data/mfulghieri/Nyx/my_worksite/analysis/minisb/outputs/log.out
#SBATCH -e /data/mfulghieri/Nyx/my_worksite/analysis/minisb/outputs/err.err


python3 /data/mfulghieri/Nyx/my_worksite/analysis/minisb/analyze_msb.py \
       --plotfiles "/data/mfulghieri/Nyx/outputs/mini_sb/msb_13217" \
       --save "/data/mfulghieri/Nyx/my_worksite/analysis/minisb/outputs"





