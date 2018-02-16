#!/bin/bash
#SBATCH -D /Users/yangjl/Documents/Github/pvpDiallel
#SBATCH -o /home/jolyang/Documents/pvpDiallel/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/pvpDiallel/slurm-log/error-%j.txt
#SBATCH -J tw
set -e
set -u

GenSel4R slurm-scripts/tw_trun.inp > slurm-scripts/tw_trun.log
python /home/jolyang/bin/send_email.py -s slurm-scripts/tw_trun.inp
