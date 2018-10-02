#!/bin/sh
#SBATCH --job-name=1147-8D.map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=160:00:00
#SBATCH --mem=60G
#SBATCH --error=1147-8D.err
#SBATCH --output=1147-8D.out
#SBATCH --licenses=common
#SBATCH --partition=jclark
module load bowtie
/home/jyanglab/gxu6/software/Bismark/bismark --bowtie2 -p 16 -N 1 -L 25 -o /lustre/work/jyanglab/shared/sh_gen /common/jyanglab/gxu6/B73_v4 -1 /common/jyanglab/gxu6/msfs_teo/Basic_data/WGBS/landrace/1147-8D_1.fq.gz -2 /common/jyanglab/gxu6/msfs_teo/Basic_data/WGBS/landrace/1147-8D_2.fq.gz > 1147-8D.log