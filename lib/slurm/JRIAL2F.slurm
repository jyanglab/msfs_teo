#!/bin/sh
#SBATCH --job-name=JRIAL2F.map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=160:00:00
#SBATCH --mem=60G
#SBATCH --error=JRIAL2F.err
#SBATCH --output=JRIAL2F.out
#SBATCH --licenses=common
#SBATCH --partition=jclark
module load bowtie
/home/jyanglab/gxu6/software/Bismark/bismark --bowtie2 -p 16 -N 1 -L 25 -o /lustre/work/jyanglab/shared/sh_gen /common/jyanglab/gxu6/msfs_teo/Basic_data/TEO_PseudoRef/JRIAL2F -1 /lustre/work/jyanglab/gxu6/msfs_teo/Basic_data/WGBS/teo/JRIAL2F_R1.fastq.gz -2 /lustre/work/jyanglab/gxu6/msfs_teo/Basic_data/WGBS/teo/JRIAL2F_R2.fastq.gz > JRIAL2F.log