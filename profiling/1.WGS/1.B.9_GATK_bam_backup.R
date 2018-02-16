### Jinliang Yang
### June 7th, 2016


library(farmeR)

files <- list.files(path="largedata/gatk_vcf", pattern="RG.bam$")

# /iplant/home/yangjl/methyl_teo20/bam_dedup
run_iput(files, jobs1cpu=18, localpwd="largedata/gatk_vcf")

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemh --mem 8196 --ntasks=1 /home/jolyang/Documents/Github/methylation/slurm-script/run_iput_array.sh


### extract and gz fastq files
# for i in *.lz4; do lz4 -d $i; done
# for i in *.fastq; do gzip --fast $i; done
# for i in *.lz4; do rm $i; done

