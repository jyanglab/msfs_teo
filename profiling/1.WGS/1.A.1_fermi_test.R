### Jinliang Yang
### use Li Heng's De novo assembly based variant calling pipeline
### testing

library(farmeR)
fq1 <- list.files(path="~/dbcenter/BMfastq", pattern="_1.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="~/dbcenter/BMfastq", pattern="_2.fastq.gz$", full.names = TRUE)
sampleid <- read.table("~/dbcenter/BMfastq/sampleid.txt", header=TRUE)

fq <- data.frame(fq1=fq1[c(1,17)], fq2=fq2[c(1,17)], out=gsub(".sra.*", "", fq1[c(1,17)]),
                 sample=c("Mo17", "B73"))
fq$out <- gsub("BMfastq/", "BMfastq/fermi/", fq$out)


run_fermikit(fq,kitpath="$HOME/bin/fermikit/fermi.kit",
             ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
             s='2.3g', l=75,
             email="yangjl0930@gmail.com", 
             runinfo = c(TRUE, "bigmemh", 16))


###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemh --mem=128000 --ntasks=16 slurm-script/run_fermikit_array.sh

run_fermikit_vcfcall(bamdir="$HOME/dbcenter/BMfastq/fermi", kitpath="$HOME/bin/fermikit/fermi.kit",
                     ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
                     email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemh", 1))
