### Jinliang Yang
### 10-16-2016
### compute IBD among 20 teo lines


library(farmeR)
### split into chrs
cmd1 <- "cd largedata/gatk_vcf"
#cmd2 <- "bcftools index -f JRI20_bi_snps_annot.vcf.gz"
cmd4 <- c()
for(i in 1:9){
    tmp <- paste0("bcftools filter -r ", i,  " JRI20_bi_snps_annot.vcf.gz",
                  " -Ou -o ", "vcf_bychr/chr", i, "_JRI20_bi_snps_annot.vcf.gz")
    cmd4 <- c(cmd4, tmp)
}

set_farm_job(slurmsh = "slurm-script/bcf2plink.sh",
             shcode = c(cmd1, cmd4), wd = NULL, jobid = "split-chr",
             email = "yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", "8"))
###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 64G --time 2:00:00 --ntasks=8 slurm-script/bcf2plink.sh


#########################
library(farmeR)


cmd <- c("module load java", "cd largedata/gatk_vcf/vcf_bychr")
for(i in 1:9){
    cmd1 <- paste0("bcftools convert chr", i, "_JRI20_bi_snps_annot.vcf.gz -Ov -o chr", i, "_JRI20_bi_snps_annot.vcf")
    cmd2 <-  paste0("java -Xmx60g -jar ~/bin/beagle.22Feb16.8ef.jar nthreads=8 ",
                   "gt=chr", i, "_JRI20_bi_snps_annot.vcf ibd=true out=chr", i, "_ibd")
    cmd <- c(cmd, cmd1, cmd2)
}


set_farm_job(slurmsh = "slurm-script/run_beagleibd.sh",
             shcode = cmd, wd = NULL, jobid = "ibd",
             email = "yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", "8", "60G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 60G --ntasks=8 --time=96:00:00 slurm-script/run_beagleibd.sh

#No genetic map is specified: using 1 cM = 1 Mb
ibd <- read.table("largedata/gatk_vcf/JRI20_ibd.ibd")
write.table(ibd, "largedata/gatk_vcf/JRI20_ibd_partial.ibd", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

