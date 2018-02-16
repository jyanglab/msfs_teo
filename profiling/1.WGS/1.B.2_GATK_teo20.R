### Jinliang Yang
### GATK best practice to call SNP and InDels for teo20
### 4/6/2016


bamfile <- list.files(path="/group/jrigrp4/teosinte-parents/20parents/bams",
                      pattern="sorted.*bam$", full.names=TRUE)
inputdf <- data.frame(
    bam=bamfile, 
    out=gsub(".*20parents/bams","/home/jolyang/Documents/Github/methylation/largedata/gatk_vcf", bamfile), 
    group="1", 
    sample=gsub(".*sorted.|_index.*", "", bamfile),
    PL="illumina", LB="lib1", PU="unit1")
inputdf$out <- gsub("sorted.|_index.*", "", inputdf$out)

###########
library(farmeR)
run_GATK(inputdf[1:20,], 
         ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
         gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
         picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
         minscore = 5, markDup=TRUE, addRG=TRUE, 
         realignInDels=FALSE, indels.vcf="indels.vcf",
         recalBases=FALSE, dbsnp.vcf="dbsnp.vcf", 
         email="yangjl0930@gmail.com",
         runinfo = c(FALSE, "bigmemm", 16))

###############
gvcf <- list.files(path="largedata/gatk_vcf", pattern="RG.bam$", full.names=TRUE)


bamfile <- list.files(path="/home/jolyang/Documents/Github/methylation/largedata/gatk_vcf",
                      pattern="RG.bam$", full.names=TRUE)
inputdf <- data.frame(
    bam=bamfile, 
    out=bamfile, 
    group="1", 
    sample=gsub(".*/|.sorted.*", "", bamfile),
    PL="illumina", LB="lib1", PU="unit1")
inputdf$out <- gsub(".sorted.*", "", inputdf$out)

###########
library(farmeR)
run_GATK(inputdf[-1,], 
         ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
         gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
         picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
         minscore = 5, markDup=FALSE, addRG=FALSE, 
         realignInDels=FALSE, indels.vcf="indels.vcf",
         recalBases=FALSE, dbsnp.vcf="dbsnp.vcf", 
         email="yangjl0930@gmail.com",
         runinfo = c(FALSE, "bigmemm", 16))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 140000 --ntasks=16 slurm-script/run_gatk_array.sh



bamfile <- list.files(path="/group/jrigrp4/teosinte-parents/20parents/bams",
                      pattern="sorted.*bam$", full.names=TRUE)
inputdf <- data.frame(
    bam=bamfile, 
    out=gsub(".*20parents/bams","/home/jolyang/Documents/Github/methylation/largedata/gatk_vcf", bamfile), 
    group="1", 
    sample=gsub(".*sorted.|_index.*", "", bamfile),
    PL="illumina", LB="lib1", PU="unit1")
inputdf$out <- gsub("sorted.|_index.*", "", inputdf$out)

rg <- list.files(path="/home/jolyang/Documents/Github/methylation/largedata/gatk_vcf",
                      pattern="RG.bam$", full.names=TRUE)
rg <- gsub(".*/|.sorted.*", "", rg)

library(farmeR)
run_GATK(inputdf[!(inputdf$sample %in% rg),], 
         ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
         gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
         picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
         minscore = 5, markDup=TRUE, addRG=TRUE, 
         realignInDels=FALSE, indels.vcf="indels.vcf",
         recalBases=FALSE, dbsnp.vcf="dbsnp.vcf", 
         email="yangjl0930@gmail.com",
         runinfo = c(TRUE, "bigmemm", 16))








