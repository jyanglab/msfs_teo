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

gvcf <- list.files(path="largedata/gatk_vcf",
                   pattern="g.vcf$", full.names=TRUE)
sample1 <- gsub(".*/|.g.vcf", "", gvcf)

inputdf <- subset(inputdf, !(sample %in% sample1))

###########
library(farmeR)
run_GATK(inputdf, 
         ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
         gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
         picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
         minscore = 5, markDup=TRUE, addRG=TRUE, 
         realignInDels=FALSE, indels.vcf="indels.vcf",
         recalBases=FALSE, dbsnp.vcf="dbsnp.vcf", 
         email="yangjl0930@gmail.com",
         runinfo = c(TRUE, "bigmemm", 32))


############







