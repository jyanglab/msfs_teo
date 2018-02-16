### Jinliang Yang
### June 7th, 2016


source("lib/run_pseudoref.R")
library(farmeR)
inputdf <- data.frame(input.vcf="largedata/gatk_vcf/JRIAL1A_joint_call.filtered_snps.vcf",
                      out.fa="largedata/pgenome/mysample.fa")
run_pseudoref(inputdf,
              ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
              gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
              email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemh", 8))



### run bsmap
library(farmeR)
inputdf <- data.frame(fq1="/group/jrigrp4/BS_teo20/fastq/JRA1_CTTGTA_R1.fastq.gz",
                      fq2="/group/jrigrp4/BS_teo20/fastq/JRA1_CTTGTA_R2.fastq.gz",
                      out="$HOME/Documents/Github/methylation/largedata/bsmap/JRA1")

runa_bsmap(inputdf, ref.fa="$HOME/Documents/Github/methylation/largedata/pgenome/mysample.fa",
           picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
           email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 8))
