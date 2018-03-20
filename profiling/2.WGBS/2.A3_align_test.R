### Jinliang Yang
### Updates: June 24th, 2016

library("farmeR")

########### alignment
fq1 <- list.files(path="largedata/fastq", pattern="R1.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="largedata/fastq", pattern="R2.fastq.gz$", full.names = TRUE)

#bamfiles <- list.files(path="/group/jrigrp4/BS_teo20/WGBS/BSM", pattern="bam$", full.names = TRUE)

### note: for alignment, "bam" col should not present.
inputdf <- data.frame(fq1 = fq1[1],  fq2 = fq2[1], outbase = "JRA1")

### AGPv2
run_bismark(inputdf, genome = "/home/jolyang/dbcenter/AGP/AGPv2",
            outdir = "/home/jolyang/Documents/Github/methylation/largedata/bismark", 
            N = 1, align = TRUE,
            email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 16))

### pesudo-ref
### note: for alignment, "bam" col should not present.
inputdf <- data.frame(fq1 = fq1[1], fq2 = fq2[1], outbase = "pg_JRA1")
run_bismark(inputdf, genome = "/home/jolyang/Documents/Github/methylation/largedata/pgenome",
            outdir = "/home/jolyang/Documents/Github/methylation/largedata/bismark", 
            N = 1, align = TRUE,
            email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 32))
