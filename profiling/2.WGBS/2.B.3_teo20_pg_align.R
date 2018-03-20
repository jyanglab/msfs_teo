### Jinliang Yang
### July 25th, 2016

library(farmeR)

idtab <- read.csv("data/teo20_ids.csv")
idtab$idchar <- gsub("-", "", idtab$idchar)
idtab$pid <- paste0("JR", idtab$plate)

########### Note, move fq files to /group/jrigrp7
########### alignment
fq1 <- list.files(path="/group/jrigrp3/jyang/methylation/largedata/wgbs_fq", pattern="R1_trimmed.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="/group/jrigrp3/jyang/methylation/largedata/wgbs_fq", pattern="R2_trimmed.fastq.gz$", full.names = TRUE)

#bamfiles <- list.files(path="/group/jrigrp4/BS_teo20/WGBS/BSM", pattern="bam$", full.names = TRUE)
### note: for alignment, "bam" col should not be present.
inputdf <- data.frame(fq1 = fq1,  fq2 = fq2, outbase = gsub(".*/|_.*", "", fq1), out2=gsub(".*/|_.*", "", fq2))

if(sum(as.character(inputdf$outbase) != as.character(inputdf$out2)) > 0  ){stop("!!! PE not right !!!")}

inputdf <- merge(inputdf, idtab[, c("idchar", "pid")], by.x="outbase", by.y="pid")
inputdf$genome <- paste0("$HOME/Documents/Github/methylation/largedata/wgbs_pgen/", inputdf$idchar)

### AGPv2
run_bismark(inputdf[-1,], genome = NULL,
            outdir = "/home/jolyang/Documents/Github/methylation/largedata/wgbs_align", 
            N = 1, align = TRUE,
            email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 16))

