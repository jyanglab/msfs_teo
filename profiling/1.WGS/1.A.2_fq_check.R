### Jinliang Yang
### use Li Heng's De novo assembly based variant calling pipeline
### 3/22/2016

fq <- list.files(path="/group/jrigrp/teosinte-parents/seq-merged", pattern="fastq.gz$", full.names = TRUE)
df <- data.frame(fq=fq, out=paste0("largedata/fqchk/", gsub(".*/", "", fq), ".qc"))

### run quality checking
run_fastq_qc(df, q=20, email=NULL, runinfo = c(FALSE, "bigmemh", 1))
    
    
### get summary stat
files <- list.files(path="largedata/fqchk", pattern="qc", full.names=TRUE)
res <- get_qc(files, genomesize=2500000000)
write.table(res, "cache/teo20_fastq_qc.csv", sep=",", row.names=FALSE, quote=FALSE)

### plot results

res <- read.csv("cache/teo20_fastq_qc.csv")

idx <- seq(from=1, to =40, by=2)
dp <- res[idx,]
dp$depth <- 2*dp$depth
dp$file <- gsub(".*Sample_|_index.*", "", dp$file)

par(mfrow=c(1,2))
barplot(dp$avgQ, names.arg=dp$file, col="#8b2323", las=2, main="Avg Base Quality (N=20)")
barplot(dp$depth, names.arg=dp$file, col="#458b74", las=2, main="Read Depth (N=20)")


