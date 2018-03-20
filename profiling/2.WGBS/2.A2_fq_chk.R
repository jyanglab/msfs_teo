### Jinliang Yang
### check the quality of WGBS data
### 3/22/2016

library("farmeR")

### rename gzipped files
files <- list.files(path="largedata/fastq/", pattern="fastq.lz4.fastq.gz", full.names=T)
file.rename(from=files, to=gsub("fastq.lz4.", "", files))

files <- list.files(path="largedata/fastq/", pattern="gz.qc", full.names=T)
file.rename(from=files, to=gsub(".gz.qc", "fastq.gz.qc", files))

### get information
fq <- list.files(path="largedata/fastq", pattern="fastq.gz$", full.names = TRUE)
df <- data.frame(fq=fq, out=paste0("largedata/fastq/", gsub(".*/", "", fq), ".qc"))

### run quality checking
run_fastq_qc(df, q=20, email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 2))
    
    
### get summary stat
files <- list.files(path="largedata/fastq", pattern="qc$", full.names=TRUE)
res <- get_qc(files, genomesize=2500000000)
write.table(res, "cache/wgbs_teo20_fastq_qc.csv", sep=",", row.names=FALSE, quote=FALSE)

### plot results

res <- read.csv("cache/wgbs_teo20_fastq_qc.csv")

idx <- seq(from=1, to =40, by=2)
dp <- res[idx,]
dp$depth <- 2*dp$depth
dp$file <- gsub(".*fastq/|_.*", "", dp$file)

par(mfrow=c(1,2))
barplot(dp$avgQ, names.arg=dp$file, col="#bdb76b", las=2, main="Avg Base Quality (N=20)")
barplot(dp$depth, names.arg=dp$file, col="#698b69", las=2, main="Read Depth (N=20)")


