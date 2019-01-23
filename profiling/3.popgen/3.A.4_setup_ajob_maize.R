### Jinliang Yang
### 01-22-2019
### set up array job

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

source("lib/doimpute.R")
#source("lib/mcmcbc.R")
library("data.table")
library("plyr")

files <- list.files(path="largedata/mr_100bp", pattern="elite_", full.names=TRUE)

mx <- fread(files[JOBID], data.table=FALSE)
out <- doimpute(mx, ncols=2:17, binsize=1000)

outfile <- gsub("matrix", "matrix_imp", files[JOBID])
fwrite(out, outfile, sep=",", row.names=FALSE, quote=FALSE)
#####

