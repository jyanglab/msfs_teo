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

files <- list.files(path="largedata/mr_100bp", pattern="elite_rm_", full.names=TRUE)

mx <- fread(files[JOBID], data.table=FALSE)
out <- doimpute(mx, ncols=2:15, binsize=1000)

outfile <- gsub("matrix", "matrix_imp", files[JOBID])
fwrite(out, outfile, sep=",", row.names=FALSE, quote=FALSE)
#####


### manullay did the following for CpG!
#out1 <- doimpute(mx, ncols=2:7, binsize=1000)

# outfile <- gsub("matrix", "matrix_imp2-7", files[JOBID])
# fwrite(out1, outfile, sep=",", row.names=FALSE, quote=FALSE)
#####

# out1 <- doimpute(mx, ncols=8:15, binsize=1000)

# outfile <- gsub("matrix", "matrix_imp8-15", files[JOBID])
# fwrite(out1, outfile, sep=",", row.names=FALSE, quote=FALSE)
#####

#id,CpG_A632,CpG_B37,CpG_B73,CpG_B97,CpG_CML322,CpG_HP301,CpG_IL14,CpG_LH123HT,CpG_LH82,CpG_Mo17,CpG_Oh43,CpG_P39,CpG_PHZ51,CpG_Tx303
#uid,chr,pos,CpG_A632,CpG_B37,CpG_B73,CpG_B97,CpG_CML322,CpG_HP301
#CpG_IL14,CpG_LH123HT,CpG_LH82,CpG_Mo17,CpG_Oh43,CpG_P39,CpG_PHZ51,CpG_Tx303

#imp1 <- fread("largedata/mr_100bp/elite_rm_CpG_matrix_imp2-7.csv", data.table=FALSE)
#imp2 <- fread("largedata/mr_100bp/elite_rm_CpG_matrix_imp8-15.csv", data.table=FALSE)


#out <- merge(imp1, imp2[, -2:-3], by="uid")
#fwrite(out, "largedata/mr_100bp/elite_rm_CpG_matrix_imp.csv", sep=",", row.names=FALSE, quote=FALSE)

#cg <- fread("largedata/mr_100bp/elite_rm_CpG_matrix.csv", data.table=FALSE)
#cg <- cg[!duplicated(cg$uid), ]
#fwrite(cg, "largedata/mr_100bp/elite_rm_CpG_matrix_dedup.csv", sep=",", row.names=FALSE, quote=FALSE)


#cg <- fread("largedata/mr_100bp/elite_CHG_matrix_imp.csv", data.table=FALSE)

#idx <- which(duplicated(imp1$uid))
#imp1[idx[1:10], ]

