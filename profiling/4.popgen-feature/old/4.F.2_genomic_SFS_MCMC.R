### Jinliang Yang
### 10-14-2016
### run MCMC for comet on chr10

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)

source("lib/mplots.R")
source("lib/mcmcbc.R")

sfs <- read.csv("cache/gbm_sfs_genomic.csv")

# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res <- MCMCBC(my_sfs=sfs$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=1:40,
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)

#####
out <- gsub("cache", "largedata", files[JOBID])
out <- gsub("csv", "RData", out)
save(list="res", file=out)
### plot trace and posteriors

sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=1:40)
mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8))
# posterior mu [ 1.09605214547552e-06 ], nu [ 4.26242752079954e-08 ] and s [ 1.63352605892327e-05 ]
# posterior mu [ 3.00331157824573e-09 ], nu [ 1.26187710000336e-08 ] and s [ 1.12704938374162e-09 ]


sfs <- read.csv("cache/non-gbm_sfs_genomic.csv")
# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res <- MCMCBC(my_sfs=sfs$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=1:40,
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)

save(list="res", file="cache/res_non-gbs_sfs_genomic.RData")

ob <- load("cache/res_non-gbs_sfs_genomic.RData")

sfsplot(res,burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot="plotmean", Ne=150000, k=0:39)
mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8))
# posterior mu [ 1.09605214547552e-06 ], nu [ 4.26242752079954e-08 ] and s [ 1.63352605892327e-05 ]
# posterior mu [ 3.00331157824573e-09 ], nu [ 1.26187710000336e-08 ] and s [ 1.12704938374162e-09 ]


