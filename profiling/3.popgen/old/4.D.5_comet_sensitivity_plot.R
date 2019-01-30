### Jinliang Yang
### 10-14-2016
### run MCMC for comet on chr10

source("lib/mplots.R")
source("lib/mcmcbc.R")


ob <- load("largedata/sfs_chr10_comet_blocks_0.15.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
# posterior mu [ 1.09605214547552e-06 ], nu [ 4.26242752079954e-08 ] and s [ 1.63352605892327e-05 ]

ob <- load("largedata/sfs_chr10_comet_blocks_0.2.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
# posterior mu [ 4.64966149155177e-06 ], nu [ 2.30681386298064e-05 ] and s [ 2.42128603993789e-10 ]

ob <- load("largedata/sfs_chr10_comet_blocks_0.25.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))

ob <- load("largedata/sfs_chr10_comet_blocks_0.3.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
#posterior mu [ 1.53776232028128e-06 ], nu [ 1.24345606499969e-08 ] and s [ 2.09458340198109e-05 ]

ob <- load("largedata/sfs_chr10_comet_blocks_0.33.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
#posterior mu [ 1.18009029537551e-06 ], nu [ 3.33944547293781e-08 ] and s [ 1.72989849299552e-05 ]

ob <- load("largedata/sfs_chr10_comet_blocks_0.35.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
#posterior mu [ 1.02609402217551e-06 ], nu [ 5.61823394416025e-08 ] and s [ 1.52640860569889e-05 ]

ob <- load("largedata/sfs_chr10_comet_blocks_0.4.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
#posterior mu [ 8.27072344436003e-07 ], nu [ 1.37035497080222e-07 ] and s [ 1.17645674736642e-05 ]

ob <- load("largedata/sfs_chr10_comet_blocks_0.45.RData")
mplot(res, burnin=0.25, rates=c(1E8,1E8,1E8))
#posterior mu [ 7.45926304020419e-07 ], nu [ 2.30426848262742e-07 ] and s [ 9.50313895797221e-06 ]





