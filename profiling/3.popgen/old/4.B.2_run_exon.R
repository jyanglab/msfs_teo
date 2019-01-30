### Jinliang Yang
### June 20th, 2016

source("lib/mplots.R")
source("lib/mcmcbc.R")


sfs <- read.csv("cache/sfs_exon_cg.csv")

pdf(file="largedata/sfs_exon.pdf", width=8, height=8)
plot(sfs$Freq~(c(0:40)), pch=19, cex=2, ylab="counts", xlab="number of chromosomes", cex.lab=1.5)
dev.off()

# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res <- MCMCBC(my_sfs=sfs$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)


#####
save(list="res", file="largedata/res_k40_exon5.RData")
### plot trace and posteriors

ob <- load("largedata/res_k40_exon5.RData")
mplot(res, burnin=0.2, rates=c(1E8,1E8,1E8))
plot(ltrace)

### plot obs and post SFS
sfsplot(res, k=0:40)

##################################################
## Gene Body Methylation
sfs <- read.csv("cache/sfs_gbM_98k.csv")
plot(sfs$Freq~(c(0:40)), pch=19, cex=2, ylab="counts", xlab="number of chromosomes", cex.lab=1.5)

#as before, decreasing these values will increase acceptance rates and vice versa.
res <- MCMCBC(my_sfs=sfs$Freq, rates=c(1E6,1E7,1E5), sd=c(0.1,0.1,0.1), k=0:40,
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)
#posterior mu [ 8.4268710137098e-07 ], nu [ 1.42911732377031e-07 ] and s [ 1.01127321650216e-05 ]
mplot(res, burnin=0.25, rates=c(1E7,1E8,1E6))
sfsplot(res, k=0:40)

### plot obs and post SFS
## non-Gene Body Methylation
sfs <- read.csv("cache/sfs_ngbM_550k.csv")
plot(sfs$Freq~(c(0:40)), pch=19, cex=2, ylab="counts", xlab="number of chromosomes", cex.lab=1.5)

#as before, decreasing these values will increase acceptance rates and vice versa.
res <- MCMCBC(my_sfs=sfs$Freq, rates=c(1E6,1E7,1E5), sd=c(0.05,0.05,0.05), k=0:40,
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)

mplot(res, burnin=0.25, rates=c(1E6,1E7,1E5))
#posterior mu [ 8.45688646566231e-07 ], nu [ 9.33374738987097e-08 ] and s [ 2.56883378284637e-06 ]
sfsplot(res, k=0:40)
