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

#files <- list.files(path="cache/mcmc_res", pattern="new_msfs", full.names=TRUE)
df <- read.csv("slurm-script/slurm_rates_parameters.csv")
sfs <- read.csv(as.character(df$file[JOBID]))

set.seed(1234567)
# If acceptance too high, increase these values to explore wider space. If acceptance too low, decrease.
res <- MCMCBC(my_sfs=sfs$Freq, rates=c(df$rates[JOBID],df$rates[JOBID],df$rates[JOBID]), 
              sd=c(df$sd[JOBID], df$sd[JOBID], df$sd[JOBID]), k=0:(nrow(sfs)-1),
              conditional=FALSE, Ne=150000, ngen=100000, verbose=TRUE)

#####
#out <- gsub("cache", "largedata", files[JOBID])
#out <- gsub("csv", "RData", df$file[JOBID])
save(list="res", file=df$out[JOBID])
### plot trace and posteriors


