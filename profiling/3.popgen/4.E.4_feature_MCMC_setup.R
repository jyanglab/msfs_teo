### Jinliang Yang
### 10-11-2016
### purpose: run BSmoothing



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:3, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/4.popgen/4.E.3_feature_MCMC.R",
           arrayshid = "slurm-script/run_rcode_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemh", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 16G --ntasks=2 --time 12:00:00 slurm-script/run_rcode_array.sh

run_Rcodes(inputdf=data.frame(file=1:3, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/4.popgen/4.E.3_feature_MCMC_gbM.R",
           arrayshid = "slurm-script/run_rcode_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --time=4:00:00 slurm-script/run_rcode_array.sh