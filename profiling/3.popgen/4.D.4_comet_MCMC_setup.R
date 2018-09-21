### Jinliang Yang
### 10-11-2016
### purpose: run BSmoothing



library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:9, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/4.popgen/4.D.3_comet_MCMC.R",
           arrayshid = "slurm-script/run_rcode_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 2, "16G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --mem 16G --ntasks=2 --time 12:00:00 slurm-script/run_rcode_array.sh
