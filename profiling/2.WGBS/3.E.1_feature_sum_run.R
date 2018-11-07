### Jinliang Yang
### 12-01-2016
### purpose: run BSmoothing

pwd <- list.files(path=c("largedata/COMET", "largedata/COMET_CHG", "largedata/COMET_CHH"), 
                  pattern="^J", full.names = TRUE)
chr <- "chr1.txt"
infile <- paste(pwd, chr, sep="/")


df <- data.frame(infile=rep(infile, each=1), output=c(1,2,3), context=rep(c("CHG", "CHH", "CG"), each=20))
df$output <- paste0("largedata/lcache/", df$context, "_chr1_fea_", gsub(".*JR|\\/chr.*", "", df$infile), ".csv")
# col, infile="largedata/COMET/chr1.txt"
# col: output="cache/CG_chr1_TE_class1.csv"
write.csv(df, "largedata/run_gene_df.csv")

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:60, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.COMET/3.E.1_feature_summary.R",
           arrayshid = "slurm-script/run_fea_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --time 24:00:00 --mem 8G --ntasks=1 --exclude=bigmem1,bigmem2 slurm-script/run_fea_array.sh
