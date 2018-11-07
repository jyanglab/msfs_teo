### Jinliang Yang
### 12-01-2016
### purpose: run BSmoothing


cl <- c("open", "closed")
pwd <- list.files(path=c("largedata/COMET", "largedata/COMET_CHG", "largedata/COMET_CHH"), 
                  pattern="^J", full.names = TRUE)
chr <- "chr1.txt"
infile <- paste(pwd, chr, sep="/")


df <- data.frame(infile=rep(infile, each=2), type=cl, output=c(1,2,3), 
                 context=rep(c("CHG", "CHH", "CG"), each=40))
df$output <- paste0("largedata/lcache/", df$context, "_chr1_chromtin_", df$type, "_",
                    gsub(".*JR|\\/chr.*", "", df$infile), ".csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
write.csv(df, "largedata/run_chromatin_df.csv")

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:120, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.COMET/3.E.3_feature_chromatin.R",
           arrayshid = "slurm-script/run_chromatin_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --time 48:00:00 --exclude=bigmem1,bigmem2 --mem 8G --ntasks=1 slurm-script/run_chromatin_array.sh
