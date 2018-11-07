### Jinliang Yang
### 12-01-2016
### purpose: run BSmoothing


cl <- c("low", "high")
### additional file path
af <- c("cache/geneset_rnaseq.csv", "cache/geneset_protein.csv")
pwd <- list.files(path=c("largedata/COMET", "largedata/COMET_CHG", "largedata/COMET_CHH"), 
                  pattern="^J", full.names = TRUE)
chr <- "chr1.txt"
infile <- paste(pwd, chr, sep="/")


df <- data.frame(infile=rep(infile, each=4), af=rep(af, each=2), type=cl, 
                 output=rep(c("RNA", "protein"), each=2), 
                 context=rep(c("CHG", "CHH", "CG"), each=80))
df$output <- paste0("largedata/lcache/", df$context, "_chr1_exp_", 
                    df$output, "_", df$type, "_",
                    gsub(".*JR|\\/chr.*", "", df$infile), ".csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
df <- subset(df, af %in% "cache/geneset_protein.csv")
write.csv(df, "largedata/run_exp_df.csv")

library("farmeR")
run_Rcodes(inputdf=data.frame(file=1:120, out=1), outdir="slurm-script", cmdno=1,
           rcodes = "profiling/3.COMET/3.E.4_walley_data.R",
           arrayshid = "slurm-script/run_exp_array.sh",
           email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1, "8G"))

###>>> In this path: cd /home/jolyang/Documents/Github/methylation
###>>> RUN: sbatch -p bigmemm --time 48:00:00 --exclude=bigmem1,bigmem2 --mem 8G --ntasks=1 slurm-script/run_exp_array.sh
