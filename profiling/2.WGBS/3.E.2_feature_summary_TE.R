### Jinliang Yang
### 02-03-2017
### purpose: using smoothed data to summarize methylation levels

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

JOBID <- as.numeric(as.character(args[1]))
print(JOBID)


library("data.table")
library("GenomicRanges")
library("tidyr")
source("lib/find_cval_gff.R")


#########################################
runit <- function(df, id){
    ######## format GFF and repeat files
    # AGPv2 three annotation files
    
    repeats <- fread("~/dbcenter/AGP/AGPv2/repeats/ZmB73_5a_MTEC_repeats.gff", header=FALSE, data.table=FALSE)
    names(repeats) <- c("seqname", "source", "feature", "start", "end", "score",
                        "strand", "frame", "attribute")
    repeats$class <- gsub(";.*", "", repeats$attribute)
    gff <- repeats[, -3]
    names(gff)[ncol(gff)] <- "feature"
    gff <- subset(gff, strand %in% "+")
    
    mytype <- as.character(df$type[id])
    myout <- as.character(df$output[id]) 
    ####### smoothed files CG
    out <- find_cval_gff(infile=as.character(df$infile[id]), gff, TE=TRUE,
                         features=c(mytype, "up1k", "down1k"))
    write.csv(out, myout)
}


## control elements
df <- read.csv("largedata/run_TE_df.csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
runit(df, id=JOBID)










