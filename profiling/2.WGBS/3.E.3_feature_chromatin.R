### Jinliang Yang
### 02-14-2017
### purpose: chromatin status

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
    
    bed <- fread("largedata/AP.bfthresh1.1.MNaseHS.Ranges.dat", header=FALSE, data.table=FALSE)
    #names(repeats) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    names(bed) <- c("seqname", "start", "end")
    bed$feature <- "open"
    
    newbed <- data.frame()
    for(i in 1:10){
        sub <- subset(bed, seqname == i)
        sub <- sub[order(sub$start), ]
        sub$ns <- c(1, sub$end + 1)[-(nrow(sub)+1)]
        sub$ne <- sub$start
        newbed <- rbind(newbed, sub[, c("seqname", "ns", "ne", "feature")])
    }
    newbed$feature <- "closed"
    names(newbed)[2:3] <- c("start", "end")
    
    gff <- rbind(bed, newbed)
    gff$strand <- "+"
    
    
    mytype <- as.character(df$type[id])
    myout <- as.character(df$output[id]) 
    ####### smoothed files CG
    out <- find_cval_gff(infile=as.character(df$infile[id]), gff=subset(gff, feature %in% mytype), TE=TRUE,
                         features=c(mytype, "up1k", "down1k"))
    write.csv(out, myout)
}


## control elements
df <- read.csv("largedata/run_chromatin_df.csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
runit(df, id=JOBID)










