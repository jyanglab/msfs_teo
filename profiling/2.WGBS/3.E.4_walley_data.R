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
    
    mytype <- as.character(df$type[id])
    myout <- as.character(df$output[id]) 
    myaf <- as.character(df$af[id]) 
    ######## format GFF and repeat files
    # AGPv2 three annotation files
    # additional data file
    d <- read.csv(myaf)
    #protein <- read.csv("cache/geneset_protein.csv")
    
    ######## format GFF and repeat files
    gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    info <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt", data.table=FALSE)
    info <- subset(info, is_canonical %in% "yes")
    
    ### we only consider positive strand!!!!
    gff <- subset(gff, strand %in% "+")
    
    gene <- subset(gff, feature %in% "gene")
    gene$id <- gsub("ID=|;Name=.*", "", gene$attribute)
    message(sprintf("###>>> get [total=%s] canonical gene on positive strand!",
                    nrow(gene) ))
    gene <- merge(gene, d, by.x="id", by.y="tracking_id")
    
    med <- median(log2(gene$m))
    gene[gene$m < med, ]$feature <- "low"
    gene[gene$m >= med, ]$feature <- "high"
    
    ####### smoothed files CG
    out <- find_cval_gff(infile=as.character(df$infile[id]), gff=subset(gene, feature %in% mytype), TE=TRUE,
                         features=c(mytype, "up1k", "down1k"))
    write.csv(out, myout)
}


## control elements
df <- read.csv("largedata/run_exp_df.csv")
# col, pwd="largedata/COMET"
# col: type="class=I"
# col: output="cache/CG_chr1_TE_class1.csv"
runit(df, id=JOBID)










