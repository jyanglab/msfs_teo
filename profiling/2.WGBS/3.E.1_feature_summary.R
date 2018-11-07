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

runit_gene <- function(df, id){
    mypwd <- as.character(df$pwd[id])
    #mytype <- as.character(df$type[id])
    myout <- as.character(df$output[id]) 
    
    ######## format GFF and repeat files
    gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    info <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt", data.table=FALSE)
    info <- subset(info, is_canonical %in% "yes")
    
    ### we only consider positive strand!!!!
    gff <- subset(gff, strand %in% "+")
    
    
    exs <- subset(gff, feature %in% "exon")
    exs$id <- gsub("Parent=|;Name=.*", "", exs$attribute)
    exs <- subset(exs, id %in% info$transcript_id)
    # get the first exon
    exs1 <- exs[order(exs$seqname, exs$start), ]
    exs1 <- exs1[!duplicated(exs1$id), ]
    exs1$feature <- "exon1st"
    
    exs2 <- exs[order(exs$seqname, exs$start, decreasing = TRUE), ]
    exs2 <- exs2[!duplicated(exs2$id), ]
    exs2$feature <- "exonlast"
    message(sprintf("###>>> get [total=%s] canonical exons: [1st=%s] and [last=%s]",
                    length(unique(exs$id)), nrow(exs1), nrow(exs2) ))
    
    ins <- subset(gff, feature %in% "intron")
    ins$id <- gsub("Parent=|;Name=.*", "", ins$attribute)
    ins <- subset(ins, id %in% info$transcript_id)
    ### get first intron
    ins1 <- ins[order(ins$seqname, ins$start), ]
    ins1 <- ins1[!duplicated(ins1$id), ]
    ins1$feature <- "intron1st"
    
    ins2 <- ins[order(ins$seqname, ins$start, decreasing = TRUE), ]
    ins2 <- ins2[!duplicated(ins2$id), ]
    ins2$feature <- "intronlast"
    message(sprintf("###>>> get [total=%s] canonical Introns: [1st=%s] and [last=%s]",
                    length(unique(ins$id)), nrow(ins1), nrow(ins2) ))
    
    
    gen <- subset(gff, feature %in% "gene")
    gen$id <- gsub("ID=|;Name=.*", "", gen$attribute)
    message(sprintf("###>>> get [total=%s] canonical gene on positive strand!",
                    nrow(gen) ))
    
    newgff <- rbind(gen, exs1, exs2, ins1, ins2)
    
    ####### smoothed files CG
    #pwd1 <- list.files(path=mypwd, pattern="^J", full.names = TRUE)
    
    #res1 <- data.frame()
    #for(i in 1:length(pwd1)){
    #    out1 <- find_cval_gff(infile=paste0(pwd1[i], "/chr1.txt"), gff=newgff, TE=FALSE,
    #                          features=c("exon1st", "exonlast", "gene", "intron1st", "intronlast", "up1k", "down1k"))
    #    res1 <- rbind(res1, out1)
    #}
    out <- find_cval_gff(infile=as.character(df$infile[id]), gff=newgff, TE=FALSE,
                         features=c("exon1st", "exonlast", "gene", "intron1st", "intronlast", "up1k", "down1k"))
    write.csv(out, myout)
}


## control elements
df <- read.csv("largedata/run_gene_df.csv")
# col, pwd="largedata/COMET"
# col: output="cache/CG_chr1_TE_class1.csv"
runit_gene(df, id=JOBID)










