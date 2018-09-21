### Jinliang Yang
### June 16th, 2016
### purpose: get SFS for features

library(GenomicFeatures)
library("data.table")

res <- read.csv("cache/stat_exon_mean_var.csv")


RAdata <- fread("largedata/vcf_files/teo20_RA_exon.txt")


getsfs <- function(RAdata, context="CG", cols=3:22, BINSIZE=100, geneids=subset(res, mm > 0.6)$geneid){
    
    ### features
    gene <- get_feature(gff="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", features="gene")
    subgen <- subset(gene, attribute %in% geneids)
    gr0 = with(subgen, GRanges(seqnames=seqname, IRanges(start=start, end=end), strand=strand, geneid=attribute ))
    
    ### methylation
    cg <- RAdata[ V2 == context] #9270420      22
    cg$chr <- gsub("_.*", "", cg$V1)
    cg$pos <- as.numeric(as.character(gsub(".*_", "", cg$V1)))
    
    ###########
    gr1 = with(cg, GRanges(seqnames=chr, IRanges(start=pos, end=pos), snpid=V1, context=V2))
    ex1 = findOverlaps(gr0, gr1)
    #ranges(gr0)[queryHits(ex1)] = ranges(gr1)[subjectHits(ex1)]
    
    ids <- gr1$snpid[subjectHits(ex1)]
    
    subcg <- cg[V1 %in% ids]
    #gr2 <- as.data.frame(gr1)
    #gr2 <- as.data.table(gr2)
    
    subcg[subcg=="."] <- "NA"
    subcg <- as.data.frame(subcg)
    subcg[, cols] <- apply(subcg[, cols], 2, as.numeric)
    
    
    subcg$bin <- paste(subcg$chr, round(subcg$pos/BINSIZE, 0), sep="_")
    #.SD is a data.table and holds all the values of all columns, except the one specified in by.
    subcg <- as.data.table(subcg)
    res <- subcg[, lapply(.SD, mymean), by=bin, .SDcols = paste0("V", cols)]
    
    res <- as.data.frame(res)
    f <- apply(res[, -1], 1, getcount)
    sfs <- table(f)
}

mymean <- function(x){
    return(mean(x, na.rm=T))
}

getcount <- function(x, mmin=0.3, mmax=0.7){
    n0 <- sum(x < mmin)
    n2 <- sum(x > mmax)
    n1 <- sum(x >= mmin & x <= mmax)
    return(2*n2+n1)
}

###
get_feature <- function(gff="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", features="gene"){
    gff <- fread(gff)
    gff <- as.data.frame(gff)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    fe <- subset(gff, feature %in% features)
    fe$attribute <- gsub("ID=", "", fe$attribute)
    fe$attribute <- gsub(";.*", "", fe$attribute)
    fe <- subset(fe, seqname %in% 1:10)
    
    fe$strand <- as.character(fe$strand)
    return(fe)
}

sfs <- getsfs(context="CHG", cols=3:22, BINSIZE=100)
write.table(sfs, "cache/sfs_test.csv", sep=",", row.names=FALSE, quote=FALSE)

write.table(sfs, "cache/sfs_gbM_98k.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(sfs, "cache/sfs_ngbM_550k.csv", sep=",", row.names=FALSE, quote=FALSE)

########################
cg <- as.data.frame(cg)
sub <- subset(cg, chr == 10)

gbm <- read.csv("cache/sfs_gbM_98k.csv")
s1 <- ggplot(gbm, aes(x=f, y=Freq)) + geom_bar(stat="identity", position = "dodge") +
    #theme_bw() +
    xlab("") +
    ylab("Frequency") +
    scale_fill_manual(values="red") +
    
    theme(legend.position="top", axis.text=element_text(size=15), axis.title=element_text(size=15) )
#########
s1

ngbm <- read.csv("cache/sfs_ngbM_550k.csv")
s2 <- ggplot(ngbm, aes(x=f, y=Freq)) + geom_bar(stat="identity", position = "dodge") +
    #theme_bw() +
    xlab("") +
    ylab("Frequency") +
    scale_fill_manual(values="red") +
    theme(legend.position="top", axis.text=element_text(size=15), axis.title=element_text(size=15) )
#########
s2
    