library("GenomicRanges")
#library("data.table")
library(plyr)


find_overlap_win_gff <- function(df, gff){
    
    ## df: chr, start, end, ratio (mehtylation ratio) [data.frame]
    ## gff: require seqname, start, end, exid (transcript id), exonid=exonid [data.frame]
    
    df$chr <- gsub("Chr", "", df$chr)
    
    grc <- with(df, GRanges(seqnames=chr, IRanges(start=start, end=end),  value=ratio))
    grf <- with(gff, GRanges(seqnames=seqname, IRanges(start=start, end=end), txid=txid, exonid=exonid))
    
    #######
    tb <- findOverlaps(query=grf, subject=grc)
    tb <- as.matrix(tb)
    
    out1 <- as.data.frame(grf[tb[,1]])
    out2 <- as.data.frame(grc[tb[,2]])
    
    out <- cbind(out1, out2[, "value"])
    names(out)[8] <- "value"
    
    out$exonid <- paste(out$txid, out$exonid, sep="_")
    
    mtx <- ddply(out, .(txid), summarise,
                 mean = mean(value, na.rm = TRUE),
                 var = var(value, na.rm=TRUE))
    mexon <- ddply(out, .(exonid), summarise,
                 mean = mean(value, na.rm = TRUE),
                 var = var(value, na.rm=TRUE))
    
    myout <- list(mtx, mexon)
    return(myout)
}

