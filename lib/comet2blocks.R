### Jinliang Yang
### 10-12-2016
### Chop COMET into shared blocks
### updated 12-06-2016: no gap and dynamic cutoffs


comet2blocks <- function(files, chri=10, collapse=TRUE, verbose=TRUE, cutoff=c(0.02, 0.60)){
    
    chrlen <- read.table("~/dbcenter/AGP/AGPv2/ZmB73_RefGen_v2.fasta.fai", header=FALSE)
    names(chrlen)[1:2] <- c("chr", "len")
    cometls <- list()
    for(i in 1:length(files)){
        df <- fread(files[i], data.table=FALSE)
        # https://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
        #library(mixtools)
        if(is.null(cutoff)){
            res <- normalmixEM(x=df$meth, mu=c(0.01, 0.5, 0.8), k =3)
            #lambda The final mixing proportions.
            #mu The final mean parameters.
            #sigma The final standard deviations. If arbmean = FALSE, then only the smallest standard
            #deviation is returned. See scale below
            
            #lambda <- summary(res)
            mycutoff <- quantile(df$meth, cumsum(res$lambda))
        }else{
            mycutoff <- cutoff
        }
        sid <- gsub(".*/|_.*", "", files[i])
        
        if(verbose) message(sprintf("[comet2blocks] sample [%s], recalculating low, med, high, using cutoff [ <= %s and >= %s] ... ",
                        sid, round(mycutoff[1], 3), round(mycutoff[2], 3)))
        
        df <- subset(df, chr == chri)
        df$level <- "med"
        df[df$meth <= mycutoff[1], ]$level <- "low"
        df[df$meth >= mycutoff[2], ]$level <- "high"
        
       
        if(verbose) message(sprintf("[comet2blocks] close gaps ..."))
        ### first block 
        
        df$start[1] <- 1
        df$stop[nrow(df)] <- subset(chrlen, chr == chri)$len
        df$start[2:nrow(df)] <- df$stop[1:(nrow(df)-1)] + 1
        df$newsize <- df$stop - df$start + 1
        
        if(collapse == TRUE){
            if(verbose) message(sprintf("[comet2blocks] condense for [chr %s] ... ", chri))
            cometls[[sid]] <- condense(df)
        }else{
            cometls[[sid]] <- df
        }
        
    }
    
    #message(sprintf("[comet2blocks] chop into shared blocks ... ")
    res <- chop2blocks(cometls, chri)
    res <- apply(res, 2, as.character)
    res[res=="high"] <- 2
    res[res=="med"] <- 1
    res[res=="low"] <- 0
    res[res=="a"] <- -9
    return(res)
}



## condense COMET by chr
condense <- function(df){
    
    df <- df[order(df$start), ]

    df$level <- as.character(df$level)
    df$level0 <- c("no", df$level[-nrow(df)])
    df$eval <- 1
    df[df$level == df$level0, ]$eval <- 0
    df$block <- cumsum(df$eval)
    
    out <- ddply(df, .(block, level), summarise,
                 start=min(start),
                 stop=max(stop),
                 meth=mean(meth))
    return(out)
}

## chop into shared blocks
chop2blocks <- function(cometls, chri){
    
    ##### determine the break points
    out <- lapply(1:length(cometls), function(x){
        return(unique(c(cometls[[x]]$start, cometls[[x]]$stop+1)))
    })
    bp <- c()
    for(i in 1:length(out)){
        bp <- c(bp, out[[i]])
    }
    bp <- sort(unique(bp))
    message(sprintf("[chop2blocks] chop [chr%i] into [ %s ] shared blocks ... ", chri, length(bp)))
    
    ### use genomicranges to find overlaps
    gr <- GRanges(seqnames=Rle(paste0("chr", chri)),
                  ranges=IRanges(start=bp[-length(bp)], end=bp[-1]-1),
                  strand = Rle(strand("+")) )
    
    out <- data.frame(chr=seqnames(gr), bid=paste(start(gr), end(gr), sep="_"))
    for(j in 1:length(cometls)){
        message(sprintf("[chop2blocks] working on [ sample: %s] ... ", names(cometls[j])))
        ### use genomicranges to find overlaps
        gr1 <- GRanges(seqnames=paste0("chr", chri),
                       ranges=IRanges(start=cometls[[j]]$start, end=cometls[[j]]$stop),
                       strand = Rle(strand("+")),
                       score = cometls[[j]]$meth,
                       level = cometls[[j]]$level)
        mygr <- gr
        idx <- findOverlaps(query=gr1, subject=mygr, select="all")
        #mygr$score <- -9
        #mygr[subjectHits(idx), ]$score <- gr1[queryHits(idx), ]$score
        mygr$level <- "a"
        mygr[subjectHits(idx), ]$level <- gr1[queryHits(idx), ]$level
        
        mx <- data.frame( bid=paste(start(mygr), end(mygr), sep="_"), level=mcols(mygr)$level)
        names(mx)[2] <- names(cometls[j])
        out <- merge(out, mx, by="bid", sort=FALSE)
    }
    return(out)
}



