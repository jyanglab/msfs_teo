mycut <- function(x){
    t <- levels(cut(range(x), breaks=10, format="d", dig.lab = 9))
    return(t)
}

find_cval_gff <- function(infile="largedata/COMET/JRA1/chr1.txt", TE=FALSE, gff,
                          features=c("exon", "intron", "up1k", "down1k", "gene", "CDS")){
    
    ## gff: require seqname, start, end, feature, strand [data.frame]
    
    df <- fread(infile)
    df$V1 <- gsub("chr", "", df$V1)
    
    #myquery$chr <- gsub("chr", "", myquery$chr)
    myout <- data.frame()
    grc <- with(df, GRanges(seqnames=V1, IRanges(start=V2, end=V3),  value=V4))
    for(fea in features){
        
        if (fea == "up1k"){
            ### upstream 1kb regions
            if(TE){
                p <- subset(gff, strand %in% "+")
            }else{
                p <- subset(gff, feature %in% "gene" & strand %in% "+")
            }
            
            p$end <- p$start - 1
            p$start <- p$start - 1001
            
            #if(TE){
            #    m <- subset(gff, strand %in% "-")
            #}else{
            #    p <- subset(gff, feature %in% "CDS" & strand %in% "+")
            #}
            
            #m$start <- m$end + 1
            #m$end <- m$end + 1001
            
            #f <- rbind(p, m)
            f <- p
            #grf <- with(onek, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
        } else if(fea == "down1k"){
            ### upstream 1kb regions
            if(TE){
                p <- subset(gff, strand %in% "+")
            }else{
                p <- subset(gff, feature %in% "gene" & strand %in% "+")
            }
           
            p$start <- p$end + 1
            p$end <- p$end + 1001
            
            #if(TE){
            #    m <- subset(gff, strand %in% "-")
            #}else{
            #    m <- subset(gff, feature %in% "gene" & strand %in% "-")
            #}
            
            #m$end <- m$start - 1001
            #m$start <- m$start - 1
            
            #f <- rbind(p, m)
            f <- p
            #grf <- with(onek, GRanges(seqnames=seqname, IRanges(start=start, end=end) ))
        }else{
            f <- subset(gff, feature %in% fea)
        }
        
        ###### prepare normalized bins
        out <- t(apply(as.matrix(f[, c("start", "end")]), 1, mycut))
        out <- as.data.frame(out)
        out$seqname <- f$seqname
        
        lout <- gather(out, key="var", value="range", V1:V10)
        lout$start <-  as.numeric(gsub("\\(|\\,.*", "", lout$range))
        lout$end <- as.numeric(gsub(".*\\,|\\]", "", lout$range))
        lout <- lout[order(lout$start), ]
        
        lout$start <- round(lout$start, 0)
        lout$start <- lout$start + 1
        lout$end <- round(lout$end, 0)
        
        
        grf <- with(lout, GRanges(seqnames=seqname, IRanges(start=start, end=end), v=var))
        
        message(sprintf("###>>> compute for [sample=%s] and [feature=%s]", infile, fea))
        for(v in paste0("V", 1:10)){    
            mygrf <- grf[ grf$v %in% v ]
            #######
            tb <- findOverlaps(query=mygrf, subject=grc)
            tb <- as.matrix(tb)
            
            out <- as.data.frame(grc[unique(tb[,2])])
            
            tem <- data.frame(file=infile, feature=fea, mc=mean(out$value), vc=var(out$value), var=v)
            myout <- rbind(myout, tem)
        }
        
    }
    return(myout)
}

