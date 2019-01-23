
doimpute <- function(mx, ncols=2:21, binsize=1000){
  
  mx$chr <- gsub("_.*", "", mx$uid)
  mx$pos <- gsub(".*_", "", mx$uid)
  outdf <- mx[, c("uid", "chr", "pos")]
  
  for(i in ncols){
    sid <- names(mx)[i]
    out <- impute_onesample(mx=mx[, c("uid", sid)], sid, binsize)
    outdf <- merge(outdf, out, by="uid")
  }
  return(outdf)
}


impute_onesample <- function(mx, sid, binsize){
  # 1st impuate: missing data imputation using window=binsize:
  nanum1 <- sum(is.na(mx[,2]))
  if(nanum1 > 0){
    mx$chr <- gsub("_.*", "", mx$uid)
    mx$pos <- gsub(".*_", "", mx$uid)
    mx$bin <- paste(mx$chr, round(as.numeric(as.character(mx$pos))/binsize,0), sep="_")
    df1 <- mimpute(mymx=mx[, c("bin", names(mx)[2])])
  }
  
  # 2nd round of imputation using window=10xbinsize
  nanum2 <- sum(is.na(df1[,2]))
  if(nanum2 > 0){
    df1$chr <- gsub("_.*", "", df1$bin)
    df1$pos <- gsub(".*_", "", df1$bin)
    df1$bin2 <- paste(df1$chr, round(as.numeric(as.character(df1$pos))/10,0), sep="_")
    
    mymx2 <- df1[, c("bin2", "mean")]
    names(mymx2) <- c("bin", "mr")
    df2 <- mimpute(mymx=mymx2)
  }
  
  # 3rd round of imputation using upside down fillin  
  nanum3 <- sum(is.na(df2[,2]))
  if(nanum3 > 0){
    df2$chr <- gsub("_.*", "", df2$bin)
    df2$pos <- gsub(".*_", "", df2$bin)
    #df2$bin2 <- paste(df2$chr, round(as.numeric(as.character(df2$pos))/10,0), sep="_")
    df2 <- df2[order(df2$chr, df2$pos),]
    
    idx <- which(is.na(df2$mean))
    for(i in idx){
      df2[idx, ]$mean <- df2[idx -1, ]$mean
    }
    nanum4 <- sum(is.na(df2$mean))
    if(nanum4 > 0){
      stop("### ERROR ### fillin failed!")
    }
  }
  
  ##### ----------------------------- #########
  # backward NA replacement:
  df12 <- merge(df1, df2[, 1:2], by.x="bin2", by.y="bin", all.x=TRUE)
  df12[is.na(df12$mean.x),]$mean.x <- df12[is.na(df12$mean.x),]$mean.y
  
  mx2 <- merge(mx[, c("uid", "bin", sid)], df12[, c("bin", "mean.x")], by="bin", all.x=TRUE)
  mx2[is.na(mx2[,sid]), sid] <- mx2[is.na(mx2[,sid]),]$mean.x
  mx2 <- mx2[, c("uid", sid)]
  
  message(sprintf("###>>> finished imputing [sample:%s], NA 1:[%s], 2:[%s], 3:[%s], 4:[%s]", sid, nanum1, nanum2, nanum3, nanum4))
  return(mx2)
}

mimpute <- function(mymx=mx[, c("bin", names(mx)[2])]){
  #sum(is.na(mymx[,2]))
  #mymx <- as.data.frame(mymx)
  names(mymx)[2] <- "mr"
  df <- ddply(mymx, .(bin), summarise, 
              mean=mean(mr, na.rm=TRUE))
  return(df)
}

