# Jinliang Yang
# Purpose: quick plot of GWAS results
# start: 2.11.2012
# updated: 7/14/2014
# add the RNA-seq LM regression data

# location: 129.186.85.7
methy_plot <- function(res=m1, cex=.9, pch=16, col=rep(c("slateblue", "cyan4"), 5), 
                         GAP=5e+06, yaxis=NULL,
                         col2plot="ModelFreq", ... ){
  
  source("~/Documents/Github/zmSNPtools/Rcodes/newpos.R")
  res <- newpos(res, GAP = GAP, version = "v3")
  chrtick <- chrline_tick(GAP = GAP, version = "v3")
  
  sub <- seq(1, 16389921, by=1000000)
  res1 <- subset(res, chr == 10)
  BINSIZE = 1000
  res1$bin <- round(res1$pos/BINSIZE, 0)
  library(plyr)
  res1$B73_ratio <- as
  mymean <- function(x){
      return(mean(x, na.rm=TRUE))
  }
  tab <- ddply(res1, .(bin), summarise,
               ratio = mean(B73_ratio, na.rm=TRUE)
              )
  
  #### setup the cavon
  if(is.null(yaxis)){
    plot(x=res1$pos, y=res1$B73_ratio,  type="h", xaxt="n", xlab="", 
         xlim=c(0, max(chrtick$chrlines)), ylim=c(0, 1.3) )
    
    makeTransparent = function(..., alpha=0.5) {
        
        if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
        
        alpha = floor(255*alpha)  
        newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
        
        .makeTransparent = function(col, alpha) {
            rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
        }
        
        newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
        
        return(newColor)
        
    }  
    makeTransparent("blue", alpha=0.1)
    plot(x=tab$bin*1000, y= (1-tab$ratio),  type="h", xaxt="n", xlab="",  ylim=c(0, 1.3), cex=0.1, col="#0000FF19")
    ,
         ...)
  }else{
    plot(x=res, y=-1000,  type="h", xaxt="n", yaxt="n", xlab="",
         xlim=c(0, max(chrtick$chrlines)),
         ...)
    axis(side=2, at=yaxis, labels=yaxis)
  }
  axis(side=1, at=chrtick$ticks, labels=c("chr1", "chr2", "chr3", "chr4", "chr5", 
                                          "chr6", "chr7", "chr8", "chr9", "chr10"))
  abline(v=chrtick$chrlines, col="grey")
  
  for(i in 1:10){
    points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
         pch = pch, col=col[i], cex=cex);
  }  
}
