#' \code{Stacking plot}
#'
#'
#' @param df An input data.frame, with columns of chr and pos.
#' @param chrlen The chr length file, tab delimited.
#' @param centcsv Centeromere file, csv.
#' @param plotcol The column name to plot.
#' @param dorescale Whether to do rescale or not. Default=FALSE.
#' @param ... Other parameters will pass to points().
#'
#' @return return a plot.
#'
#' @examples
#'
#' @export
plot_stack <- function(df, chrlen="data/ZmB73_RefGen_v2.length",
                       centcsv=NULL, plotcol="Dm", dorescale=FALSE, ...){
    
    #source("~/Documents/Github/zmSNPtools/Rcodes/rescale.R")
    
    df <- df[order(df$chr, df$pos), ]
    #### read chr length
    cl <- read.table(chrlen, header=FALSE)
    names(cl) <- c("chrom", "BP")
    
    plot(c(0, max(cl$BP)), c(10,110), type="n",
         xlab="Physical position (kb)", ylab="", yaxt="n", bty="n")
    
    axis(side=2, tick =FALSE, las=1, at=c(105, 95, 85, 75, 65, 55, 45, 35, 25, 15),
         labels=paste("Chr", 1:10, sep=""))
    #### chromosome
    for (i in 1:10){
        lines(c(0, cl[i,]$BP), c(105-10*(i-1), 105-10*(i-1)), lwd=2, col="grey")
        #lines (c(centromere[i,]$Start,centromere[i,]$End),
        #       c(105-10*i, 105-10*i),lwd=5, col="tomato")
    }
    ### core plot
    cols <- rep(c("slateblue", "cyan4"), 5)
    for (chri in 1:10){
        mytab <- subset(df, chr == chri)
        if(dorescale){
            mytab$res <- rescale(mytab[, plotcol], c(-5, 5))
        }else{
            mytab$res <- mytab[, plotcol]
        }
        
        points(mytab$pos, 105 - 10*(chri-1) + mytab$res, pch=19, col=cols[chri], ...)
    }

    if(!is.null(centcsv)){
        #"data/AGPv2_centromere.csv"
        cen <- read.csv(centcsv)
        for(i in 1:nrow(cen)){
            rect(xleft=cen$start[i]*1000000, ybottom= 110 - 10*(cen$Chromosome[i] -1), 
                 xright=cen$end[i]*1000000, ytop= 110 - 10*(cen$Chromosome[i]), border = "red")
        }
    }
}


#http://cran.r-project.org/doc/contrib/Lemon-kickstart/rescale.R
# linearly transforms a numeric object to fit a different range.
# could use a few more sanity checks
#' @export
rescale <- function(x,newrange) {
    if(nargs() > 1 && is.numeric(x) && is.numeric(newrange)) {
        # if newrange has max first, reverse it
        if(newrange[1] > newrange[2]) {
            newmin<-newrange[2]
            newrange[2]<-newrange[1]
            newrange[1]<-newmin
        }
        xrange<-range(x)
        if(xrange[1] == xrange[2]) stop("can't rescale a constant vector!")
        mfac<-(newrange[2]-newrange[1])/(xrange[2]-xrange[1])
        invisible(newrange[1]+(x-xrange[1])*mfac)
    }
    else {
        cat("Usage: rescale(x,newrange)\n")
        cat("\twhere x is a numeric object and newrange is the min and max of the new range\n")
    }
}
