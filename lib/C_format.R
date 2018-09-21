library("data.table")
library(plyr)
C_format <- function(inputdf, myr,  minsite=10, verbose=FALSE){
    #files <- list.files(path=inputdir, full.names=TRUE)
    #df1 <- df2 <- df3 <- data.frame()
    inputdf$file <- as.character(inputdf$file)
    for(i in myr){
    #res <- lapply(myr, function(i){
        onegene <- try(fread(inputdf$file[i], header=FALSE), silent = TRUE)
        if (!inherits(onegene, 'try-error')){
            onegene <- as.data.frame(onegene)
            onegene$miss <- apply(onegene, 1, function(x) sum(x=="."))
            onegene <- subset(onegene, miss==0)
            
            #chg <- subset(onegene, V3 == "CHG")
            #chh <- subset(onegene, V3 == "CHH")
            
            cg <- subset(onegene, V3 == "CG")
            if(nrow(cg) > minsite){
                outcg <- paste0(gsub("input", "CG", inputdf$file[i]), "_cg")
                write.table(cg[, c(1:2, 4:43)], outcg, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
                tem1 <- data.frame(file=outcg, sites=nrow(cg))
                write.table(tem1, paste0(outcg,".len"), sep="\t", row.names=FALSE, quote=FALSE)
            }
            
            chg <- subset(onegene, V3 == "CHG")
            if(nrow(chg) > minsite){
                outchg <- paste0(gsub("input", "CHG", inputdf$file[i]), "_chg")
                write.table(chg[, c(1:2, 4:43)], outchg, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
                tem2 <- data.frame(file=outchg, sites=nrow(chg))
                write.table(tem2, paste0(outchg,".len"), sep="\t", row.names=FALSE, quote=FALSE)
                
            }
            
            chh <- subset(onegene, V3 == "CHH")
            if(nrow(chh) > minsite){
                outchh <- paste0(gsub("input", "CHH", inputdf$file[i]), "_cg")
                write.table(chh[, c(1:2, 4:43)], outchh, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
                tem3 <- data.frame(file=outchh, sites=nrow(chh))
                write.table(tem3, paste0(outchh,".len"), sep="\t", row.names=FALSE, quote=FALSE)
            }
            if(verbose) print(i)
        }
    }
}    

##get command line args
    
options(echo=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
JOBID <- as.numeric(as.character(args[1]))


cmdno = 100
files <- list.files(path="largedata/Dm/region_input", full.names=TRUE)
inputdf <- data.frame(file=files, out="n")
tot <- ceiling(nrow(inputdf)/cmdno)

srow <- cmdno*(JOBID-1) + 1
if(JOBID == tot){
    erow <- nrow(inputdf)
}else{
    erow <- cmdno*JOBID
}

res <- C_format(inputdf, myr=srow:erow,  minsite=10, verbose=TRUE)

