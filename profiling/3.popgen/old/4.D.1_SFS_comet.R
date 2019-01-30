### Jinliang Yang
### 10-16-2016
### purpose: get SFS for features


library("data.table")
get_comet_sfs <- function(cometfile="largedata/COMET/CG_COMET/comet_blocks.csv", 
                          outfile="cache/sfs_cg_comet.csv",
                          cols=3:22){
    comet <- fread(cometfile, data.table=FALSE)
    
    comet$miss <- apply(comet[, cols], 1, function(x){
        return(sum(x < 0))
    })
    comet <- subset(comet, miss == 0)
    ###########
    f <- apply(comet[, cols], 1, function(x){
        n2 <- sum(x == 2)
        n1 <- sum(x == 1)
        return(2*n2 + n1)
    })
    sfs <- table(f)
    write.table(sfs, outfile, sep=",", row.names=FALSE, quote=FALSE)
    return(sfs)
}

###############
files <- list.files(path="largedata/COMET/CG_COMET", pattern="chr10_comet", full.names = TRUE)
for(i in 1:length(files)){
    outf <- gsub(".*/", "cache/sfs_", files[i])
    sfs <- get_comet_sfs(cometfile=files[i], 
                         outfile=outf,
                         cols=3:22)
}
    
############
sfs_files <- list.files(path="cache", pattern="sfs_chr10_comet", full.names = TRUE)
par(mfrow=c(3, 3))
for(i in 1:length(files)){
    sfs <- read.csv(sfs_files[i])
    f <- gsub(".*_|.csv", "", sfs_files[i])
    plot(sfs$Freq, main=f)
}


