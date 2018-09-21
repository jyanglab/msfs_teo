run_mcmc_bygeneset <- function(JOBID, typefile="largedata/type.csv", geneset, cutoff, outid){
    ## geneset: [data.frame], cols: geneid, value
    ## cutoff: cutoff for the values
    ### outid: prefix for outid
    
    ###### main codes:    
    a <- read.csv(typefile)
    mya <- a[JOBID, ]
    
    ### type
    if(mya$type == "CG"){
        df <- fread("largedata/lcache/SFS_comet_blocks_CG.csv", data.table=FALSE)
    }else if(mya$type == "CHG"){
        df <- fread("largedata/lcache/SFS_comet_blocks_CHG.csv", data.table=FALSE)
    }else if(mya$type == "CHH"){
        df <- fread("largedata/lcache/SFS_comet_blocks_CHH.csv", data.table=FALSE)
    }else{
        stop("### type error!")
    }
    df$chr <- gsub("chr", "", df$chr)
    df$cid <- paste(df$chr, df$bid, sep="_")
    
    ### length quantile
    qt <- quantile(df$length)
    if(mya$length == 1){
        df <- subset(df, length <= qt[2])
    }else if(mya$length == 2){
        df <- subset(df, length > qt[2] & length <= qt[3])
    }else if(mya$length == 3){
        df <- subset(df, length > qt[3] & length <= qt[4])
    }else if(mya$length == 4){
        df <- subset(df, length > qt[4])
    }else if(mya$length == 0){
        df <- df
    }
    
    
    if(mya$TE == "yes"){
        repeats <- fread("~/dbcenter/AGP/AGPv2/repeats/ZmB73_5a_MIPS_repeats.gff", header=FALSE, data.table=FALSE)
        names(repeats) <- c("seqname", "source", "repeat", "start", "end", "score",
                            "strand", "frame", "attribute")
        
        repeats$feature <- gsub(".*type=|;name=.*", "", repeats$attribute)
        gff <- repeats
    }else if(mya$TE == "no"){
        
        if(!is.null(mya$gset)){
            ######## format GFF and repeat files
            gff <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=TRUE, data.table=FALSE)
            names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                            "strand", "frame", "attribute")
            
            if(mya$gset == "above"){
                #res <- read.csv("cache/stat_exon_mean_var.csv")
                gff$geneid <- gsub(";.*|_.*", "", gff$attribute)
                gff$geneid <- gsub(".*=", "", gff$geneid)
                ####
                g1 <- subset(geneset, value > cutoff)
                gff <- subset(gff, geneid %in% g1$geneid)
                
            }else if(mya$gset == "below"){
                #res <- read.csv("cache/stat_exon_mean_var.csv")
                gff$geneid <- gsub(";.*|_.*", "", gff$attribute)
                gff$geneid <- gsub(".*=", "", gff$geneid)
                ####
                g2 <- subset(geneset, value <= cutoff)
                gff <- subset(gff, geneid %in% g2$geneid)
            }
        }
    }
    
    ##########
    out <- run_overlap_MCMC(df, gff, fea=mya$fs, runid=mya$id, outdir="largedata/lcache/", outid=outid)
    
}