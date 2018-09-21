runmcmc_te_geneset <- function(JOBID, inputdf="largedata/type.csv"){
    ## JOBID: control the array job number. [num, 1]
    ## inputdf: [data.frame, cols=comet_file:(comet blocks file, i.e. largedata/lcache/SFS_comet_blocks_CG.csv), 
    #                    length:(blocks quantile len, 0,1,2,3,4), TE:(yes, no), 
    #                    outRD:(chr, "out.RData"),
    #                    geneset_file:(chr, csv, cols: geneid, value, i.e. "geneset.tex"),
    #                    
    #                    optional:
    #                    feature:(chr, "exon", "intron", "up1k", "gene"),
    #                    gset:(chr, "below" or "above", "wholeset"),
    #                    cutoff:(num, num, 0.5)]
    ## geneset: gene set to determine which set to cal sfs. [data.frame, cols: geneid, value]
    ## cutoff: cutoff for the values to get up set and down set of the gene. [num, 0.5]
    
    ###### main codes:
    a <- read.csv(inputdf, header=TRUE)
    mya <- a[JOBID, ]
    
    outfile <- as.character(mya$outRD)
    if(!is.null(mya$geneset_file)){
        geneset <- read.csv(as.character(mya$geneset_file), header=TRUE)
    }
    if(!is.null(mya$cutoff)){
        cutoff <- mya$cutoff
    }
    if(!is.null(mya$gset)){
        gset <- as.character(mya$gset)
    }
    
    
    ### type: comet context
    #df <- fread("largedata/lcache/SFS_comet_blocks_CG.csv", data.table=FALSE)
    df <- fread(as.character(mya$comet_file), data.table=FALSE)
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
    
    out <- get_overlap_sfs(df, gff, fea=mya$feature)
    
    
    ##########
    res <- MCMCBC(my_sfs=out$Freq, rates=c(1E8,1E8,1E8), sd=c(0.05,0.05,0.05), k=0:40,
                  conditional=FALSE, Ne=150000, ngen=100000, verbose=FALSE)
    d <- mplot(res, burnin=0.1, rates=c(1E8,1E8,1E8), doplot=FALSE)
    # posterior mu [ 2.03493482965863e-06 ], nu [ 1.98123661448284e-07 ] and s [ 1.89713888507979e-05 ]
    save(list=c("d", "res"), file=outfile)
}