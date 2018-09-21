
re_calc <- function(myfile){
    ### First, SNP table:
    h <- read.table("largedata/gatk_vcf/JRI20_bi_snps_annot.header", header=T)
    ### I do not need to merge because the orders are exactly the same.
    idtab <- read.csv("data/teo20_ids.csv")
    n <- c("chr", "pos", "ref", "alt", paste0("JR", idtab$plate))
    snpdt <- fread("largedata/gatk_vcf/JRI20_bi_snps_annot.txt")
    names(snpdt) <- n
    
    res <- data.frame()
    dt <- fread(myfile)
    ### loop through chr and then 
    for(j in 1:10){
        sid <- gsub(".*\\/|_pe.*", "", myfile)
        message(sprintf("###>>> re-cal sample [ ID=%s ] chr [ %s ]", sid, j))
        chr <- dt[V1 == j]
        chr[, snpid := paste(V1, V2, sep="_")]
        
        
        sub <- snpdt[, c("chr", "pos", "ref", "alt", sid),  with=FALSE]
        sub <- sub[chr == j]
        names(sub)[5] <- "SAMPLE"
        sub <- sub[SAMPLE == "Y" | SAMPLE == "R"]
        #sub$snpid <- paste(sub$chr, sub$pos, sep="_")
        sub[, snpid := paste(chr, pos, sep="_")]
        
        #V4: methylated C, V5: unmethylated C
        # for C/T and G/A sites, we simply assume half of the observed counts were from the T variant. 
        # Therefore, we divided the total observed counts by 2. 
        out <- merge(chr, sub, by.x="snpid", by.y="snpid", all.x=TRUE)
        out[, tot := V4 + V5]
        out[!is.na(SAMPLE), tot := ceiling((V4 + V5)/2)]
        message(sprintf("###>>> summary: YR sites/tot [%s], ratio >1 [ %s ]", 
                        round(nrow(sub)/nrow(out), 3), 
                        round(nrow(out[V4 > tot & (SAMPLE == "Y" | SAMPLE == "R"),])/nrow(out), 3)))
        
        tem <- data.frame(file=sid, chr=j, YRs=nrow(sub), tot=nrow(out))
        res <- rbind(res, tem)
        
        out[V4 > tot & (SAMPLE == "Y" | SAMPLE == "R"), tot := V4]
        
        out[order(V2)]
        outf <- paste0(myfile, ".rc")
        message(sprintf("###>>> writing [ chr%s ] to [ %s ] ...", j, outf))
        if(j == 1){
            write.table(out[, c("snpid", "V1", "V2", "V3", "V4", "V5", "tot"), with=FALSE], 
                        outf, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE)
        }else{
            write.table(out[, c("snpid", "V1", "V2", "V3", "V4", "V5", "tot"), with=FALSE], 
                        outf, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
        }
    }
    return(res)
}