### Jinliang Yang
### 10-16-2016
### purpose: get SFS for features


get_gff_features <- function(gff_file="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff"){
    library(data.table)
    gff <- fread(gff_file, data.table=FALSE)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    ### get canonical genes on chr1:10
    info <- fread("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt", data.table=FALSE)
    info <- subset(info, chromosome %in% paste0("chr", 1:10) & is_canonical == "yes" )
    message(sprintf("###>>> Number of canonical transcripts on chr1-10: [ %s ]", nrow(info)))
    
    gene <- subset(gff, feature=="gene")
    gene$attribute <- gsub("ID=", "", gene$attribute)
    gene$attribute <- gsub(";.*", "", gene$attribute)
    gene <- subset(gene, attribute %in% info$gene_id)
    
    
    exon <- subset(gff, feature=="exon")
    exon$attribute <- gsub("Parent=", "", exon$attribute)
    exon$attribute <- gsub(";.*", "", exon$attribute)
    exon <- subset(exon, attribute %in% info$transcript_id)
    
    intron <- subset(gff, feature=="intron")
    intron$attribute <- gsub("Parent=", "", intron$attribute)
    intron$attribute <- gsub(";.*", "", intron$attribute)
    intron <- subset(intron, attribute %in% info$transcript_id)
    
    out <- list()
    out[['gene']] <- gene
    out[['exon']] <- exon
    out[['intron']] <- intron
    return(out)
    
}

#######
f <- get_gff_features(gff_file="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff")

save(list="f", file="largedata/AGPv2_features.RData")



