### Jinliang Yang
### 11-21-2016
### purpose: get SFS from genomic data


get_gene_meth <- function(gff_file="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff",
                          methra_file = "cache/stat_exon_mean_var.csv"){
    library(data.table)
    gff <- fread(gff_file, data.table=FALSE)
    names(gff) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")
    
    gene <- subset(gff, feature=="gene")
    gene$attribute <- gsub("ID=", "", gene$attribute)
    gene$attribute <- gsub(";.*", "", gene$attribute)
    
    ### extract mehtylation ratio of the population
    res <- read.csv(methra_file)
    #gbm <- subset(res, mm > cutoff)
    #ngbm <- subset(res, mm <= cutoff)
    
    gene <- merge(gene, res, by.x="attribute", by.y="geneid")
    
    out <- list()
    out[['gene']] <- gene
    #out[['exon']] <- exon
    #out[['intron']] <- intron
    return(out)
    
}

####

out <- get_gene_meth(gff_file="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff",
                     methra_file = "cache/stat_exon_mean_var.csv")
res <- out[[1]]
gbM <- subset(res, mm > 0.6)
ngbM <- subset(res, mm <= 0.6)

write.table(gbM[, c("seqname", "start", "end")], "largedata/gatk_vcf/gbm_gene_pos.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(ngbM[, c("seqname", "start", "end")], "largedata/gatk_vcf/non-gbm_gene_pos.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# bcftools filter --regions-file file -f '%CHROM\t%POS\t%REF\t%ALT[\t%IUPACGT]\n'  file.vcf.gz > bisnp.txt

### convert to PLINK
cmd1 <- "cd largedata/gatk_vcf"
cmd2 <- "bcftools index -f JRI20_bi_snps_annot.vcf.gz"
### extract regions of interest
cmd3 <- paste("bcftools filter --regions-file gbm_gene_pos.txt JRI20_filtered_snps_annot.bcf.gz",
              "-Ou -o JRI20_gbm.bcf")
cmd4 <- paste("bcftools filter --regions-file non-gbm_gene_pos.txt JRI20_filtered_snps_annot.bcf.gz",
              "-Ou -o JRI20_non-gbm.bcf")

### compute frq
cmd5 <- paste("plink -bcf JRI20_gbm.bcf --keep-allele-order --make-bed --out gbm_JRI20", 
              "--allow-extra-chr --freq")
cmd6 <- paste("plink -bcf JRI20_non-gbm.bcf --keep-allele-order --make-bed --out non-gbm_JRI20", 
              "--allow-extra-chr --freq")


set_farm_job(slurmsh = "slurm-script/bcf2plink.sh",
             shcode = c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6), wd = NULL, jobid = "bcf",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))


####################
library("data.table")
frq <- fread("largedata/gatk_vcf/gbm_JRI20.frq", header=TRUE)
frq1 <- subset(frq, NCHROBS == 40)
frq1$count <- frq1$NCHROBS*frq1$MAF
out1 <- table(frq1$count)

write.table(as.data.frame(out1), "cache/gbm_sfs_genomic.csv", sep=",", row.names=FALSE, quote=FALSE)


frq2 <- fread("largedata/gatk_vcf/non-gbm_JRI20.frq", header=TRUE)
frq2 <- subset(frq2, NCHROBS == 40)
frq2$count <- frq2$NCHROBS*frq2$MAF
out2 <- table(frq2$count)

write.table(as.data.frame(out2), "cache/non-gbm_sfs_genomic.csv", sep=",", row.names=FALSE, quote=FALSE)




