### Jinliang Yang
### August 8th, 2016


### filter to biallelic loci
# bcftools view JRI20_filtered_snps_annot.bcf.gz -m2 -M2 -v snps -Oz -o JRI20_bi_snps_annot.vcf.gz
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%IUPACGT]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.txt
# bcftools query -f 'chr\tpos\tref\talt[\t%SAMPLE]\n' JRI20_bi_snps_annot.vcf.gz > JRI20_bi_snps_annot.header

library("data.table")
h <- read.table("largedata/gatk_vcf/JRI20_bi_snps_annot.header", header=T)
snpdt <- fread("largedata/gatk_vcf/JRI20_bi_snps_annot.txt")
names(snpdt) <- names(h)

library("Biostrings")
fa <- readDNAStringSet(filepath = "~/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa", format="fasta")
#bck <- fa
uniqueLetters(fa)
#[1] "A" "C" "G" "T" "S" "Y" "N"
# alphabetFrequency(fa)
# length of each chromosome
width(fa)
names(fa)

####### make pseudo-ref using an R package
library(pseudoRef)
myr <- data.frame(from=c("M", "Y", "R", "K"), to=c("C", "C", "G", "G"))
out <- pseudoRef(fa, snpdt, sidx=5:ncol(snpdt), arules=myr, outdir="largedata/wgbs_pgen")
 
getsum <- function(x){
    apply(x, 2, sum)
}   
out2 <- sapply(out, getsum)  
write.table(out2, "cache/pgen_report.csv", sep=",", quote=FALSE)    


####### make pseudo-ref IUPAC
out <- pseudoRef(fa, snpdt, sidx=5:ncol(snpdt), arules=NULL, 
                 outdir="/group/jrigrp7/jyang/methylation/largedata/pgen")

geneid <- read.table("/home/jolyang/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt", header=T)
#geneid <- subset(geneid, is_canonical == "yes")

g <- c("GRMZM2G161472","GRMZM2G067591", "GRMZM2G081554","GRMZM2G044481","GRMZM2G068808",
       "AC218998.2FG011","GRMZM2G391312","AC214360.3FG001","GRMZM2G093603","GRMZM2G016922",
       "GRMZM2G093526","GRMZM2G028306","GRMZM2G127087")

gid <- subset(geneid, gene_id %in% g)


