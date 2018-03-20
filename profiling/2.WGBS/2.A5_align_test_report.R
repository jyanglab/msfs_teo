### Jinliang Yang
### July 15th, 2016

library(farmeR)

files <- list.files(path = "largedata/bismark", pattern = "PE_report.txt", full.names = TRUE)
features <- c("Sequence pairs analysed in total:\t",
              "Number of paired-end alignments with a unique best hit:\t",
              "Sequence pairs with no alignments under any condition:\t",
              "Sequence pairs did not map uniquely:\t",
              #"Sequence pairs which were discarded because genomic sequence could not be extracted:",
              "Total number of C's analysed:\t",
              "Total methylated C's in CpG context:\t",
              "Total methylated C's in CHG context:\t",
              "Total methylated C's in CHH context:\t",
              "Total unmethylated C's in CpG context:\t",
              "Total unmethylated C's in CHG context:\t",
              "Total unmethylated C's in CHH context:\t"
              )
res <- get_file2tab(files, features, replace=T )

names(res) <- c("totpe", "hit1", "hit0", "hitN", "totC","mCG", "mCHG", "mCHH", "unCG", "unCHG", "unCHH")
row.names(res) <- gsub(".*/|_PE_.*", "", files) 

res <- as.data.frame(apply(res, 2, as.numeric))
res$seqdp <- (res$totpe)*200/(2500*1e6)
res$mr1 <- with(res, round(hit1/totpe,3)) #uniquely mapped rate
res$mr0 <- with(res, round(hit0/totpe,3)) #non-mapping rate
res$mrN <- with(res, round(hitN/totpe,3)) #multiple-mapping rate

res$cg <- with(res, round(mCG/(mCG+unCG), 3))
res$chg <- with(res, round(mCHG/(mCHG+unCHG), 3))
res$chh <- with(res, round(mCHH/(mCG+unCHH), 3))
res$seqid <- gsub(".*/|_PE_.*", "", files) 

write.table(res, "cache/comp2_bismap_stat.csv", sep=",", row.names=FALSE, quote=FALSE)

#########
library("data.table")

bi <- read.csv("cache/comp2_bismap_stat.csv")







