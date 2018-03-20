### Jinliang Yang
### check the quality of WGBS data
### 3/22/2016

library("farmeR")

### list files and run QC
fq1 <- list.files(path="largedata/wgbs_fq", pattern="R1.fastq.gz$", full.names=T)
fq2 <- list.files(path="largedata/wgbs_fq", pattern="R2.fastq.gz$", full.names=T)


inputdf <- data.frame(fq1=fq1, fq2=fq2, out1= gsub("_R1", "_R1_trimmed",  fq1), 
                      out2=gsub("_R2", "_R2_trimmed",  fq2))
run_cutadapt(inputdf[2:20, ], ad1="AGATCGGAAGAGC", ad2="AGATCGGAAGAGC", q=20,
             email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 1))


###########################
files <- list.files(path = "largedata/wgbs_fq", pattern = "testout", full.names = TRUE)
features <- c("Command line parameters: -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ",
              "Total read pairs processed:",
              "Read 1 with adapter:",
              "Read 2 with adapter:",
              "Pairs written",
              "Total basepairs processed:",
              "Quality-trimmed:",
              "Total written"
)
res <- get_file2tab(files, features, replace=T)

row.names(res) <- 1:nrow(res)
names(res) <- c("files", "totpe", "r1ap", "r2ap", "outpe","totbp", "qtrimmed", "totout")
res$files <- gsub(".*/", "", res$files)

res$outpe <- gsub(".*\\)\\:", "", res$outpe)
res$totout <- gsub(".*\\)\\:", "", res$totout)
write.table(res, "cache/maize_bismap_trimming_stat.txt", sep="\t", row.names=FALSE, quote=FALSE)

######## methylation
stat <- read.delim("cache/maize_bismap_trimming_stat.txt", header=TRUE)
#stat$totpe <- gsub("\\s+|,|\\(.*|bp", "", stat$totpe)

res <- apply(stat[,-1], 2, function(x) gsub("\\s+|,|\\(.*|bp", "",x))
res <- as.data.frame(apply(res, 2, as.numeric))
res$files <- gsub("_.*", "", stat$files)

res$aptrimmed <- res$totbp - res$totout - res$qtrimmed

library(ggplot2)
library(tidyr)

lres <- gather(res[, 6:9], type, reads, c(1,2,4))
lres <- lres[order(lres$files, lres$type, decreasing = TRUE),]

lres$type <- factor(lres$type, levels = c("aptrimmed", "qtrimmed",  "totout"), 
                    labels=c("Adapter", "Q<20", "remaning"), ordered=TRUE)

theme_set(theme_grey(base_size = 18)) 
s <- ggplot(lres, aes(x=files, y=reads, fill = type)) + 
    #opts(axis.text.x=theme_text(angle=90)) +
    geom_bar(stat="identity") +
    labs(x="Accession ID", y="Base-pairs", fill="Trimming") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

pdf("graphs/SFig_wgbs_trimming.pdf", width=10, height=5)
s
dev.off()
