### Jinliang Yang
### Gene-body Methylation

res <- read.csv("cache/stat_exon_mean_var.csv")

hist(res$mm, breaks=50, main="Avg. Levels of gbM", col="#cdb79e",
     xlab="Methylation Levels")
abline(v = 0.6, lty=2, lwd=3)

subset(res, mm > 0.6) %>% nrow
subset(res, mm <= 0.6) %>% nrow

####
gbM <- subset(res, mm > 0.6)
ngbM <- subset(res, mm <= 0.6)





res <- MCMCBC(my_sfs, sites, rates, sd, k, burnin,
                   conditional=FALSE, Ne, ngen,
                   verbose=TRUE)
    