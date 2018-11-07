### Jinliang Yang
### 02-03-2017
### purpose: plot with the summarized methylation levels

##get command line args

pwd <- list.files(path=c("largedata/lcache"),  pattern="chr1_fea", full.names = TRUE)
length(pwd)
out <- data.frame()
for(i in 1:length(pwd)){
    tem <- read.csv(pwd[i])
    tem$fid <- pwd[i]
    out <- rbind(out, tem)
}
out$context <- gsub("cache\\/|_.*", "", out$fid)
write.table(out, "cache/chrp_feas.csv", sep=",", row.names=FALSE)


#############
d <- read.csv("cache/chrp_feas.csv")
d$x <- as.numeric(as.character(gsub("V", "", d$var)))
d$feature <- as.character(d$feature)

d$context <- gsub("largedata\\/|\\/JR.*", "", d$file)
d[d$context == "COMET",]$context <- "CG"
d$context <- gsub("COMET_", "", d$context)



d[d$feature == "exon1st", ]$x <- d[d$feature == "exon1st", ]$x + 10
d[d$feature == "intron1st", ]$x <- d[d$feature == "intron1st", ]$x + 20
d[d$feature == "intronlast", ]$x <- d[d$feature == "intronlast", ]$x + 30
d[d$feature == "exonlast", ]$x <- d[d$feature == "exonlast", ]$x + 40
d[d$feature == "down1k", ]$x <- d[d$feature == "down1k", ]$x + 50

d <- subset(d, !(feature %in% "gene"))






library(ggplot2)

plot_eff <- function(){
    
    #######
    theme_set(theme_grey(base_size = 18)) 
    
    fsize=18
    p1 <- ggplot(d, aes(x=x, y=mc, colour=factor(context)) )+
        labs(colour="context") +
        theme_bw() +
        xlab("") +
        ylab("Methylation Level") +
        scale_color_manual(values=c("#8b2323", "#E69F00", "#56B4E9")) +
        #scale_color_manual(values=context) +
        #scale_linetype_manual(values=lty1) +
        guides(size=FALSE) +
        geom_smooth(method="loess", span = 0.08, size=2) +
        #geom_smooth(span = 0.1) +
        geom_vline(xintercept=c(10, 20, 30, 40, 50), col="black") +
        #facet_grid(~ feature) +
        scale_x_continuous(breaks=c(5,15,25,35, 45, 55),
                           labels=c("Upstream 1kb", "First Exon",
                                  "First Intron", "Last Intron",
                                  "Last Exon", "Downstream 1kb")) +
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text.x = element_text(angle = 20, hjust = 0.7),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize) )
    #abline(v=c(10, 20, 30, 40, 50), col="black")
    p1
    
}

########
pdf("graphs/Fig2b_var.pdf", width=8, height =5)
plot_eff()
dev.off()
