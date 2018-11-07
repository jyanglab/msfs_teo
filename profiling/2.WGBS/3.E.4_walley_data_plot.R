### Jinliang Yang
### 02-03-2017
### purpose: plot with the summarized methylation levels

##get command line args

pwd <- list.files(path=c("largedata/lcache"),  pattern="chr1_exp", full.names = TRUE)
length(pwd)
out <- data.frame()
for(i in 1:length(pwd)){
    tem <- read.csv(pwd[i])
    tem$fid <- pwd[i]
    out <- rbind(out, tem)
}
out$context <- gsub(".*lcache\\/|_.*", "", out$fid)
out$type <- gsub(".*_chr1_exp_|_.*", "", out$fid)
table(out$context)

write.table(out, "cache/chrp_exp.csv", sep=",", row.names=FALSE)


#############
d <- read.csv("cache/chrp_exp.csv")
d$x <- as.numeric(as.character(gsub("V", "", d$var)))
d$feature <- as.character(d$feature)

d$status <- gsub(".*RNA_|.*protein_|_.*", "", d$fid)


table(d$status)

d[d$feature == "down1k", ]$x <- d[d$feature == "down1k", ]$x + 10
d[d$feature == "up1k", ]$x <- d[d$feature == "up1k", ]$x - 10






library(ggplot2)

plot_eff <- function(outfile, getpdf){
    
    #######
    theme_set(theme_grey(base_size = 18)) 
    
    fsize=18
    p1 <- ggplot(subset(d, type %in% "RNA"), aes(x=x, y=mc, colour=factor(context)) )+
        labs(colour="context") +
        theme_bw() +
        xlab("") +
        ylab("Methylation Level") +
        #scale_color_manual(values=context) +
        #scale_linetype_manual(values=lty1) +
        guides(size=FALSE) +
        geom_smooth(method="loess", span = 0.08, size=2) +
        #geom_smooth(span = 0.1) +
        geom_vline(xintercept=c(1, 10), col="black") +
        facet_grid(~ status) +
        scale_color_manual(values=c("#8b2323", "#E69F00", "#56B4E9")) +
        #ylim(0, 0.05) +
        scale_x_continuous(breaks=c(-5, 5, 15),
                           labels=c("Upstream 1kb", "Expression", "Downstream 1kb")) +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
              axis.text.x = element_text(angle = 20, hjust = 0.7),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize) )
    #abline(v=c(10, 20, 30, 40, 50), col="black")
    p1
    
    p1 <- ggplot(subset(d, type %in% "protein"), aes(x=x, y=mc, colour=factor(context)) )+
        labs(colour="context") +
        theme_bw() +
        xlab("") +
        ylab("Methylation Level") +
        #scale_color_manual(values=context) +
        #scale_linetype_manual(values=lty1) +
        guides(size=FALSE) +
        geom_smooth(method="loess", span = 0.08, size=2) +
        #geom_smooth(span = 0.1) +
        geom_vline(xintercept=c(1, 10), col="black") +
        facet_grid(~ status) +
        scale_color_manual(values=c("#8b2323", "#E69F00", "#56B4E9")) +
        #ylim(0, 0.05) +
        scale_x_continuous(breaks=c(-5, 5, 15),
                           labels=c("Upstream 1kb", "Expression", "Downstream 1kb")) +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
              axis.text.x = element_text(angle = 20, hjust = 0.7),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize) )
    #abline(v=c(10, 20, 30, 40, 50), col="black")
    p1
    
}

########
p <- plot_eff(outfile="graphs/Fig2b_var.pdf", getpdf)
p
