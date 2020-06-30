library(qtl)
d <- read.csv("Studer_2011_NatGenet_tb1_data-110706_rqtl.csv")
d <- d[-1:-2, ]

d$Year <- as.numeric(as.factor(d$Year))
d$ID <- as.numeric(as.factor(d$ID))
d$Family <- as.numeric(as.factor(d$Family))
#fit <- lm(TILL ~ Year*ID*Family, data=d)
x <- as.matrix(d[, c("Year", "ID","Family")])

tb1 <- read.cross("csv", file="Studer_2011_NatGenet_tb1_data-rqtl_input.csv",
                  genotype=c("MM","IM", "II"), alleles=c("I","M"), estimate.map=FALSE, na.strings=".", crosstype = "f2")


tb1 <- calc.genoprob(tb1, step=0, off.end=0, error.prob=0.01)
out.em <- scanone(tb1, pheno.col="TILL", method="em", addcovar=x)

plot(out.em, ylab="LOD score",xlim=c(270.35,270.56),ylim=c(150,210),xlab="Chr1 (Mb)",col="blue",cex.axis=1.2,cex.lab=1.2)
abline(v=c(270.495,270.513,270.5536),lty=2,lwd=2,col="red")

effectplot(tb1, mname1="GS1")
effectplot(tb1, mname1="GS2")
effectplot(tb1, mname1="GS3")
effectplot(tb1, mname1="GS4")
effectplot(tb1, mname1="GS6")

out2 <- scantwo(tb1)
tiff("twoD_interaction.tiff",res=600,units = "mm",height = 140,width = 160)
par(mar=c(2,5,1,5))

plot(out2,xlab="Chr1 (Mb)",ylab="Chr1 (Mb)")
dev.off()
summary(out2, allpairs=T)



library("lme4")
d <- read.csv("Studer_2011_NatGenet_tb1_data-110706.csv")
d[,10:24] <- apply(d[,10:24], 2, as.character)
d[d=="MM"] <- 0
d[d=="IM"] <- 1
d[d=="II"] <- 2
d$Year <- as.numeric(as.factor(d$Year))
d$Line <- as.numeric(as.factor(d$Line))
d$Family <- as.numeric(as.factor(d$Family))

d$TILL <- as.numeric(as.character(d$TILL))
m1 <- lmer(TILL ~ GS1 + GS8 + GS7 + GS4 + GS6 + GS5 + GS2 + GS3 + (1 | Year)+ (1 | Family)+ (1 |Line), data=d)
anova(m1)

fit1 <- lm(TILL ~ as.factor(Year) + GS1 + GS8 + GS7 + GS4 + GS6 + GS5 + GS2 + GS3 + GS6:GS2 + GS6:GS3 + GS3:GS4 + GS4:GS4+ GS2:GS3, data=d)

summary(fit1)


d$LBIL <- as.numeric(as.character(d$LBIL))

fit2 <- lm(LBIL ~ as.factor(Year) + GS1 + GS8 + GS7 + GS4 + GS6 + GS5 + GS2 + GS3 + GS6:GS2 + GS6:GS3 + GS5:GS2 + GS5:GS3, data=d)
summary(fit2)
