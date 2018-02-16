### Jinliang Yang
### use Li Heng's De novo assembly based variant calling pipeline
### 3/22/2016

fq1 <- list.files(path="~/dbcenter/BMfastq", pattern="_1.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="~/dbcenter/BMfastq", pattern="_2.fastq.gz$", full.names = TRUE)

sampleid <- read.table("~/dbcenter/BMfastq/sampleid.txt", header=TRUE)

bam <- list.files(path="~/dbcenter/BMfastq/bam", pattern="bam$", full.names = TRUE)
inputdf <- data.frame(fq1=fq1[c(1,7)], fq2=fq2[c(1,7)], out="mysample",
                      group=c("g1", "g2"), sample=c("Mo17", "B73"), PL="illumina", 
                      LB=c("lib1", "lib2"), PU=c("unit1", "unit2"),
                      bam=bam[c(1,3)])
inputdf$out <- gsub(".sra.*", "", inputdf$fq1)
inputdf$out <- gsub("BMfastq", "BMfastq/bam", inputdf$out)

###########
library(farmeR)
run_GATK(inputdf, 
         ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
         gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
         picardpwd="$HOME/bin/picard-tools-2.1.1/picard.jar",
         minscore = 5, markDup=TRUE, addRG=TRUE,
         realignInDels=FALSE, indels.vcf="indels.vcf",
         recalBases=FALSE, dbsnp.vcf="dbsnp.vcf", 
         email="yangjl0930@gmail.com",
         runinfo = c(FALSE, "bigmemm", 16))

####### joint variant calling
gvcf <- list.files(path="~/dbcenter/BMfastq/bam", pattern="g.vcf$", full.names = TRUE)
outvcf <- "~/dbcenter/BMfastq/bam/joint_call.vcf"
run_GATK_JointGenotype(
    gvcf, outvcf,
    ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
    gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
    includeNonVariantSites = FALSE,
    hardfilter= TRUE,
    snpflt="\"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"",
    indelflt="\"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"",
    email="yangjl0930@gmail.com",
    runinfo = c(TRUE, "bigmemh", 4)
)


