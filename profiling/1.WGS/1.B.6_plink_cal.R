### Jinliang Yang
### July 15th, 2016

### calculate other things

cmd1 <- "cd largedata/gatk_vcf"
cmd2 <- "plink --bfile JRI20_snp --threads 6 --allow-extra-chr --out JRI20_snp --genome --homozyg"

set_farm_job(slurmsh = "slurm-script/plink_stat.sh",
             shcode = c(cmd1, cmd2), wd = NULL, jobid = "pstat",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "16"))

## ROH: detected no runs of homozygocity
cmd <- "bcftools roh JRI20_filtered_snps_annot.bcf.gz"
cmd <- "plink --bfile JRI20_snp --threads 4 --allow-extra-chr --maf 0.05 --homozyg-snp 50 --homozyg-kb 500 --out JRI20_snp"


### inbreeding
# small sample, inaccurate!
# F coefficient estimates (i.e. ([observed hom. count] - [expected count]) / ([total observations] - [expected count]))
het <- read.table("largedata/gatk_vcf/JRI20_snp.het", header=T)
ibc <- read.table("largedata/gatk_vcf/JRI20_snp.ibc", header=T)

### Handy

cmd <- "plink --bfile JRI20_snp --threads 6 --allow-extra-chr --out JRI20_snp --hardy midp"

hwe <- fread("largedata/gatk_vcf/JRI20_snp.hwe")
res <- hwe[P < 1e-5]
write.table(res, "cache/JRI20_hwe_pe5.csv", sep=",", row.names=FALSE, quote=FALSE)


res <- read.csv("cache/JRI20_hwe_pe5.csv")
res <- subset(res, P < 1e-6)
res$log10p <- -log10(res$P)

res <- subset(res, CHR %in% 1:10)
res$chr <- as.numeric(as.character(res$CHR))
res$pos <- as.numeric(as.character(gsub(".*_", "", res$SNP)))

plot_mht(res = res, cex = 0.9, pch = 16, col = rep(c("slateblue","cyan4"), 5), 
         GAP = 5e+06, yaxis = NULL, col2plot = "log10p")


### identify Beagle IBD
shcode = c("module load java", 
           "cd largedata/gatk_vcf/",
           "java -Xmx60g -jar ~/bin/beagle.22Feb16.8ef.jar gt=JRI20_filtered_snps_annot.vcf.gz ibd=true out=JRI20_snp")
set_array_job(shid = "slurm-script/run_beagle.sh",
              shcode = shcode, arrayjobs = "1", wd = NULL,
              jobid = "getibd", email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))

