### Jinliang Yang
### July 13th, 2016

### transform VCF to PLINK
library("farmeR")

### annotate VCF file and remove unknown regions
cmd1 <- "cd largedata/gatk_vcf "
cmd2 <- paste("bcftools annotate --set-id +'%CHROM\\_%POS'", 
              "JRI20_joint_call.filtered_indels.vcf.gz",
              "-Ob -o JRI20_filtered_indel_annot.bcf.gz")
cmd3 <- "tabix JRI20_filtered_indel_annot.bcf.gz"

set_farm_job(slurmsh = "slurm-script/bcf_annot_indel.sh",
             shcode = c(cmd1, cmd2, cmd3), wd = NULL, jobid = "bcf",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))

### convert to PLINK
cmd1 <- "cd largedata/gatk_vcf"
cmd2 <- paste("plink -bcf JRI20_filtered_indel_annot.bcf.gz --keep-allele-order --make-bed --out JRI20_indel", 
              "--allow-extra-chr --freq --missing --het --ibc")
set_farm_job(slurmsh = "slurm-script/bcf2plink.sh",
             shcode = c(cmd1, cmd2), wd = NULL, jobid = "bcf",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))

### plot
library(data.table)

miss <- fread("largedata/gatk_vcf/JRI20.lmiss")
frq <- fread("largedata/gatk_vcf/JRI20.frq")
#G0	Missing genotype count (so C1 + C2 + 2 * G0 is constant on autosomal variants)

pdf("graphs/teo20_miss_maf.pdf", height=5, width=10)
par(mfrow=c(1, 2))
hist(miss$F_MISS, main="Missing rate (N=20)", xlab="missing", col="#cdaa7d")
hist(frq$MAF, main="SFS (N=20)", xlab="Non-ref allele frq", col="#cdaa7d")
dev.off()



