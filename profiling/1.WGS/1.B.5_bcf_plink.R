### Jinliang Yang
### July 13th, 2016

### transform BCF to PLINK

library("farmeR")
cmd1 <- "cd largedata/gatk_vcf "
cmd2 <- "bgzip JRI20_joint_call.filtered_snps.vcf -@ 4"
cmd3 <- "tabix -p vcf JRI20_joint_call.filtered_snps.vcf.gz"
cmd4 <- "bcftools convert JRI20_joint_call.filtered_snps.vcf.gz -Ou -o JRI20_joint_call.filtered_snps.bcf"

set_farm_job(slurmsh = "largedata/GenSel/CL_test.sh",
             shcode = "sh largedata/myscript.sh", wd = NULL, jobid = "myjob",
             email = NULL, runinfo = c(TRUE, "bigmemh", "1"))

### annotate VCF file and remove unknown regions
cmd1 <- "cd largedata/gatk_vcf "
#"-e 'CHROM ~ \"UNKNOWN\" | CHROM ~ \"mitochondrion\" | CHROM ~ \"chloroplast\" '",
cmd2 <- paste("bcftools annotate --set-id +'%CHROM\\_%POS'", 
              "JRI20_joint_call.filtered_snps.bcf",
              "-Ob -o JRI20_filtered_snps_annot.bcf.gz")
#cmd3 <- "bgzip JRI20_filtered_snps_annot.bcf -@ 8"
cmd3 <- "tabix JRI20_filtered_snps_annot.bcf.gz"

set_farm_job(slurmsh = "slurm-script/bcf_annot.sh",
             shcode = c(cmd1, cmd2, cmd3), wd = NULL, jobid = "bcf",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))


### convert to PLINK
cmd1 <- "cd largedata/gatk_vcf"
cmd2 <- paste("plink -bcf JRI20_filtered_snps_annot.bcf.gz --keep-allele-order --make-bed --out JRI20", 
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



