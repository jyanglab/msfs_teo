### Jinliang Yang
### July 13th, 2016

### transform SNP to BCF 
library("farmeR")
cmd1 <- "cd largedata/gatk_vcf "
cmd2 <- "bgzip JRI20_joint_call.filtered_snps.vcf -@ 4"
cmd3 <- "tabix -p vcf JRI20_joint_call.filtered_snps.vcf.gz"
cmd4 <- "bcftools convert JRI20_joint_call.filtered_snps.vcf.gz -Ou -o JRI20_joint_call.filtered_snps.bcf"

set_farm_job(slurmsh = "largedata/GenSel/CL_test.sh",
             shcode = "sh largedata/myscript.sh", wd = NULL, jobid = "myjob",
             email = NULL, runinfo = c(TRUE, "bigmemh", "1"))

### transform InDels to BCF 
cmd1 <- "cd largedata/gatk_vcf "
cmd2 <- "bgzip JRI20_joint_call.filtered_indels.vcf -@ 8"
cmd3 <- "tabix -p vcf JRI20_joint_call.filtered_indels.vcf.gz"
cmd4 <- "bcftools convert JRI20_joint_call.filtered_indels.vcf.gz -Ou -o JRI20_joint_call.filtered_indels.bcf"
cmd5 <- "bcftools stats JRI20_joint_call.filtered_indels.bcf > JRI20_filtered_indels.stat"
cmd6 <- "plot-vcfstats --main-title InDels -p indel_plots/ JRI20_filtered_indels.stat"

set_farm_job(slurmsh = "slurm-script/vcf2bcf.sh",
             shcode = c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6), wd = NULL, jobid = "myjob",
             email = NULL, runinfo = c(TRUE, "bigmemm", "8"))


### check BCF summary and plot
cmd1 <- "bcftools stats JRI20_joint_call.filtered_indels.bcf > JRI20_filtered_indels.stat"
cmd2 <- "plot-vcfstats --main-title teo20 -p plotall/ JRI20_snps.stat JRI20_filtered_indels.stat"

