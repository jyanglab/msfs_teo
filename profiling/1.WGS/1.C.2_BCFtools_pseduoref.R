### Jinliang Yang
### July 13th, 2016


library(farmeR)

### test
"samtools faidx $HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa 1:1-10000 | bcftools consensus -i JRI20_joint_call.filtered_snps.vcf.gz -o out1.fa"
"samtools faidx $HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa 1:1-10000 | bcftools consensus JRI20_joint_call.filtered_snps.vcf.gz -o out2.fa"
"samtools faidx $HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa 1:1-10000 > out0.fa"

"bcftools consensus -i -s NA001 -f in.fa in.vcf.gz > out.fa"


cmd <- paste("bcftools consensus -i -f $HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
             "largedata/gatk_vcf/JRI20_joint_call.filtered_snps.vcf.gz", 
             "-o largedata/pgenome/teo20.fa")

set_farm_job(slurmsh = "slurm-script/bcf2cons.sh",
             shcode = cmd, wd = NULL, jobid = "cns",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))

### checking results
"samtools faidx teo20.fa 1:1-10000"


