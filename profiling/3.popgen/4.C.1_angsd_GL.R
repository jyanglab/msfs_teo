### Jinliang Yang
### June 30, 2016

library(farmeR)

bamfile <- list.files(path="/group/jrigrp4/teosinte-parents/20parents/bams",
                      pattern="sorted.*bam$", full.names=TRUE)
write.table(bamfile, "largedata/bamlist.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### got Genotype likelihood
run_ANGSD(type = "GL", shfile = "slurm-script/run_angsd_gl.sh",
          bamlist = "largedata/bamlist.txt", outfile = "largedata/teo20", 
          ref = "$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa", 
          anc = "$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa", 
          email = "yangjl0930@gmail.com",
          runinfo = c(TRUE, "bigmemh", 16))

#ngsF --n_ind 20 --n_sites 10000000 --glf teo20.glf --out teo20 --approx_EM --fast_lkl

### do theta
run_ANGSD(type = "theta", shfile = "slurm-script/run_angsd_theta.sh",
          bamlist = "largedata/bamlist.txt", outfile = "largedata/teo20", 
          ref = "$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa", 
          anc = "$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
          nInd=20, minInd=5,
          email = "yangjl0930@gmail.com",
          runinfo = c(TRUE, "bigmemm", 32))