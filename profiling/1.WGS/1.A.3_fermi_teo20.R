### Jinliang Yang
### use Li Heng's De novo assembly based variant calling pipeline
### 3/22/2016

fq <- list.files(path="/group/jrigrp/teosinte-parents/seq-merged", pattern="fastq.gz$", full.names = TRUE)


for(i in 1:length(fq)){
    shid <- paste0("slurm-script/run_fqchk_", i, ".sh")
    command <- paste0("seqtk fqchk -q 20 ", fq[i], " > ", paste0("largedata/fqchk/", gsub(".*/", "", fq[i]), ".qc"))
    cat(command, file=shid, sep="\n", append=FALSE)
}
shcode <- paste("sh slurm-script/run_fqchk_$SLURM_ARRAY_TASK_ID.sh", sep="\n")

set_array_job2(shid="slurm-script/run_fqchk_run.sh", shcode=shcode,
              arrayjobs="1-40", wd=NULL, jobid="myjob", email=NULL,
              run = c(TRUE, "med", "5200", "2"))


fq1 <- list.files(path="/group/jrigrp/teosinte-parents/seq-merged", pattern="_1.fastq.gz$", full.names = TRUE)
fq2 <- list.files(path="/group/jrigrp/teosinte-parents/seq-merged", pattern="_2.fastq.gz$", full.names = TRUE)

fq <- data.frame(fq1=fq1, fq2=fq2, out="",
                 sample=gsub(".*Sample_|_index.*", "", fq1), sample2=gsub(".*Sample_|_index.*", "", fq2))
fq$out <- paste0("$HOME/Documents/Github/methylation/largedata/fermi/", fq$sample)

library(farmeR)
run_fermikit(fq, kitpath="$HOME/bin/fermikit/fermi.kit",
             ref.fa="$HOME/dbcenter/AGP/AGPv2/Zea_mays.AGPv2.14.dna.toplevel.fa",
             s='2.3g', l=100,
             email="yangjl0930@gmail.com", 
             runinfo = c(TRUE, "bigmemm", 16))


