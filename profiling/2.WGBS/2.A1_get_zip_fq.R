### Jinliang Yang
### Updates: June 24th, 2016

library("farmeR")

### prepare genome
cmd1 <- "module load bismark/0.14.3"
cmd2 <- "module load bowtie2/2.2.5"
cmd3 <- "bismark_genome_preparation --bowtie2 $HOME/Documents/Github/methylation/largedata/pgenome"

set_farm_job(slurmsh = "slurm-script/pgenome.sh",
             shcode = c(cmd1, cmd2, cmd3), wd = NULL, jobid = "pgenome",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "8"))


### fastq file
### extract and gz fastq files
# for i in *.lz4; do lz4 -d $i; done
# for i in *.fastq; do gzip --fast $i; done
# for i in *.lz4; do rm $i; done

### prepare fastq files
cmd0 <- "cd largedata/fastq"
cmd1 <- "for i in *.lz4; do"
cmd2 <- " lz4 -d $i > $i.fastq;"
cmd3 <- " rm $i;"
cmd4 <- " pigz -p 16 $i.fastq"
cmd5 <- "done"
    
set_farm_job(slurmsh = "slurm-script/gzip.sh",
             shcode = c(cmd0, cmd1, cmd2, cmd3, cmd4, cmd5), wd = NULL, jobid = "pigz",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", "16"))


### get it by array
files <- list.files(path="/group/jrigrp4/BS_teo20/fastq", pattern="lz4", full.names=TRUE)
for(i in 1:length(files)){
    shid <- paste0("slurm-script/runmv_", i, ".sh")
    command <- paste("mv", files[i], "/home/jolyang/Documents/Github/methylation/largedata/fastq")
    cat(command, file=shid, sep="\n", append=FALSE)
}

set_array_job(shid="slurm-script/runmv.sh", shcode="sh slurm-script/runmv_$SLURM_ARRAY_TASK_ID.sh",
              arrayjobs= paste0("1-", length(files)), wd=NULL, jobid="mv", 
              email="yangjl0930@gmail.com", run = c(TRUE, "bigmemm", "2"))
