### Jinliang Yang
### July 25th, 2016


library(farmeR)

idtab <- read.csv("data/teo20_ids.csv")
idtab$idchar <- gsub("-", "", idtab$idchar)
idtab$pid <- paste0("JR", idtab$plate)

#####################################################################################################
# bcftools view -r 1:1-1000 JRI20_filtered_snps_annot.bcf.gz --samples JRIAL2A
### checking results
"samtools faidx JRIAL2A.fa 1:1-1000"
# Note: heter=> change to alt, multi=>change to dominant alt, missing=> not change

########## preparing genome
for(i in 1:nrow(idtab)){
    shid <- paste0("slurm-script/run_pg", i, ".sh")
    cmd1 <- paste0("bismark_genome_preparation --bowtie2 largedata/wgbs_pgen/", idtab$idchar[i])
    cat(cmd1, file=shid, sep="\n", append=FALSE)
}    

cmd <- c("module load bismark/0.14.3", "module load bowtie2/2.2.5", "sh slurm-script/run_pg$SLURM_ARRAY_TASK_ID.sh")
set_array_job(shid="slurm-script/run_pg.sh", shcode=cmd,
              arrayjobs="1-20", wd=NULL, jobid="pgjob", email="yangjl0930@gmail.com",
              run = c(TRUE, "bigmemm", "4"))

