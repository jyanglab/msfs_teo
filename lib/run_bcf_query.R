
run_bcf_query <- function(
    inputdf, outdir="largedata/Dm/input_gene", cmdno=100,
    arrayshid = "slurm-script/run_bcf_query_array.sh",
    email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){
    
    runinfo <- get_runinfo(runinfo)
    #### create dir if not exist
    dir.create("slurm-script", showWarnings = FALSE)
    dir.create(outdir, showWarnings = FALSE)
    
    ### setup shell id
    ####
    tot <- ceiling(nrow(inputdf)/cmdno)
    for(j in 1:tot){
        
        srow <- cmdno*(j-1) + 1
        if(j == tot){
            erow <- nrow(inputdf)
        }else{
            erow <- cmdno*j
        }
        
        shid <- paste0(outdir, "/run_bcf_", j, ".sh")
        
        ##chr:start-end
        sh1 <- paste("cd", outdir)
        sh2 <- paste0("bcftools query -f \'%CHROM\\t%POS\\t%CO[\\t%GT]\\n\'", 
                      " -r ", inputdf$seqname[srow:erow], ":", inputdf$start[srow:erow], "-",
                      inputdf$end[srow:erow],
                      " ~/Documents/Github/methylation/largedata/vcf_files/teo20_methratio.bcf",
                      " -o ", inputdf$attribute[srow:erow])
        sh3 <- paste0('sed -i \"s/\\//\\t/g\" ', inputdf$attribute[srow:erow])
        cat(paste("### run GenSel", Sys.time(), sep=" "),
            c(sh1, sh2, sh3),
            file=shid, sep="\n", append=FALSE)
    }
    
    shcode <- paste0("sh ", outdir, "/run_bcf_$SLURM_ARRAY_TASK_ID.sh")
    set_array_job(shid=arrayshid, shcode=shcode, arrayjobs=paste("1", tot, sep="-"),
                  wd=NULL, jobid="bcf-query", email=email, runinfo=runinfo)
}

