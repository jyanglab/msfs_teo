
run_dm_test <- function(
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
        
        shid <- paste0(outdir, "/run_dm_", j, ".sh")
        
        ##chr:start-end
        sh1 <- paste("cd", outdir)
        
        sh2 <- paste0("Dm_methylation.pl -input ", inputdf$locus_name[srow:erow],
                      " -output ", inputdf$out[srow:erow], " -length ", inputdf$len[srow:erow],
                      " -alpha ", inputdf$alpha_value[srow:erow])
        
        cat(paste("### run Dm_methylation.pl", Sys.time(), sep=" "),
            c(sh1, sh2),
            file=shid, sep="\n", append=FALSE)
    }
    
    shcode <- paste0("sh ", outdir, "/run_dm_$SLURM_ARRAY_TASK_ID.sh")
    set_array_job(shid=arrayshid, shcode=shcode, arrayjobs=paste("1", tot, sep="-"),
                  wd=NULL, jobid="dm-test", email=email, runinfo=runinfo)
}

