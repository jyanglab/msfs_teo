
run_Rcodes <- function(
    inputdf, outdir, cmdno=100,
    rcodes = "lib/C_format.R",
    arrayshid = "slurm-script/run_bcf_query_array.sh",
    email=NULL, runinfo = c(FALSE, "bigmemh", 1)
){
    
    runinfo <- get_runinfo(runinfo)
    #### create dir if not exist
    dir.create("slurm-script", showWarnings = FALSE)
    dir.create(outdir, showWarnings = FALSE)
    
    tot <- ceiling(nrow(inputdf)/cmdno)
    for(j in 1:tot){
        
        shid <- paste0(outdir, "/run_rcode_", j, ".sh")
        
        ##chr:start-end
        #sh1 <- paste("cd", outdir)
        sh <- paste0('R --no-save "--args "', j, ' < ', rcodes)
        
        cat(paste("### run Rcode", Sys.time(), sep=" "),
            sh,
            file=shid, sep="\n", append=FALSE)
    }
    
    shcode <- paste0("sh ", outdir, "/run_rcode_$SLURM_ARRAY_TASK_ID.sh")
    set_array_job(shid=arrayshid, shcode=shcode, arrayjobs=paste("1", tot, sep="-"),
                  wd=NULL, jobid="rcode", email=email, runinfo=runinfo)
}

