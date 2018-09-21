#' \code{Run Pseudo-ref job on farm}
#'
#' GATK Best Practices: recommended workflows for variant discovery analysis.
#'
#' see more detail about GATK:
#' \url{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}
#'
#' idxing:
#' bwa index Zea_mays.AGPv2.14.dna.toplevel.fa
#'
#' module load java/1.8
#' module load bwa/0.7.9a
#'
#' local programs:
#' bwa Version: 0.7.5a-r405
#' picard-tools-2.1.1
#' GenomeAnalysisTK-3.5/
#'
#' @param inputdf An input data.frame for fastq files. Must contains fq1, fq2, out (and/or bam).
#' If inputdf contained bam, bwa alignment will be escaped.
#' Additional columns: group (group id), sample (sample id), PL (platform, i.e. illumina),
#' LB (library id), PU (unit, i.e. unit1). These strings (or info) will pass to BWA mem through -R.
#'
#' @param ref.fa The full path of genome with bwa indexed reference fasta file.
#' @param gatkpwd The absolute path of GenomeAnalysisTK.jar.
#' @param email Your email address that farm will email to once the jobs were done/failed.
#' @param runinfo Parameters specify the array job partition information.
#' A vector of c(FALSE, "bigmemh", "1"): 1) run or not, default=FALSE
#' 2) -p partition name, default=bigmemh and 3) --cpus, default=1.
#' It will pass to \code{set_array_job}.
#'
#' @return return a batch of shell scripts.
#'
#' @examples
#' inputdf <- data.frame(input.vcf="test_1.vcf", out.fa="mysample.fa")
#'
#'
#' @export
run_pseudoref <- function(inputdf,
                     ref.fa="~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta",
                     gatkpwd="$HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar",
                     email=NULL, runinfo = c(FALSE, "bigmemh", 1)){
    
    ### determine memory based on partition
    runinfo <- get_runinfo(runinfo)
    
    #### create dir if not exist
    dir.create("slurm-script", showWarnings = FALSE)
    for(i in 1:nrow(inputdf)){
        
        shid <- paste0("slurm-script/run_pseudoref_", i, ".sh")
        ### header of the shell code
        cat("### GATK pipeline created by farmeR",
            paste("###", format(Sys.time(), "%a %b %d %X %Y")),
            paste(""),
            paste0("java -Xmx", floor(as.numeric(runinfo[4])/1024), "g ",
                   "-jar ", gatkpwd, " \\"),

            "-T FastaAlternateReferenceMaker \\",
            paste0("-R ", ref.fa, " \\"),
            paste0("-o ", inputdf$out.fa[i], " \\"),
            #"-L input.intervals \\",
            paste0("-V ", inputdf$input.vcf[i], " \\"),
            
            file=shid, sep="\n", append=FALSE)
    }
    
    shcode <- paste("module load java/1.8", "module load bwa/0.7.9a",
                    "sh slurm-script/run_pseudoref_$SLURM_ARRAY_TASK_ID.sh", sep="\n")
    set_array_job(shid="slurm-script/run_pseudoref_array.sh",
                  shcode=shcode, arrayjobs=paste("1", nrow(inputdf), sep="-"),
                  wd=NULL, jobid="pref", email=email, runinfo=runinfo)
    #  sbatch -p bigmemh --mem 32784 --ntasks=4  slurm-script/run_gatk_array.sh
}





