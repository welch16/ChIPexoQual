##' Samples \code{sample.reads} from the ChIP-exo experiment and creates a list
##' of \code{ExoData} objects
##'
##'
##' @param file a character value with location of the bam file with the aligned
##' reads.
##' @param reads a \code{GAlignments} object with the aligned reads of a ChIP-exo
##' sample. It is meant to be used instead of \code{file}.
##' @param sample.depth a numeric vector with the number of reads to be sampled.
##' @param height a numeric value indicating the value used to slice the coverage
##' of the experiment into a set of regions.
##' @param mc.cores a numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##' @param nregions a numeric value indicating the number of regions sampled to 
##' estimate the quality parameter distributions. The default value is 1e3.
##' @param ntimes a numeric value indicating the number of times that regions are 
##' sampled to estimate the quality parameter distributions. The default value
##' is 1e2.
##' @param save.reads a logical value to indicate if the reads are stored in the
##' \code{ExoData} object. The default value is \code{FALSE}.
##' @param verbose a logical value indicating if the user want to receive progress
##' details. The default value is FALSE.
##' @return It returns an \code{ExoData} object with the regions obtained after
##' partitioning the genome and the summary statistics for each region. If the
##' \code{save.reads} parameter is \code{TRUE} then it contains a \code{GRanges}
##' object with the reads of the ChIP-exo experiment.
##'
##' @examples
##' 
##' files = list.files(system.file("extdata",package = "ChIPexoQualExample"),
##'     full.names = TRUE)
##' sample.depth = seq(1e5,2e5,5e4)
##' ExoDataSubsampling(file = files[5],sample.depth = sample.depth)
##'
##' @rdname ExoDataSubsampling
##' @export
ExoDataSubsampling = function(file = NULL, reads = NULL,
                              sample.depth = NULL,
                              height = 1, 
                              nregions = 1000,
                              ntimes = 1000,
                              verbose = TRUE,
                              save.reads = FALSE,
                              mc.cores = getOption("mc.cores",2L))
{
    ExoList = NULL
    if(!is.null(file) & !is.null(reads)){
        stop("Both 'file' and 'reads' are available, can't use both.")
    }else if(is.null(file) & is.null(reads)  ){
        stop("Both 'file' and 'reads' are NULL")
    }else{
        if(is.null(reads)){
            stopifnot(file.exists(file))
            stopifnot(!is.null(sample.depth))
            stopifnot(is.numeric(sample.depth))
            reads = readGAlignments(file,param = NULL)
            
        }
        depth = length(reads)
        if(any(sample.depth > depth)){
            stop("There are values in 'sample.depth' greater than the
                experiments depth, please consider smaller values 
                in 'sample.depth'")
        }
        
        auxOpt = getOption("scipen")
        options(scipen = 999)
        ExoList = lapply(sample.depth,
                 function(n,reads,height,nregions,ntimes,
                    verbose,save.reads,mc.cores ){
                    if(verbose){
                        message("Creating 'ExoData' for ",
                                prettyNum(n,big.mark = ","),
                                " reads")
                    }
                    sampleReads = sample(reads,n)
                    ExoData(reads = sampleReads,
                        height = height,nregions = nregions,
                        ntimes = ntimes,verbose = verbose,
                        save.reads = save.reads,
                        mc.cores = mc.cores)
                    },reads,height,nregions,ntimes,
            verbose,save.reads,mc.cores)
        names(ExoList) = paste(prettyNum(as.character(sample.depth),big.mark = ","),
                    "reads")
        options(scipen = auxOpt)
    }
    ExoList
}