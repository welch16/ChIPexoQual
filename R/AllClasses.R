##' @importFrom biovizBase flatGrl
NULL


##' @rdname ExoData
##' @export
setClass("ExoData",
         contains = "GRanges",
         representation = representation(
             file = "character",
             cover = "RleList",
             nreads = "numeric",
             reads = "GRanges",
             param_dist = "list"
         ),
         prototype = prototype(
             file = "",
             cover = RleList(),
             nreads = 0L,
             reads = GRanges(),
             param_dist = list()
         ))

setValidity("ExoData",
            function(object){
                all(names(mcols(object)) == c("f","r","fpos","rpos","d","u",
                                   "ARC","URC","FSR","M","A","label"))
            }
        )

##' ExoData object and constructors
##'
##' \code{ExoData} is a subclass of \code{GenomicRanges}, used to asses the 
##' quality of ChIP-exo/nexus sample.
##'
##' @param file a character value with location of the bam file with the aligned
##' reads.
##' @param height a numeric value indicating the value used to slice the coverage
##' of the experiment into a set of regions.
##' @param mc.cores a numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##' @param nregions a numeric value indicating the number of regions sampled to 
##' estimate the quality parameter distributions. The default value is 1e3.
##' @param ntimes a numeric value indicating the number of times that regions are 
##' sampled to estimate the quality parameter distributions. The default value
##' is 1e3.
##' @param save_reads a logical value to indicate if the reads are stored in the
##' \code{ExoData} object. The default value is \code{FALSE}.
##' @return \code{ExoData} returns a \code{ExoData} object which contains the
##' aggregated coverage of the experiment, the set of islands and a collection
##' of summary statistics used to asses the quality of a ChIP-exo/nexus sample.
##' @aliases ExoData ExoData-class
##'
##' @docType class
##'
##' @examples
##'
##' a = 1
##'
##' @rdname ExoData
##' @export
ExoData = function(file , height = 1 ,mc.cores = getOption("mc.cores",2L),
                   save_reads = FALSE,nregions = 1e3,ntimes = 1e3)
{
    
    stopifnot(is.character(file),file.exists(file))
    stopifnot(is.numeric(height),height >= 1)
    stopifnot(is.logical(save_reads))
    
    param_fwd = ScanBamParam(flag = scanBamFlag(isMinusStrand = FALSE))
    param_bwd = ScanBamParam(flag = scanBamFlag(isMinusStrand = TRUE))    
    
    freads = as(readGAlignments(file,param = param_fwd),"GRanges")
    breads = as(readGAlignments(file,param = param_bwd),"GRanges")
    
    cover = coverage(freads) + coverage(breads)
    
    rlist = slice(cover,lower = height,rangesOnly = TRUE)

    freads = split(freads,as.character(seqnames(freads)))
    breads = split(breads,as.character(seqnames(breads)))
    
    freads = GRangesList(freads[names(rlist)])
    breads = GRangesList(breads[names(rlist)])

    stats = mcmapply(calculate_summary,rlist,freads,breads,
                 mc.cores = mc.cores , SIMPLIFY = FALSE)
    regions = as(rlist,"GRanges")
    
    mcols(regions) = do.call(rbind,stats)
    nreads = sum(vapply(freads,length,1)) + sum(vapply(breads, length, 1))
        
    if(save_reads){
        freads = biovizBase::flatGrl(freads)
        mcols(freads) = NULL
        breads = biovizBase::flatGrl(breads)
        mcols(breads) = NULL
        reads = c(freads,breads)
    }else{
        reads = GRanges()
    }
    
    
    DT = data.table(as.data.frame(mcols(regions)))
    DT = DT[,list(d,u)]
    DT = DT[,w := width(regions)]

    param_list = lapply(seq_len(ntimes),.calculate_param_dist,DT,nregions)
                          
    param = rbindlist(param_list)

    rm(param_list,DT)
    param_dist = list("beta1" = param[term == "u",(estimate)] ,
                      "beta2" = param[term == "w",(estimate)])
    
    new("ExoData",regions,file = file,cover = cover,
        nreads = nreads,reads = reads,
        param_dist = param_dist)
}

# ##' ChIP-exo experiment class description.
# ##' 
# ##' Contains the reads of a ChIP exo experiment.
# ##' 
# ##' @slot file Character vector wih the name of the bam file of a ChIP-exo experiment.
# ##' @slot reads A \code{GRanges} object with the reads of the experiment.
# ##' @slot depth A numeric value with the experiment's depth.
# ##' @slot height Tunning parameter to partition the experiment into regions.
# ##' @slot summary_stats A \code{data.table} with a collection of summary statistics of the ChIP-exo experiment.
# ##' 
# ##' @name ChIPexo_experiment-class
# ##' @rdname ChIPexo_experiment-class
# ##' @exportClass ChIPexo_experiment
# setClass("ChIPexo_experiment",
#     representation = representation(
#       file = "character",
#       reads = "GRanges",
#       depth = "numeric",
#       height = "numeric",
#       summary_stats = "data.table"
#     ),
#     prototype = prototype(
#       file = "",
#       reads = GRanges(),
#       depth = 0,
#       height = 1,
#       summary_stats = data.table()
#     )
# )
