##' @importFrom biovizBase flatGrl
##' @importFrom Rsamtools ScanBamParam
##' @importFrom Rsamtools scanBamFlag
##' @importFrom GenomicAlignments readGAlignments 
##' @importFrom GenomicRanges GRanges
##' @importFrom IRanges slice
##' @importFrom data.table ":=" rbindlist
##' @import  methods
##' @import  GenomeInfoDb
##' @import BiocParallel
NULL


##' @rdname ExoData
##' @export
setClass("ExoData",
         contains = "GRanges",
         representation = representation(
             cover = "RleList",
             reads = "GRanges",
             param_dist = "list"
         ),
         prototype = prototype(
             cover = RleList(),
             reads = GRanges(),
             param_dist = list()
         ))

setValidity("ExoData",
            function(object){
                all(names(mcols(object)) == c("f","r","fpos","rpos","d","u",
                                   "ARC","URC","FSR","M","A"))
            }
        )

##' ExoData object and constructors
##'
##' \code{ExoData} is a subclass of \code{GenomicRanges}, used to asses the 
##' quality of ChIP-exo/nexus sample.
##'
##' @param file a character value with location of the bam file with the aligned
##' reads.
##' @param reads a \code{GAlignments} object with the aligned reads of a ChIP-exo
##' sample. It is meant to be used instead of \code{file}.
##' @param height a numeric value indicating the value used to slice the coverage
##' of the experiment into a set of regions.
##' @param mc.cores a numeric value with the number of cores to use,
##' i.e. at most how many child processes will be run simultaneously.
##' @param nregions a numeric value indicating the number of regions sampled to 
##' estimate the quality parameter distributions. The default value is 1e3.
##' @param ntimes a numeric value indicating the number of times that regions are 
##' sampled to estimate the quality parameter distributions. The default value
##' is 1e2.
##' @param save_reads a logical value to indicate if the reads are stored in the
##' \code{ExoData} object. The default value is \code{FALSE}.
##' @param verbose a logical value indicating if the user want to receive progress
##' details. The default value is FALSE.
##' @return It returns an \code{ExoData} object with the regions obtained after
##' partitioning the genome and the summary statistics for each region. If the
##' \code{save_reads} parameter is \code{TRUE} then it contains a \code{GRanges}
##' object with the reads of the ChIP-exo experiment.
##' @aliases ExoData ExoData-class
##'
##' @docType class
##'
##' @examples
##' 
##' files = list.files(system.file("extdata",package = "ChIPexoQualExample"),
##'     full.names = TRUE)
##' ExoData(files[5],mc.cores = 2L)
##' 
##'
##' @rdname ExoData
##' @export
ExoData = function(file = NULL, reads = NULL , height = 1 ,mc.cores = getOption("mc.cores",2L),
                   save_reads = FALSE,nregions = 1e3,ntimes = 1e2,verbose = FALSE)
{
    
    if(!is.null(file) & !is.null(reads)){
        stop("Both 'file' and 'reads' are available, can't use both.")
    }
    if(is.null(file) & is.null(reads)){
        stop("Both 'file' and 'reads' are NULL")
    }

    if(!is.null(file))stopifnot(is.character( file),file.exists(file))
    
    if(!is.null(reads))stopifnot(class(reads) %in% c("GAlignments","GRanges"))
    
    stopifnot(is.numeric(height),height >= 1)
    stopifnot(is.logical(save_reads))
    
    if(verbose){
        if(is.null(file))message("Using 'reads' argument")
        else message("Using 'file' argument")
    }
    if(!is.null(file)){
        if(verbose) message("Creating ExoData object using aligned reads in ",file)
    }
    if(verbose) message("Keeping reads in object: ",ifelse(save_reads,"Yes","No"))
    if(verbose) message("Loading experiment reads")
    
    if(!is.null(file)){
        param_fwd = ScanBamParam(flag = scanBamFlag(isMinusStrand = FALSE))
        param_bwd = ScanBamParam(flag = scanBamFlag(isMinusStrand = TRUE))    
        
        freads = readGAlignments(file,param = param_fwd)
        breads = readGAlignments(file,param = param_bwd)
    }else{
        freads = subset(reads,as.character(strand(reads)) == "+")
        breads = subset(reads,as.character(strand(reads)) == "-")
        file = ""
    }
    
    if(class(freads) == "GAlignments")freads = as(freads,"GRanges")
    if(class(breads) == "GAlignments")breads = as(breads,"GRanges")
    
    cover = coverage(freads) + coverage(breads)
    
    rlist = slice(cover,lower = height,rangesOnly = TRUE)
    
    if(any(vapply(rlist,length,0L) ==  0)){
        chr = names(which(vapply(rlist,length,0L) > 0))
        rlist = rlist[chr]
    }else{
        chr = names(rlist)
    }
    

    freads = split(freads,as.character(seqnames(freads)))
    breads = split(breads,as.character(seqnames(breads)))
    
    freads = freads[chr]
    breads = breads[chr]
    
#     freads = GRangesList(freads[chr])
#     breads = GRangesList(breads[chr])

    if(verbose) message("Calculating summary statistics")
    
    if(Sys.info()[["sysname"]] == "Windows"){
        snow = SnowParam(workers = mc.cores, type = "SOCK")
        stats = bpmapply(.calculate_summary,rlist,freads,breads,
                         BPPARAM = snow,SIMPLIFY = FALSE)       
    }else{
        stats = bpmapply(.calculate_summary,rlist,freads,breads,
            BPPARAM = MulticoreParam(workers = mc.cores),
            SIMPLIFY = FALSE)       
    }
    
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
    d = NULL; u = NULL; w = NULL
    
    if(verbose) message("Calculating quality scores distribution")
    DT = data.table(as.data.frame(mcols(regions)))
    DT = DT[,list(d,u)]
    DT = DT[,w := width(regions)]

    param_list = lapply(seq_len(ntimes),.calculate_param_dist,DT,nregions)
                          
    param = rbindlist(param_list)

    estimate = NULL; term = NULL
    
    rm(param_list,DT)
    param_dist = list("beta1" = param[term == "u",(estimate)] ,
                      "beta2" = param[term == "w",(estimate)])
    
    metadata(regions) = list("file"=file,"nreads"=nreads,
                             "ntimes"=ntimes,"nregions"=nregions)
    if(verbose) message("Done!")
    new("ExoData",regions,cover = cover,
        reads = reads,
        param_dist = param_dist)
}
