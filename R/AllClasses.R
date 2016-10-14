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
             paramDist = "list"
         ),
         prototype = prototype(
             cover = RleList(),
             reads = GRanges(),
             paramDist = list()
         ))

setValidity("ExoData",
            function(object){
                all(names(mcols(object)) == c("fwdReads","revReads",
                                              "fwdPos","revPos",
                                              "depth","uniquePos",
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
##' @param save.reads a logical value to indicate if the reads are stored in the
##' \code{ExoData} object. The default value is \code{FALSE}.
##' @param verbose a logical value indicating if the user want to receive progress
##' details. The default value is FALSE.
##' @return It returns an \code{ExoData} object with the regions obtained after
##' partitioning the genome and the summary statistics for each region. If the
##' \code{save.reads} parameter is \code{TRUE} then it contains a \code{GRanges}
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
ExoData = function(file = NULL, reads = NULL , height = 1 ,
                   mc.cores = getOption("mc.cores",2L),
                   save.reads = FALSE,nregions = 1e3,ntimes = 1e2,
                   verbose = TRUE)
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
    stopifnot(is.logical(save.reads))
    
    if(verbose){
        if(is.null(file))message("Using 'reads' argument")
        else message("Using 'file' argument")
    }
    if(!is.null(file)){
        if(verbose) message("Creating ExoData object using aligned reads in ",file)
    }
    if(verbose) message("Keeping reads in object: ",ifelse(save.reads,"Yes","No"))
    if(verbose) message("Loading experiment reads")
    
    if(!is.null(file)){
        paramFwd = ScanBamParam(flag = scanBamFlag(isMinusStrand = FALSE))
        paramRev = ScanBamParam(flag = scanBamFlag(isMinusStrand = TRUE))    
        
        fwdReads = readGAlignments(file,param = paramFwd)
        revReads = readGAlignments(file,param = paramRev)
    }else{
        fwdReads = subset(reads,as.character(strand(reads)) == "+")
        revReads = subset(reads,as.character(strand(reads)) == "-")
        file = ""
    }
    
    if(class(fwdReads) == "GAlignments")fwdReads = as(fwdReads,"GRanges")
    if(class(revReads) == "GAlignments")revReads = as(revReads,"GRanges")
    
    cover = coverage(fwdReads) + coverage(revReads)
    
    rlist = slice(cover,lower = height,rangesOnly = TRUE)
    
    if(any(vapply(rlist,length,0L) ==  0)){
        chr = names(which(vapply(rlist,length,0L) > 0))
        rlist = rlist[chr]
    }else{
        chr = names(rlist)
    }
    
    fwdReads = split(fwdReads,as.character(seqnames(fwdReads)))
    revReads = split(revReads,as.character(seqnames(revReads)))
    
    fwdReads = fwdReads[chr]
    revReads = revReads[chr]
    
    if(verbose) message("Calculating summary statistics")
    
    if(Sys.info()[["sysname"]] == "Windows"){
        snow = SnowParam(workers = mc.cores, type = "SOCK")
        stats = bpmapply(calculateSummary,rlist,fwdReads,revReads,
                         BPPARAM = snow,SIMPLIFY = FALSE)       
    }else{
        stats = bpmapply(calculateSummary,rlist,fwdReads,revReads,
            BPPARAM = MulticoreParam(workers = mc.cores),
            SIMPLIFY = FALSE)       
    }

    regions = as(rlist,"GRanges")
    mcols(regions) = do.call(rbind,stats)
    nreads = sum(vapply(fwdReads,length,1)) + sum(vapply(revReads, length, 1))
    
    if(save.reads){
        fwdReads = biovizBase::flatGrl(fwdReads)
        mcols(fwdReads) = NULL
        revReads = biovizBase::flatGrl(revReads)
        mcols(revReads) = NULL
        reads = c(fwdReads,revReads)
    }else{
        reads = GRanges()
    }
    depth <- uniquePos <- width <-  NULL
    
    if(verbose) message("Calculating quality scores distribution")
    DT = data.table(as.data.frame(mcols(regions)))
    DT = DT[,list(depth,uniquePos)]
    DT = DT[,width := width(regions)]

    paramList = lapply(seq_len(ntimes),calculateParamDist,DT,nregions)
                          
    param = rbindlist(paramList)

    estimate <- term <- NULL
    
    rm(paramList,DT)
    paramDist = list("beta1" = param[term == "uniquePos",(estimate)] ,
                      "beta2" = param[term == "width",(estimate)])
    
    metadata(regions) = list("file"=file,"nreads"=nreads,
                             "ntimes"=ntimes,"nregions"=nregions)
    if(verbose) message("Done!")
    new("ExoData",regions,cover = cover,
        reads = reads,
        paramDist = paramDist)
}
