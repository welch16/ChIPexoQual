##' @importFrom stats quantile
##' @importFrom S4Vectors subset mcols DataFrame "metadata<-" metadata
##' @importFrom scales alpha trans_format math_format
##' @importFrom viridis viridis
##' @importFrom RColorBrewer brewer.pal
##' @importFrom dplyr desc arrange
##' @importFrom hexbin hexbin
##' @import ggplot2
##' @import GenomicRanges
NULL

##' Calculates the FSR distribution of all regions with depth > value
##'  
##' Calculates the FSR distribution of all regions with depth > value
##'  
##' @param value a numeric value with the depth lower bound.
##' @param DF a \code{DataFrame} with depth and FSR.
##' @param quantiles a numeric vector with the quantiles used to summarize
##' the FSR distribution.
##' 
##' @return a \code{DataFrame} with three columns depth, quantiles and FSR.
##' 
##' @rdname .filterQuantiles
##' @name .filterQuantiles
##' 
.filterQuantiles <- function(value,DF,quantiles)
{
    depth <- NULL
    dist <- S4Vectors::subset(DF,depth > value)
    quant <- quantile(dist$FSR,probs = quantiles)
    names(quant) <- NULL
    DataFrame(depth = value ,quantiles , FSR = quant)
}    

##' Calculates the region composition probabilities for all 
##' regions with depth > value
##'  
##' Calculates the region composition probabilities for all 
##' regions with depth > value
##'  
##' @param value a numeric value with the depth lower bound.
##' @param DF a \code{DataFrame} with depth and FSR.
##' 
##' @return a \code{DataFrame} with three columns d, label and prob.
##' 
##' @rdname .filterRegionComp
##' @name .filterRegionComp
##' 
.filterRegionComp <- function(value,DF)
{
    depth <- NULL
    DF <- S4Vectors::subset(DF,depth > value)
    lab <- c("both","fwd","bwd")
    tabl <- table(fwd = factor(DF$fwdReads > 0,levels = c(TRUE,FALSE)) ,
                 bwd = factor(DF$revReads > 0,levels = c(TRUE,FALSE)))
    counts <- c(tabl[1,1],tabl[1,2],tabl[2,1])
    props <- counts / sum(counts)
    DataFrame(depth = value,lab,prob = props )
}

.generateNames <- function(names.input,nms, length){

    if(is.null(names.input)){
        if(is.null(nms)){
            nms <- paste0("Sample: ",seq_len(length))
        }
    }else{
        nms <- names.input
    }
    nms
}

.nameJoin <- function(my_list,nms){
    
    DF <- do.call(rbind,my_list)
    nms <- mapply(rep,nms,each = vapply(my_list,nrow,1L),SIMPLIFY = FALSE)
    nms <- do.call(c,nms)
    DF$sample <- nms
    data.table(as.data.frame(DF))
    
}

##' .MADataFrame
##' 
##' Returns a \code{DataFrame} with the info. to generate a MA
##' plot to analyze the strand imbalance in ChIP-exo data.
##' 
##' @param object A \code{ExoData} object.
##' 
##' @return A \code{DataFrame} with two columns: M and A.
##' @rdname .MADataFrame
##' @name .MADataFrame
##' @examples 
##' data(exoExample)
##' .MADataFrame(exoExample)
.MADataFrame <- function(object)
{
    A <- NULL
    DF <- mcols(object)[,c("M","A")]
    DF <- subset(DF,!is.infinite(A))
    DF
    
}

##' MAplot 
##' 
##' \code{MAplot} returns a \code{ggplot} object with the MA plot to 
##' analyze the strand imbalance in ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names.input a character vector with the names to use in the
##' plot. If it is empty \code{MAplot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' 
##' @return A \code{ggplot2} object with the MA plot.
##' @rdname MAplot
##' @name MAplot
##' @export
##' @examples 
##' data(exoExample)
##' MAplot(exoExample)
MAplot <- function(...,names.input = NULL)
{
    M <- A <- .x <- NULL
    args <- unlist(list(...))
    if(!is.null(names.input)){
        stopifnot(length(names.input) == length(args))
    }
    MA_list <- lapply(args,.MADataFrame)
    nsamples <- length(MA_list)
    nms <- .generateNames(names.input,names(MA_list), 
                          nsamples)
    names(MA_list) <- NULL
    MA_DF <- .nameJoin(MA_list,nms)
    r <- viridis(1e3, option = "D")
    
    theme_set(theme_bw())
    
    p <- ggplot(MA_DF,aes(M,A))+stat_binhex(bins = 70)+
        scale_fill_gradientn(colours = r,trans = 'log10',
            labels=trans_format('log10',math_format(10^.x)) )+
        theme(legend.position = "top")
    if(nsamples <= 3){
        p <- p + facet_wrap( ~ sample,nrow = 1)
    }else{
        p <- p + facet_wrap( ~ sample,ncol = 4)
    }
    p
}

##' .ARCvURCDataFrame
##' 
##' \code{.ARCvURCDataFrame} returns a \code{DataFrame} object with the ARC and
##' URC columns used to plot enrichment and library complexity in ChIP-exo data.
##' 
##' @param object a \code{ExoData} object.
##' @param both.strand A logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{DataFrame} with the ARC and URC columns of the \code{ExoData}
##' object.
##' @rdname .ARCvURCDataFrame
##' @name .ARCvURCDataFrame
## @export
##' @examples 
##' data(exoExample)
##' .ARCvURCDataFrame(exoExample,both.strand = FALSE)
##' .ARCvURCDataFrame(exoExample,both.strand = TRUE)
.ARCvURCDataFrame = function(object,both.strand)
{
    fwdReads <- revReads <-  NULL
    
    DF <- mcols(object)[,c("ARC","URC","fwdReads","revReads")]
    if(both.strand){
        DF <- subset(DF,fwdReads > 0 & revReads > 0)
    }
    DF[,c("ARC","URC")]
}


##' ARCvURCplot 
##' 
##' \code{ARCvURCplot} returns a \code{ggplot} object with the ARC vs
##' URC plot to analyze enrichment and library complexity in ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names.input a character vector with the names to use in the
##' plot. If it is empty \code{ARCvURCplot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' @param both.strand A logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{ggplot2} object with the ARC vs URC plot.
##' @rdname ARCvURCplot
##' @name ARCvURCplot
##' @export
##' @examples 
##' data(exoExample)
##' ARCvURCplot(exoExample)
ARCvURCplot <- function(...,names.input = NULL,both.strand = FALSE)
{
    
    ARC <- URC <- .x <- NULL
    
    args <- unlist(list(...))
    if(!is.null(names.input)){
        stopifnot(length(names.input) == length(args))
    }
    
    ARCvURCList <- lapply(args,.ARCvURCDataFrame,both.strand)
    nsamples <- length(ARCvURCList)
    nms <- .generateNames(names.input,names(ARCvURCList), 
                          nsamples)
    names(ARCvURCList) <- NULL
    .ARCvURCDataFrame <- .nameJoin(ARCvURCList,nms)
    r <- viridis(1e3, option = "D")
    
    theme_set(theme_bw())
    
    p <- ggplot(.ARCvURCDataFrame,aes(ARC,URC))+stat_binhex(bins = 50)+
        scale_fill_gradientn(colours = r,trans = 'log10',
                             labels=trans_format('log10',math_format(10^.x)) )+
        theme(legend.position = "top")
    if(nsamples <= 3){
        p <- p + facet_wrap( ~ sample,nrow = 1)
    }else{
        p <- p + facet_wrap( ~ sample,ncol = 4)
    }
    p <- p + xlim(0,3)+ylim(0,1)
    p
}

##' .FSRDistDataFrame
##' 
##' \code{.FSRDistDataFrame} returns a \code{DataFrame} object with the Forward 
##' Strand Ratio distribution of the \code{ExoData} object to analyze strand 
##' imbalance in ChIP-exo data.
##' 
##' @param object a \code{ExoData} object.
##' @param quantiles a numeric vector with the quantiles used to estimate the
##' FSR distribution at a given depth. 
##' @param depth.values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. 
##' @param both.strand a logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. 
##' 
##' @return A \code{DataFrame} object with the FSR distribution of the 
##' \code{ExoData} object.
##' @rdname .FSRDistDataFrame
##' @name .FSRDistDataFrame
## @export
##' @examples 
##' data(exoExample)
##' .FSRDistDataFrame(exoExample,quantiles = c(.25,.5,.75),
##'   depth.values = seq_len(25),both.strand = FALSE)
.FSRDistDataFrame <- function(object,quantiles,depth.values,
                                both.strand)
{
    fwdReads <- revReads <- NULL
    
    baseDF <- mcols(object)[,c("fwdReads","revReads","depth","FSR")]
    
    if(both.strand){
        baseDF <- subset(baseDF,fwdReads > 0 & revReads > 0)
    }
    
    DFlist <- lapply(depth.values,.filterQuantiles,
                    baseDF,quantiles)
    do.call(rbind,DFlist)
}


##' FSRDistplot 
##' 
##' \code{FSRDistplot} returns a \code{ggplot} object with the Forward 
##' Strand Ratio distribution plot to analyze strand imbalance in 
##' ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names.input a character vector with the names to use in the
##' plot. If it is empty \code{FSRDistplot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' @param quantiles a numeric vector with the quantiles used to estimate the
##' FSR distribution at a given depth. The default value is 
##' \code{c(0,.25,.5,.75,1))}
##' @param depth.values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' @param both.strand a logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{ggplot2} object with the FSR distribution plot.
##' @rdname FSRDistplot
##' @name FSRDistplot
##' @export
##' @examples 
##' data(exoExample)
##' FSRDistplot(exoExample)
FSRDistplot <- function(...,names.input = NULL,
                         quantiles = c(0,.25,.5,.75,1),
                         depth.values = seq_len(30),
                         both.strand = FALSE)
{
    depth <- FSR <- NULL
    
    args <- unlist(list(...))
    if(!is.null(names.input)){
        stopifnot(length(names.input) == length(args))
    }
    FSRList <- lapply(args,.FSRDistDataFrame,quantiles,
                      depth.values,
                      both.strand)
    nsamples <- length(FSRList)
    nms <- .generateNames(names.input,names(FSRList), 
                          nsamples)
    names(FSRList) <- NULL
    FSRDataFrame <- .nameJoin(FSRList,nms)
    
    theme_set(theme_bw())
    
    p <- ggplot(FSRDataFrame,aes(depth,FSR,colour = as.factor(quantiles)))+
        geom_line(size = 1)+
        theme(legend.position = "top")+facet_grid(sample ~ .)+
        scale_color_brewer(palette = "Dark2",name = "")+
        xlab("Minimum number of reads")+ylab("Forward Strand Ratio \n (FSR)")+
        ylim(0,1)
    p
}

##' .regionCompDataFrame
##' 
##' \code{.regionCompDataFrame} returns a \code{DataFrame} with the info. 
##' necessary to generate the Region Composition plot to analyze strand
##' imbalance in ChIP-exo data.
##' 
##' @param object a \code{ExoData} object.
##' @param depth.values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' @return A \code{DataFrame} with three columns: depth, label and propotion
##' of regions at a given depth level with the respective label.
##' @rdname .regionCompDataFrame
##' @name .regionCompDataFrame
## @export
##' @examples 
##' data(exoExample)
##' .regionCompDataFrame(exoExample,seq_len(10))
.regionCompDataFrame <- function(object,depth.values)
{
    baseDF <- mcols(object)[,c("fwdReads","revReads","depth")]
    DataFrameList <- lapply(depth.values,.filterRegionComp,
                           baseDF)
    do.call(rbind,DataFrameList)
}

##' regionCompplot 
##' 
##' \code{regionCompplot} returns a \code{ggplot} object with the 
##' Region Composition plot to analyze strand imbalance in ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names.input a character vector with the names to use in the
##' plot. If it is empty \code{regionCompplot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' @param depth.values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' 
##' @return A \code{ggplot2} object with the Region Composition plot.
##' @rdname regionCompplot
##' @name regionCompplot
##' @export
##' @examples 
##' data(exoExample)
##' regionCompplot(exoExample)
regionCompplot <- function(...,names.input = NULL,
                            depth.values = seq_len(15))
{
    lab <- depth <- prob <- NULL
    
    args <- unlist(list(...))
    if(!is.null(names.input)){
        stopifnot(length(names.input) == length(args))
    }
    regionList <- lapply(args,.regionCompDataFrame,depth.values)
    nsamples <- length(regionList)
    nms <- .generateNames(names.input,names(regionList), 
                          nsamples)
    names(regionList) <- NULL
    regionDataFrame <- .nameJoin(regionList,nms)
    r <- brewer.pal(name = "Set1",3)
    regionDataFrame <- regionDataFrame[,lab := 
        factor(lab,levels = c("both","fwd","bwd"))]
    
    theme_set(theme_bw())
    
    p <- ggplot(arrange(regionDataFrame,lab, desc(prob)),
                  aes(depth,prob,fill = lab))+
        geom_bar(stat = "identity")+
        theme(legend.position = "top")+
        facet_grid(sample ~ .)+
        scale_fill_brewer(palette = "Pastel1",name = "Strand composition")+
        xlab("Minimum number of reads")+ylab("Proportion of regions")
    p
    
}

##' paramDistBoxplot 
##' 
##' \code{paramDistBoxplot} returns a \code{ggplot} object with a 
##' boxplot comparing the \code{ntimes} estimations of the chosen 
##' parameter.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param which.param a character value with either \code{"beta1"} or 
##' \code{"beta2"} that determines which paramters in the model 
##' depth_i ~ uniquePos_i + width_i to plot. The default value is 
##' \code{"beta1"}.
##' @param names.input a character vector with the names to use in the
##' plot. If it is empty \code{paramDistBoxplot} is going to create the 
##' names as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' 
##' @return A \code{ggplot2} object with the boxplot of the chosen 
##' parameter
##' @rdname paramDistBoxplot
##' @name paramDistBoxplot
##' @export
##' @examples 
##' data(exoExample)
##' paramDistBoxplot(exoExample)
paramDistBoxplot <- function(...,names.input = NULL,
                            which.param = "beta1")
{

    lab <- depth <- prob <- NULL
    
    stopifnot(which.param %in% c("beta1","beta2"))
    args <- unlist(list(...))
    if(!is.null(names.input)){
        stopifnot(length(names.input) == length(args))
    }
    paramList <- lapply(args,paramDist)
    nsamples <- length(paramList)
    nms <- .generateNames(names.input,names(paramList), 
                          nsamples)
    names(paramList) <- NULL
    paramDataFrame <- .nameJoin(paramList,nms)
    paramDataFrame$beta2 <- - paramDataFrame$beta2
    p <- ggplot(paramDataFrame,
                aes_string(x = "sample",
                           y = which.param))+
        geom_boxplot()+
        theme_bw()+
        theme(legend.position = "top",
              axis.title.y = element_text(angle = 0),
              axis.text.x = element_text(angle = 30,hjust = 1))
    if(which.param == "beta1"){
        p <- p + ylab(expression(beta[1]))+
            geom_abline(slope = 0,intercept = 10,linetype = 2)
    }else{
        p <- p + ylab(expression(beta[2]))+
            geom_abline(slope = 0,intercept = 0,linetype = 2)
    }
    p+xlab("Sample")
}

