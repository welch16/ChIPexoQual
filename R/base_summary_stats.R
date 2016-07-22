##' @importFrom S4Vectors DataFrame
##' @importFrom broom tidy
##' @importFrom IRanges start end
##' @importFrom stats lm
NULL

##' Calculate the summary statistics for each region
##' 
##' @param region a IRanges object with the regions of a certain chromosome
##' @param freads a GRanges object with the fwd. reads of the ChIP-exo experiment
##' @param breads a GRanges object with the bwd. reads of the ChIP-exo experiment
##' 
##' @return stats a \code{DataFrame} object with the summary statistics necessary
##' for the \code{ExoData}
##' 
##' @export
##' 
##' @rdname calculate_summary
##' @name calculate_summary
##' 
##' @examples 
##' a = 1
calculate_summary <- function(region,freads,breads)
{
    ## fix formats and stuff
    width(freads) = 1; width(breads) = 1
    
    freads = ranges(freads); breads = ranges(breads)
    
    fwd = countOverlaps(region,freads)
    bwd = countOverlaps(region,breads)
    
    fpos = IRanges(unique(start(freads)),width = 1)
    bpos = IRanges(unique(end(breads)),width = 1)
    fpos = countOverlaps(region,fpos)
    bpos = countOverlaps(region,bpos)
    w = width(region)    
    
    d = fwd + bwd
    u = fpos + bpos
    
    arc = d /w 
    urc = u / d 
    fsr = fwd / d
    
    M = log2(fwd) + log2(bwd) - 2 * log2(w)
    A = log2(fwd / bwd)
    
    # lab = ifelse(fwd > 0  & bwd > 0,"both",ifelse(fwd > 0,"fwd","bwd"))
    
    stats = DataFrame("f"=fwd,"r"=bwd,"fpos"=fpos,"rpos"=bpos,"d"=d,
                 "u"=u,"ARC"=arc,"URC"=urc,"FSR"=fsr,"M"=M,"A"=A)
                 # "label"=lab)
    stats
}


##' Calculates the quality parameters of one iteration
##'  
##' This function samples \code{nregions} rows from the stat matrix and 
##' fits the linear model \code{lm(d ~ 0 + u + w)}
##'  
##' @param i a numeric value indicating the current iteration.
##' @param stats a \code{data.table} object with the response and covariates for
##' the model
##' @param nregions a numeric value indicating the number of regions sampled.
##' 
##' @return a \code{data.table} with both parameters and some extra info
##' 
##' @export
##' 
##' @rdname .calculate_param_dist
##' @name .calculate_param_dist
##' 
.calculate_param_dist = function(i,stats,nregions)
{
    dt = stats[sample(.N,nregions)]
    model = lm(d ~ 0 + u + w , data = dt)
    data.table(broom::tidy(model))
}