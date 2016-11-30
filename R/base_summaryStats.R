##' @importFrom broom tidy
##' @importFrom IRanges start end IRanges
##' @importFrom stats lm
NULL

## Calculate the summary statistics for each region
## 
## @param region a IRanges object with the regions of a certain chromosome
## @param fwdReads a GRanges object with the forward reads of the ChIP-exo 
## experiment
## @param revReads a GRanges object with the reverse reads of the ChIP-exo 
## experiment
## 
## @return stats a \code{DataFrame} object with the summary statistics 
## necessary for the \code{ExoData}.
## 
## @rdname calculateSummary
## @name calculateSummary
## 
calculateSummary <- function(region,fwdReads,revReads)
{
    ## fix formats and stuff
    width(fwdReads) <- 1
    width(revReads) <- 1
    
    fwdReads <- ranges(fwdReads)
    revReads <- ranges(revReads)
    
    fwd <- countOverlaps(region,fwdReads)
    rev <- countOverlaps(region,revReads)
    
    fpos <- IRanges(unique(start(fwdReads)),width = 1)
    rpos <- IRanges(unique(end(revReads)),width = 1)
    fpos <- countOverlaps(region,fpos)
    rpos <- countOverlaps(region,rpos)
    w <- width(region)    
    
    d <- fwd + rev
    u <- fpos + rpos
    
    arc <- d /w 
    urc <- u / d 
    fsr <- fwd / d
    
    M <- log2(fwd) + log2(rev) - 2 * log2(w)
    A <- log2(fwd / rev)

    stats <- DataFrame("fwdReads"=fwd,"revReads"= rev,
                      "fwdPos"=fpos,"revPos"=rpos,
                      "depth"=d,"uniquePos" = u,
                      "ARC"=arc,"URC"=urc,
                      "FSR"=fsr,"M"=M,"A"=A)
    stats
}


##' calculateParamDist
##'  
##' \code{calculateParamDist} calculates the quality parameters of one iteration.
##' This function samples \code{nregions} rows from the stat matrix and fits 
##' the linear model \code{lm(d ~ 0 + u + w)}
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
##' @rdname calculateParamDist
##' @name calculateParamDist
## 
calculateParamDist <- function(i,stats,nregions)
{
    dt <- stats[sample(.N,nregions)]
    model <- lm(depth ~ 0 + uniquePos + width , data = dt)
    data.table(broom::tidy(model))
}