

##' Calculates the FSR distribution of all regions with depth > value
##'  
##' Calculates the FSR distribution of all regions with depth > value
##'  
##' @param value a numeric value with the depth lower bound.
##' @param DT a \code{data.table} with depth and FSR.
##' @param quantiles a numeric vector with the quantiles used to summarize
##' the FSR distribution.
##' 
##' @return a \code{data.table} with three columns depth, quantiles and FSR.
##' 
##' @rdname .filter_quantiles
##' @name .filter_quantiles
##' 
.filter_quantiles = function(value,DT,quantiles)
{
    dist = DT[d > value,quantile(FSR,probs = quantiles)]
    dist = data.table(d = value ,quantiles , FSR = dist)
    dist
}    
