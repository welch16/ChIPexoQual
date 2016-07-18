##' @importFrom S4Vectors DataFrame
NULL


##' Calculate the summary statistics for each region
##' 
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
    
    lab = ifelse(fwd > 0  & bwd > 0,"both",ifelse(fwd > 0,"fwd","bwd"))
    
    stats = DataFrame("f"=fwd,"r"=bwd,"fpos"=fpos,"rpos"=bpos,"d"=d,
                 "u"=u,"ARC"=arc,"URC"=urc,"FSR"=fsr,"M"=M,"A"=A,
                 "label"=lab)
    return(stats)
}
