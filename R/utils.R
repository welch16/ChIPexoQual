##' @import data.table
##' @importFrom GenomicRanges GRanges
##' @importFrom IRanges IRanges
NULL

##' Converts a GRanges object to data.table
##'
##' @description Converts a GRanges object to data.table, the columns
##' are gonna have the same names as the GRanges fields
##'
##' @param gr A GRanges object
##'
##' @export
##'
##' @return A data.table object
##'
##' @rdname gr2dt
##' @name gr2dt
##'
##' @seealso \code{\link{dt2gr}}
##' @examples
##'
##' gr <- GRanges(seqnames = "chr1",
##'         ranges = IRanges(start = 3:10,width = 30),
##'         strand = "*")
##' gr2dt(gr)
##'

gr2dt <- function(gr)
{
  out <- data.table(seqnames = as.character(seqnames(gr)),
                    start = start(gr),end = end(gr),
                    strand = as.character(strand(gr)))
  extra <- mcols(gr)
  extra <- data.table(as.data.frame(extra))
  if(nrow(extra)>0){
    out <- cbind(out,extra)
  }
  return(out)
}

##################################################################################

##' Converts a data.table to GRanges
##'
##' @description Coverts a data.table into a GRanges object
##'
##' @param dt A data.table with at least the columns seqnames,start and end 
##'
##' @export
##'
##' @return A GRanges object
##'
##' @rdname dt2gr
##' @name dt2gr
##'
##' @seealso \code{\link{gr2dt}}
##' @examples
##'
##' dt <- data.table(seqnames = rep("chr1",8),start = 3:10,end = sample(11:18,8))
##' dt2gr(dt)
##' 
##' dt2 <- copy(dt); dt2[,extra := runif(8)]
##' dt2gr(dt2)
##'
##' dt3 <- copy(dt); dt3[,strand := sample(c("+","-"),8,replace = TRUE)]
##' dt2gr(dt3)
##' 
##' dt3[ , extra := runif(8)]
##' dt2gr(dt3)

dt2gr <- function(dt)
{
  stopifnot( all(c("seqnames","start","end") %in% names(dt)))

  if("strand" %in% names(dt)){
    st <- dt[,(strand)]
  }else{
    st <- "*"
  }

  gr <- GRanges(seqnames = dt[,(seqnames)],
    ranges <- IRanges(start = dt[,(start)],end = dt[,(end)]),
    strand = st)
  
  if(any(! names(dt) %in% c("seqnames","start","end","strand"))){
      extra <- dt[ ,!names(dt) %in% c("seqnames","start","end","strand"),with = FALSE]
      mcols(gr) <- as.data.frame(extra)    
  }

  return(gr)
 
}

##################################################################################

dt2ir <- function(dt)
{
  stopifnot(all(c("start","end") %in% names(dt)))

  ir <- IRanges(start = dt[,(start)],end = dt[,(end)])
  return(ir)
}
