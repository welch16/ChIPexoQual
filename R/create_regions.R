
#' Separates a set of reads into regions in the genome
#'
#' @param reads, This is a GRanges object with the reads obtained by a ChIP-seq or a ChIP-exo experiment
#'
#' @param lower, This is a positive number that indicates the lower bound used to delimite the islands
#' 
#' @return A set of regions in the genome
#'
#' @export
.create_regions <- function(reads,lower)
{
  stopifnot(class(reads)=="GRanges")
  stopifnot(lower > 0)
  cover = coverage(reads)
  islands = slice(cover,lower = lower,rangesOnly=TRUE)
  return(islands)
}

