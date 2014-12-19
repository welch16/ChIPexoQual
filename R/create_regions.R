
#' Separates a set of reads into regions in the genome
#'
#' @param reads, This is a GRanges object with the reads obtained by a ChIP-seq or a ChIP-exo experiment
#'
#' @param lower, This is a positive number that indicates the lower bound used to delimite the islands
#' 
#' @return A set of regions in the genome
#'
#' @export
create_regions <- function(reads,lower)
{
  stopifnot(class(reads)=="GRanges")
  stopifnot(lower > 0)
  cover = coverage(reads)
  islands = slice(cover,lower = lower,rangesOnly=FALSE)
  return(islands)
}



#' Find the overlaps between the reads and regions for certain chromosome
#'
#' @param chr, This is the chromosome (or genome) for which we are looking the regions
#'
#' @param chr_islands, IRanges object with a collection of regions in certain chromosome (or genome)
#'
#' @param chr_reads, GRanges object with a set of reads for certain chromosome (or genome)
#'
#' @return A list with all the overlaps between reads and island for each chromosome
#'
#' @export
reads_overlaps <- function(chr,chr_islands,chr_reads)
{
  chr_overlaps= findOverlaps(GRanges(seqnames = chr,
     ranges = chr_islands,strand = "*"),chr_reads)
  return(chr_overlaps)
}



.data.table.Hits <- function(x)data.table(as.data.frame(x))



.data.table.GRanges <- function(x)
{
  dt = data.table(seqnames = as.character( seqnames(x)),
    start = start(x),end = end(x),
    strand = as.character(strand(x)))
  return(dt)
}


# Extract the reads that overlap an island from the whole reads list for a certain chromosome (or genome)
#
# @param chr, This is the chromosome (or genome) for which we are looking the regions
#
# @param chr_reads, GRanges object with a set of reads for certain chromosome (or genome)
#
# @param chr_overlaps, A Hits object with all the overlaps between the reads and the islands or certain chromosome (or genome)
#
# @param chr_islands, IRanges object with a collection of regions in certain chromosome (or genome)
#
# @return A list with all the overlaps between reads and island for each chromosome
#
# @export
