
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

#' Function to convert to data.table
#'
#' @param x, A GRanges object
#'
#' @return A data.table object
#'
#' @export
data.table.GRanges <- function(x)data.table(as.data.frame(x))

#' Function to convert to data.table
#'
#' @param x, A Hits object
#'
#' @return A data.table object
#'
#' @export
data.table.Hits <- function(x)data.table(as.data.frame(x))




#' Extract the reads that overlap an island from the whole reads list for a certain chromosome (or genome)
#'
#' @param chr, This is the chromosome (or genome) for which we are looking the regions
#'
#' @param chr_reads, GRanges object with a set of reads for certain chromosome (or genome)
#'
#' @param chr_overlaps, A Hits object with all the overlaps between the reads and the islands or certain chromosome (or genome)
#'
#' @param chr_islands, IRanges object with a collection of regions in certain chromosome (or genome)
#'
#' @return A list with all the overlaps between reads and island for each chromosome
#'
#' @export
extract_reads <- function(chr,chr_reads,chr_overlaps,chr_islands)
{
  chr_overlaps = data.table.Hits(chr_overlaps)
  chr_reads = data.table.GRanges(chr_reads)
  setkey(chr_overlaps,queryHits)
  chr_reads$ID = 1:nrow(chr_reads)
  setkey(chr_reads,ID)
  chr_islands = data.table.GRanges(GRanges(seqnames = chr,
    ranges = ranges(chr_islands),strand = "*"))
  ov_reads = chr_reads[chr_overlaps[queryHits %in% 1:nrow(chr_islands),list(subjectHits)],]
  ov_reads[,ID := NULL]
  ov_reads[,region:=chr_overlaps[,queryHits]]
  return(ov_reads)
}
