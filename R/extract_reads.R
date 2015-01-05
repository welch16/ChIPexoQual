
#' @import data.table
NULL

#' Extract the reads that overlap an island from the whole reads list for a certain chromosome
#'
#' @param chr, This is the chromosome for which we are looking the regions
#'
#' @param chr_reads, GRanges object with a set of reads for certain chromosome
#'
#' @param chr_overlaps, A Hits object with all the overlaps between the reads and the islands or certain chromosome 
#'
#' @param chr_islands, IRanges object with a collection of regions in certain chromosome 
#'
#' @return A list with all the overlaps between reads and island for each chromosome
#'
#' @export
extract_reads <- function(chr,chr_reads,chr_overlaps,chr_islands)
{

  chr_overlaps = .data.table.Hits(chr_overlaps)
  chr_reads = .data.table.GRanges(chr_reads)
  setkey(chr_overlaps,queryHits)
  chr_reads$ID = 1:nrow(chr_reads)
  setkey(chr_reads,ID)
  chr_islands = .data.table.GRanges(GRanges(seqnames = chr,
    ranges = chr_islands,strand = "*")) # check here, possible bug
#  chr_overlaps[queryHits %in% 1:nrow(chr_islands),]
  ov_reads = chr_reads[chr_overlaps[queryHits %in% 1:nrow(chr_islands),list(subjectHits)],]
  ov_reads[,ID := NULL]
  ov_reads[,region:=chr_overlaps[,queryHits]]
  return(ov_reads)
}
