

#' This function calculates the number of reads that overlaps a set of regions
#'
#' @param chr_reads A data.table containing the reads that overlaps the set of regions, with a region key
#'
#' @param chr_nregions, An integer with the expected size of the vector, if it is NULL is going to define as the possible number of regions in chr_reads
#' @return A numeric vector with the depth for each region
#'
#' @export
depth_from_reads <- function(chr_reads,chr_nregions= NULL)
{
  if(is.null(chr_nregions)){    
    chr_depth = table(chr_reads$region)

  }else{
    chr_depth = table(factor(
      chr_reads$region,levels = 1:chr_nregions))     
  }
  dimnames(chr_depth) = NULL
  return(chr_depth)
}


