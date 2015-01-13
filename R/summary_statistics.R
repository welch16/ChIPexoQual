
#' Update the summary_table with summary statistics
#'
#' @param summary_table A data.table with the fields width, fwd_npos, bwd_npos, fwd_depth and bwd_depth
#'
#' @return Updated summary_table with compound statistics
#'
#' @export
summary_statistics <- function(summary_table)
{
  summary_table[,depth:=fwd_depth+bwd_depth]
  summary_table[,npos :=fwd_npos + bwd_npos]
  summary_table[,label:=ifelse(fwd_depth > 0 & bwd_depth >0,
      "both",ifelse(fwd_depth > 0,"fwd","bwd"))]
  summary_table[,prob:=fwd_depth/depth]
  summary_table[,fwd_dw_ratio:= fwd_depth/width]
  summary_table[,bwd_dw_ratio:= bwd_depth/width]
  summary_table[,dw_ratio:= depth/width]
  return(summary_table)
}

