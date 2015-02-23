
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
  summary_table[,pbc:=npos/depth]
  return(summary_table)
}

#' Update the summary_table including the MA-plot values
#'
#' @param summary_table A data.table with summary statistics
#'
#' @return Updated summary table including M and A values
#'
#' @export
MA_values <- function(summary_table)
{
  summary_table[,M:= log2(fwd_depth/bwd_depth)]
  summary_table[,A:=log2(fwd_depth) + log2(bwd_depth)
                    - 2*log2(width)]
  return(summary_table)
}
  
#' Returns a list (or vector when simplify=TRUE) in the following manner:
#' Let ‘i’ be the indices in 'SHIFT', ‘X_i = window(X, 1, length(X) - SHIFT[i])’,
#' and ‘Y_i = window(Y, 1 + SHIFT[i], length(Y) )’. Calculates the
#' set of ‘FUN(X_i, Y_i, ...)’ values and return the results in a convenient form
#'
#' @param SHIFT A non-negative integer vector of shift values.
#'
#' @param FUN The function, found via ‘match.fun’, to be applied to each set of
#' shifted vectors.
#'
#' @param X The Vector or R vector objects to shift.
#'
#' @param Y The Vector or R vector objects to shift.
#'
#' @param ...  Further arguments for FUN
#'
#' @param simplify A logical value specifying whether or not the result should be
#' simplified to a vector
#'
#' @export
mcshiftapply <- function(SHIFT,FUN,X,Y,mc,...,simplify=TRUE)
{
    
  FUN = match.fun(FUN)
  N =  length(X)
  
  stopifnot(length(Y) == N)
  stopifnot(length(SHIFT) > 0)
  
  infoX = getStartEndRunAndOffset(X, rep.int(1L ,length(SHIFT)), N - SHIFT)
  runStartX = infoX[["start"]][["run"]]
  offsetStartX = infoX[["start"]][["offset"]]
  runEndX = infoX[["end"]][["run"]]
  offsetEndX = infoX[["end"]][["offset"]]

  infoY = getStartEndRunAndOffset(Y, 1L + SHIFT, rep.int(N, length(SHIFT)))  
  runStartY = infoY[["start"]][["run"]]
  offsetStartY = infoY[["start"]][["offset"]]
  runEndY = infoY[["end"]][["run"]]
  offsetEndY = infoY[["end"]][["offset"]]
  
  newX = new("Rle")
  newY = new("Rle")
  
  out = mclapply(seq_len(length(SHIFT)), function(i){
    FUN(
      .Call2("Rle_window",X, runStartX[i],
        runEndX[i], offsetStartX[i], offsetEndX[i], 
        newX, PACKAGE = "S4Vectors"),
      .Call2("Rle_window", Y, runStartY[i], runEndY[i],
        offsetStartY[i], offsetEndY[i], newY, PACKAGE = "S4Vectors"))}
        ,...,mc.cores = mc)
  if(simplify)out = do.call(c,out)
  
  return(out)
}


#' Calculate the PCR bottleneck coefficient for a data set
#'
#' @param reads_table List of data.table with the experiment's separated reads by each region
#'
#' @param mcores Numeric value with the number of cores to use
#'
#' @export
pbc <- function(reads_table,mcores)
{
  fwd_tabs = mclapply(reads_table,function(x){
    table(x[strand == "+",(start)])},mc.cores = mcores)
  bwd_tabs = mclapply(reads_table,function(x){
    table(x[strand == "-",(end)])},mc.cores = mcores)
  nUniq = sum(do.call(c,
    mclapply(fwd_tabs,function(x)sum(x==1),mc.cores =mcores))) +
      sum(do.call(c,
    mclapply(bwd_tabs,function(x)sum(x==1),mc.cores = mcores)))
  nTotal = sum(do.call(c,
    mclapply(fwd_tabs,function(x)sum(x>=1),mc.cores =mcores))) +
      sum(do.call(c,
    mclapply(bwd_tabs,function(x)sum(x>=1),mc.cores = mcores)))
  return(nUniq / nTotal)       
}
