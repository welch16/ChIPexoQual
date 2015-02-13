
#' @import S4Vectors
NULL

.data.table.Hits <- function(x)data.table(as.data.frame(x))

.data.table.GRanges <- function(x)
{
  dt = data.table(seqnames = as.character( seqnames(x)),
    start = start(x),end = end(x),
    strand = as.character(strand(x)))
  return(dt)
}

.create_summary <- function(n)return(data.table(region = 1:n))
.add_variable <- function(summary_table,value,name)
{
  summary_table[,(name):=as.integer(0)]
  setkey(summary_table,region)
  summary_table[value[,.(region)],(name):=value[,.(V1)]]
  return(summary_table)
}

reads2IRanges <- function(x)IRanges(start = x$start,end = x$end)

Rle2data.table <- function(rle_data)
{
  coord = cumsum(runLength(rle_data))
  counts = runValue(rle_data)
  return(data.table(coord=coord,counts = counts))
}

normalize.tagcounts <- function(counts,depth)return(counts*1e6 / depth)

which.rank <- function(var,k=1,decreasing = TRUE)
{
  vv = sort(var,decreasing = decreasing,index.return =TRUE)
  return(vv$ix[k]) 
}

localMinima <- function(x)
{
  # Use -Inf instead if x is numeric (non-integer)
  y = diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y = cumsum(rle(y)$lengths)
  y = y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}

localMaxima <- function(x) localMinima(-x)

getStartEndRunAndOffset <- function(x, start, end)
{
    .Call2("Rle_getStartEndRunAndOffset", x, start, end, PACKAGE="S4Vectors")
}
