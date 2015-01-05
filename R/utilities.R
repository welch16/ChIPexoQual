
.data.table.Hits <- function(x)data.table(as.data.frame(x))

.data.table.GRanges <- function(x)
{
  dt = data.table(seqnames = as.character( seqnames(x)),
    start = start(x),end = end(x),
    strand = as.character(strand(x)))
  return(dt)
}

reads2IRanges <- function(x)IRanges(start = x$start,end = x$end)
