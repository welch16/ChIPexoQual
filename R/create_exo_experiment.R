
##' 
create_exo_experiment <- function(file,height,calc_summary = FALSE)
{
  stopifnot(file.exists(file))
  stopifnot(is.numeric(height))
  stopifnot(is.logical(calc_summary))
  
  rr <- readGAlignments(file,param = NULL)
  rr <- as(rr,"GRanges")
  
  out <- new("ChIPexo_experiment",file = file, reads = rr , depth = length(rr),height = height)
  
  if(calc_summary){
    out <- separate_regions(out)
  }
  
  return(out)
}


##'
separate_regions <- function(exo_expmt)
{
  
}