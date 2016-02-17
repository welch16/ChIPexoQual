
##'
##' @export
create_exo_experiment <- function(file,height = 1,calc_summary = FALSE,parallel = FALSE,mc = 0)
{
  stopifnot(file.exists(file))
  stopifnot(is.numeric(height))
  stopifnot(is.logical(calc_summary))

  rr <- readGAlignments(file,param = NULL)
  rr <- as(rr,"GRanges")
  
  out <- new("ChIPexo_experiment",file = file, reads = rr , depth = length(rr),height = height)
  
  if(calc_summary){
    out <- separate_regions(out,parallel , mc)
  }
  
  return(out)
}


##'
##' @export
separate_regions <- function(exo_expmt,parallel,mc)
{  
  exo_reads <- reads(exo_expmt)
  sqnms <- levels(seqnames(exo_reads))
  cover <- coverage(exo_reads)
  regions <- slice(cover,lower = height(exo_expmt),rangesOnly = TRUE)
  regions <- as(regions,"GRanges")
  if(!parallel){
    exo_expmt@summary_stats <- calculate_summary(regions,reads(exo_expmt))
  }else{
    regions_list <- split(regions,as.character(seqnames(regions)))
    reads_list <- reads(exo_expmt)
    reads_list <- split(reads_list,as.character(seqnames(reads_list)))
    stats_list <- mcmapply(calculate_summary,regions_list,reads_list,SIMPLIFY = FALSE,mc.cores = mc)
    exo_expmt@summary_stats <- do.call(rbind,stats_list)
  }  
  return(exo_expmt)
}




calculate_summary <- function(region,reads)
{
  ## fix formats and stuff
  ov <- findOverlaps(region,reads)
  reads <- gr2dt(reads)
  w <- width(region)    
  region <- gr2dt(region)
  region[ , width := w]
  ##
  region[, match := paste0(seqnames,":",start,"_",end)]
  reads[  subjectHits(ov), match := region[queryHits(ov), (match)] ]
  ##
  reads <- reads[!is.na(match)]
  ## get base statistics
  fwd <- reads[strand == "+",
    {
      f = length(start)
      f_pos = length(unique(start))
      list(f = f,f_pos = f_pos)
    },by = match]
  bwd <- reads[strand == "-",
    {
      r = length(end)
      r_pos = length(unique(end))
      list(r = r,r_pos = r_pos)
    },by = match]

  ## merge statistics
  stats <- merge(region,fwd, by = "match",all.x = TRUE)
  stats <- merge(stats,bwd, by = "match",all.x = TRUE)
  stats[is.na(f), f := 0]
  stats[is.na(r), r := 0]
  stats[is.na(f_pos), f_pos := 0]
  stats[is.na(r_pos), r_pos := 0]
  
  ## calculate composite stats
  stats[ , depth := f + r]
  stats[ , npos := f_pos + r_pos]
  stats[ , ave_reads := depth / width]
  stats[ , cover_rate := npos / depth]
  stats[ , fsr := f / (f + r)]

  stats[ , M := as.numeric(NA)]
  stats[ , A := as.numeric(NA)]

  stats[f > 0 & r > 0, M := log2(f) + log2(r) - 2 * log2(width)]
  stats[f > 0 & r > 0, A := log2(f/ r)]

  stats[ , strand := NULL]

  return(stats)
}
