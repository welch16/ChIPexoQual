
##'
##' @export
create_exo_experiment <- function(file,height = 1,calc_summary = FALSE)
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
##' @export
separate_regions <- function(exo_expmt)
{
  exo_reads <- reads(exo_expmt)
  sqnms <- levels(seqnames(exo_reads))
  exo_reads <- split(exo_reads,as.character(seqnames(exo_reads)))
  cover <- lapply(exo_reads,coverage)
  regions <- lapply(cover,slice , lower = height(exo_expmt),rangesOnly = TRUE)
  regions <- mapply(function(x,y)x[[y]],regions,sqnms,SIMPLIFY = FALSE)
  regions <- mapply(function(x,y)GRanges(seqnames = x,ranges = y),sqnms,regions,SIMPLIFY = FALSE)
  regions <- unlist(GRangesList(regions))
  exo_expmt@summary_stats <- calculate_summary(regions,reads(exo_expmt))
  return(exo_expmt)
}




calculate_summary <- function(region,reads,mc = detectCores())
{
  ## fix formats and stuff
  
  ov <- findOverlaps(region,reads)
  reads <- gr2dt(reads)
  w <- width(region)    
  region <- gr2dt(region)
  region[ , width := w]
  region[, match := paste0(seqnames,":",start,"_",end)]
  reads[  subjectHits(ov), match := region[queryHits(ov), (match)] ]  
  reads[,strand := ifelse(strand == "+", "F","R")]
  reads <- reads[!is.na(match)]

  ## get base statistics
  f <- reads[,sum(strand == "F"),by = match]
  setnames(f,names(f),c("match","f"))
  setkey(f,match)
  r <- reads[,sum(strand == "R"),by = match]
  setkey(r,match)
  setnames(r,names(r),c("match","r"))
  f_uniq <- reads[strand == "F",length(unique(start)),by = match]
  setnames(f_uniq,names(f_uniq),c("match","f_pos"))
  setkey(f_uniq,match)
  r_uniq <- reads[strand == "R",length(unique(end)),by = match]
  setnames(r_uniq,names(r_uniq),c("match","r_pos"))
  setkey(r_uniq,match)

  ## merge statistics
  stats <- merge(region,f,by = "match",allow.cartesian = TRUE)
  stats <- merge(stats,r,by = "match",allow.cartesian = TRUE)
  stats <- merge(stats,f_uniq,by = "match",allow.cartesian = TRUE,all = TRUE)
  stats <- merge(stats,r_uniq,by = "match",allow.cartesian = TRUE, all = TRUE)
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

  stats[f > 0 & r > 0, M := log2(f * r) - 2 * log2(width)]
  stats[f > 0 & r > 0, A := log2( f/ r)]

  stats[ , strand := NULL]

  return(stats)
}
