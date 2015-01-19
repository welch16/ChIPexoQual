

.separate_gr <- function(gr,st)
{
  gr = .data.table.GRanges(gr)
  setkey(gr,strand)
  gr = gr[strand == st,]
  gr = GRanges(seqnames = gr$seqnames,ranges =
    reads2IRanges(gr),strand = st)
  gr = sort(gr)
  gr = split(gr,as.character(seqnames(gr)))
  return(gr)
}

.add_chr <- function(chr,summary_table)
{
  return(cbind(data.table(chrID=chr),summary_table))
}

#' This functions does all the bound analysis. Briefly, it load the
#' reads, creates the regions and returns a collection of plots that
#' summarizes how the regions vary as they are filtered by depth
#'
#' @param exofile String that contains complete name of the reads file
#'
#' @param mc An integer indicating the number of cores used
#'
#' @param lowerBounds An integer vector of the bounds used in the analysis
#'
#' @export
bound_analysis <- function(exofile,mc,
  lowerBounds = c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,75,100,125,150,200,250,500,750))
{
  # Load reads
  param = ScanBamParam( what = "mapq")
  message("Loading reads")
  reads = readGAlignmentsFromBam(exofile,param = param)
  depth = length(reads)  

  # Convert to GRanges and separate it by strand and seqnames

  message("Separating reads")
  gr = as(reads,"GRanges")
  grF = .separate_gr(gr,"+")
  grR = .separate_gr(gr,"-")

  message("Creating regions")
  regions = create_regions(gr,1)

  # Extract reads per-region
  message("Calculating reads - regions overlaps")
  all_chr= sort(names(seqlengths(gr)),index.return=TRUE)
  regions = regions[all_chr$ix]
  all_chr = all_chr$x
  fwd_overlaps = mcmapply(reads_overlaps,all_chr,
    regions,grF,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
  bwd_overlaps = mcmapply(reads_overlaps,all_chr,
    regions,grR,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

  message("Extracting reads")
  fwd_reads = mcmapply(extract_reads,all_chr,grF,fwd_overlaps,
    regions, MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
  bwd_reads = mcmapply(extract_reads,all_chr,grR,bwd_overlaps,
    regions, MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
  subset_reads = mcmapply(function(x,y)rbind(x,y),
    fwd_reads,bwd_reads,MoreArgs = list(),SIMPLIFY = FALSE,
    mc.cores = mc)

  # Calculate depths per region
  chr_Nregions = mclapply(regions,length,mc.cores = mc)

  message("Calculating region depths")
  fwd_depths = mclapply(fwd_reads,
    function(reads)reads[,length(seqnames),by=region],mc.cores = mc) 
  bwd_depths = mclapply(bwd_reads,
    function(reads)reads[,length(seqnames),by=region],mc.cores = mc)

  # Calculate number of positions
  message("Calculating region's number of unique positions")
  fwd_npos = mclapply(fwd_reads,
    function(reads)reads[,length(unique(start)),by=region],
    mc.cores = mc)
  bwd_npos = mclapply(bwd_reads,
    function(reads)reads[,length(unique(end)),by=region],
    mc.cores = mc)
  widths = mclapply(regions,function(region){
    data.table(region = 1:length(region),V1 = width(region))},
    mc.cores = mc)
  
  message("Merging base statistics")  
  summary_tables = mclapply(chr_Nregions,
    .create_summary,mc.cores = mc)
  name_variables = c("width","fwd_npos","bwd_npos","fwd_depth","bwd_depth")
  variables = list(widths,fwd_npos,bwd_npos,fwd_depths,bwd_depths)
  names(variables) = name_variables
  for(i in name_variables){
    summary_tables = mcmapply(.add_variable,
      summary_tables,variables[[i]],
      MoreArgs = list(name = i),
      SIMPLIFY=FALSE,mc.cores = mc)
  }

    
  # Calculates region depth and number of unique positions
  message("Calculating summary statistics")
  summary_tables = mclapply(summary_tables,
    summary_statistics,mc.cores = mc)
  
  message("Calculating M vs A values")
  summary_tables = mclapply(summary_tables,
    MA_values,mc.cores = mc)

  summary_tables = mcmapply(.add_chr,names(summary_tables),
    summary_tables,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores=mc)

  message("Filtering summary statistics for plots")
  filtered_summary = lapply(lowerBounds,function(x,summary_tables){
    print(x)
    filtered = mclapply(summary_tables,
      function(summary_table,lower){
      return(.fn_filter(lower,summary_table,"depth"))}         
      ,x,mc.cores = mc)
    return(do.call(rbind,filtered))},summary_tables)

  plots = list()
  plots[[1]] = filter_regions_plot(lowerBounds,filtered_summary,
    "prob","fwd/(fwd + bwd)",mc)
  plots[[2]] = filter_regions_plot(lowerBounds,filtered_summary,"dw_ratio","depth/width",mc,TRUE)
  plots[[3]] = filter_label_plot(lowerBounds,filtered_summary,mc)
  plots[[4]] = filter_regions_plot(lowerBounds,filtered_summary,
     "pbc","npos/depth",mc)
  plots[[5]] = filter_MA_plot(lowerBounds,filtered_summary,mc)
  plots[[6]] = positions_reads_map(do.call(rbind,summary_tables),mp=100)
  
  out = list(plots = plots,regions = regions,depth = depth,boundRegions = table(plots[[3]]$data),subset_reads = subset_reads,summary_stats = summary_tables)

  return(out)
}
