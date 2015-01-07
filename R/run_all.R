

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
  fwd_depths = mcmapply(depth_from_reads,fwd_reads,chr_Nregions,
    SIMPLIFY=FALSE,mc.cores = mc)
  bwd_depths = mcmapply(depth_from_reads,bwd_reads,chr_Nregions,
    SIMPLIFY=FALSE,mc.cores = mc)
  chr_depths = mcmapply("+",fwd_depths,bwd_depths,
    SIMPLIFY=FALSE,mc.cores= mc)

  # Calculate ratio statistics
  message("Calculating ratio statistics")
  labels = mcmapply(function(x,y)ifelse(x >0,ifelse(y > 0,"both","fwd"),"bwd"),
    fwd_depths,bwd_depths,SIMPLIFY = FALSE,mc.cores = mc)
  fwd_strand_ratio = mcmapply("/",fwd_depths,chr_depths,SIMPLIFY=FALSE,mc.cores = mc)  
  depth_width_ratio = mcmapply(function(chr_regions,chr_depth){
    return(chr_depth/width(chr_regions))}
    ,regions,chr_depths,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)

  message("Calculating M vs A values")
  M_values = mcmapply(function(chr_regions,fwd_depth,bwd_depth){
    return(fwd_depth * bwd_depth / width(chr_regions)^2)},regions,
    fwd_depths,bwd_depths,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)
  A_values = mcmapply("/",fwd_depths,bwd_depths,MoreArgs= list(),
    SIMPLIFY=FALSE,mc.cores =mc)  
  M_values = mcmapply(log2,M_values,SIMPLIFY=FALSE,mc.cores = mc)
  A_values = mcmapply(log2,A_values,SIMPLIFY=FALSE,mc.cores=mc)  

  plots = list()
  plots[[1]] = filter_regions_plot(lowerBounds,chr_depths,fwd_strand_ratio,"fwd/(fwd + bwd)",mc)
  plots[[2]] = filter_regions_plot(lowerBounds,chr_depths,depth_width_ratio,"depth/width",mc,TRUE)
  plots[[3]] = filter_label_plot(lowerBounds,chr_depths,labels,mc)
  plots[[4]] = filter_MA_plot(lowerBounds,chr_depths,M_values,A_values,mc)
  
  out = list(plots = plots,regions = regions,depth = depth,boundRegions = table(plots[[3]]$data),
    subset_reads = subset_reads)

  return(out)
}
