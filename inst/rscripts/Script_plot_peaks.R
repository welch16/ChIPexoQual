#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_plot_peaks.R - Given a directory and a filecodename
    this script plots the peaks of a experiment with the
    following statistics

    fwd_strand_ratio = f / (f + r)
    depth_width_ratio = (f+r) / width
    M = log(f*r/ width^2)
    A = log(f/r)

where:
   f = Number of forward reads that overlap a given region
   r = Number of reverse reads that overlap a given region
   width = Width of the region
   npos = Number of unique positions in a given region

  Arguments:

    -- indir

      This is the directory where all the output files are going to be saved

    -- filecodename

      This is the code name used to save all the files that correspond to this run

    -- chr

      This is the chromosome for which the peaks are generated

    -- mcores

      This is the number of cores to be used with parallel programming
 
  Author:

    Rene Welch, Department of Statistics, University of Wisconsin - Madison

");q()}


stopifnot(length(args)==4)

indir = args[1]
filecodename = args[2]
chr = args[3]
mcores = args[4]

codedir = "/p/keles/ChIPexo/volume4/ChIPexoQual"
library(devtools)
library(scales)
library(hexbin)
library(grid)
library(data.table)
library(RColorBrewer)

load_all(codedir)

## indir = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human"
## indir = "/p/keles/ChIPexo/volume3/Analysis/Ren"
## filecodename = "ER-rep1"
## filecodename = "H3k27ac"
## chr = "chr1"
## mcores = 24

indir = file.path(indir,filecodename)

message("Loading directory: ",indir)
message("Code name: ",filecodename)
message("Number of cores: ",mcores)

figsdir = file.path(indir,"figs")
datadir = file.path(indir,"data")

message("Loading statistics, regions and reads")
load(file.path(datadir,paste0(filecodename,"_regions.RData")))
load(file.path(datadir,paste0(filecodename,"_summary_statistics.RData")))
load(file.path(datadir,paste0(filecodename,"_reads_by_region.RData")))
load(file.path(datadir,paste0(filecodename,"_depth.RData")))

# regions,depth,reads_table and summary_stats <-- load

# Get parameters
summary_stats = do.call(rbind,summary_stats)
#chr = names(regions)

# Filter statistics
summary_stats = summary_stats[chrID == chr]
summary_stats = summary_stats[npos > 2] # removing plots that we know how they look like
message("Filtering region with depth < ",median(summary_stats$depth))
summary_stats = summary_stats[depth > median(depth)] # this is a conservative filter according to MA plots
summary_stats = summary_stats[label == "both"]
setkey(summary_stats,region,chrID) # keys are region and chrID

coverage_from_reads <- function(chr,chr_reads,chr_which,strand,mc)
{
  message("Calculating ",ifelse(strand == "+","forward","backward"),
          " coverage for ",chr," islands")
  cover = mclapply( chr_which[[1]],function(i,reads,orientation){
    region_reads = reads[region == i & strand == orientation]
    iranges = reads2IRanges(region_reads)   
    return(coverage(iranges))
    },chr_reads,strand,mc.cores = mc)
  return(cover)
}

chr_which = lapply(chr,function(x,summary_stats){
  summary_stats[chrID == x,.(region)]},summary_stats)
names(chr_which) = chr

reads_table = list(reads_table[[chr]])
names(reads_table) = chr

## reads_table = mcmapply(function(reads,ww){
##   reads[region %in% ww$region]},
##   reads_table,chr_which,mc.cores = mcores,SIMPLIFY=FALSE)


fwd_cover <- mapply(coverage_from_reads,chr,reads_table,
  chr_which,MoreArgs = list("+",mcores),SIMPLIFY=FALSE)
bwd_cover <- mapply(coverage_from_reads,chr,reads_table,
  chr_which,MoreArgs = list("-",mcores),SIMPLIFY=FALSE)

Rle2data.table <- function(rle_data)
{
  coord = cumsum(runLength(rle_data))
  counts = runValue(rle_data)
  return(data.table(coord=coord,counts = counts))
}

normalize.tagcounts <- function(counts,depth)return(counts*1e6 / depth)


plot_cover <- function(fwd_cover,bwd_cover,region_start,region_end,
                       depth,annot = NULL,chr,ext)
{
  
  fwd = Rle2data.table(fwd_cover)
  bwd = Rle2data.table(bwd_cover)

  fwd$strand = "fwd"
  bwd$strand = "bwd"

  M = max(c(fwd$count,bwd$count))*1.2
  
  fwd$count = normalize.tagcounts(fwd$count,depth)
  bwd$count = -normalize.tagcounts(bwd$count,depth)

  annot =   paste0(names(annot),": ",round(annot,4),"\n")
  annot = lapply(annot,function(x)x) 
  annot = do.call(paste0,annot)

  extra = grobTree(textGrob(annot, x=0.1,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=12)))
  
  p = ggplot(data.frame(x=0,y=0))+
    geom_step(data = fwd,aes(x=coord,y=count),colour = "red",  direction = "vh")+
    geom_step(data = bwd,aes(x=coord,y=count),colour = "blue",direction = "hv")+
    theme(legend.position = "none")+annotation_custom(extra)+
    geom_abline(slope = 0,intercept = 0,linetype = 2,colour = I("black"))+
    scale_y_continuous( limits = 1.2*1e6/depth*M * c(-1,1))+ylab("Normalized counts")+
   scale_x_continuous(limits = c(region_start - ext,region_end +ext))+
   ggtitle(paste0(chr,":",region_start,"-",region_end))

  return(p)
}


filtered_regions = mcmapply(function(chr_reg,ww)chr_reg[ww[[1]]],
  list(regions[[chr]]),chr_which,SIMPLIFY=FALSE,mc.cores = mcores)
names(filtered_regions) = chr

plot_regions <- function(chr,chr_fwd,chr_bwd,chr_regions,summary_stats,depth){
  message("Generating peak plots for ",chr)
  starts = start(chr_regions)
  ends = end(chr_regions)
  summaries = summary_stats[chrID == chr, .(depth,npos,prob,dw_ratio,pbc,M,A)]
  ext = 25
  peaks = mclapply(1:nrow(summaries),function(i){
    plot_cover(chr_fwd[[i]],chr_bwd[[i]],starts[i]-ext,ends[i]+ext,
      depth,summaries[i],chr,ext)},mc.cores = mcores )
  return(peaks)
}

peaks = plot_regions(chr,fwd_cover[[chr]],bwd_cover[[chr]],filtered_regions[[1]],summary_stats,depth)

## rf = colorRampPalette(rev(brewer.pal(11,"Spectral")))
## r = rf(16)
## ggplot(summary_stats,aes(dw_ratio,pbc,colour = npos))+geom_point()+scale_color_gradientn(colours = r)

## message("Saving peaks to disc")
## save(list = "peaks",file = file.path(datadir,paste0(filecodename,"-",chr,"_peaks.RData")))



## summary_stats$idx = 1:nrow(summary_stats)
idd = sort(summary_stats$npos,decreasing =TRUE,index.return=TRUE)

message("Plotting peaks")
pdf(file = file.path(figsdir,paste0(filecodename,"-",chr,"_peaks.pdf")))
u = lapply(idd$ix,function(j)print(peaks[[j]]))
dev.off()
