#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_compare_peaks.R - Given a directory and a filecodename this
    script compare the called peaks of an experiment with the a set of
    regions with the following statistics

    fwd_strand_ratio = f / (f + r)
    depth_width_ratio = (f+r) / width
    M = log(f*r/ width^2)
    A = log(f/r)
    npos_depth_ratio = npos / (f +r)

where:
   f = Number of forward reads that overlap a given region
   r = Number of reverse reads that overlap a given region
   width = Width of the region
   npos = Number of unique positions in a given region

  Arguments:

    -- indir

      This is the directory where all the output files are going to be saved

    -- filecodename

      This is the code name used to save all the files that correspond
      to this run

    -- mcores

      This is the number of cores to be used with parallel programming

    -- peaksfilename

      This is the name of the file that contains the peaks, this is a
      data.table containing at least for columns

      chrID | peakStart | peakStop | score
       
  Author:

    Rene Welch, Department of Statistics, University of Wisconsin - Madison

");q()}


stopifnot(length(args) %in% 3:4)


indir = args[1]
filecodename = args[2]
mcores = args[3]
peaksfilename = args[4]
if(is.na(peaksfilename))
{
  peaksfilename = file.path(indir,filecodename,"peaks",paste0(filecodename,"-peaks.RData"))
}




codedir = "/p/keles/ChIPexo/volume4/ChIPexoQual"
library(devtools)
library(scales)
library(hexbin)
library(grid)
library(data.table)
library(RColorBrewer)

load_all(codedir)


## indir = "/p/keles/ChIPexo/volume3/Analysis/Ren"
## filecodename = "H3k27ac"
## mcores = 8
## peaksfilename = file.path(indir,filecodename,"peaks",paste0(filecodename,"-peaks.RData"))
## peaksfilename = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human/ER-rep1/peaks/ER-rep1-peaks.RData

indir = file.path(indir,filecodename)

message("Loading directory: ",indir)
message("Code name: ",filecodename)
message("Number of cores: ",mcores)
message("Peaks file: ",peaksfilename)

datadir = file.path(indir,"data")
figsdir = file.path(indir,"figs")
message("Loading statistics, regions and reads")
load(file.path(datadir,paste0(filecodename,"_regions.RData")))
load(file.path(datadir,paste0(filecodename,"_summary_statistics.RData")))
#load(file.path(datadir,paste0(filecodename,"_reads_by_region.RData")))
load(file.path(datadir,paste0(filecodename,"_depth.RData")))
load(peaksfilename)

# Covert peaks to IRangesList
setkey(peaks,chrID)

regions_chr = names(regions)
peaks_chr = unique(peaks[,.(chrID)][[1]])

# Filtering peaks
message("Filtering peaks with regions")
chr = intersect(regions_chr,peaks_chr)
regions = regions[regions_chr %in% chr]
summary_stats = summary_stats[names(summary_stats) %in% chr]
#reads_table = reads_table[names(reads_table) %in% chr]
peaks = peaks[chr]

# Converting peaks to IRangesList
peaks_list = mclapply(chr,function(x){
  out = peaks[x]
  setnames(out,names(out),c("chrID","start","end","score"))
  return(out)}
  ,mc.cores = mcores)

peaks_list = IRangesList(mclapply(peaks_list,
  reads2IRanges,mc.cores = mcores))
names(peaks_list) = chr

overlaps = findOverlaps(peaks_list,regions)
                      
extract_summary_stats <- function(overlap,peaks,summary_stat)
{
  setkey(summary_stat,region)
  out = summary_stat[subjectHits(overlap)]
  out[,peakID := queryHits(overlap)]
  return(out)
}

label_peaks_in_summary <- function(overlap,summary_stat)
{
  setkey(summary_stat,region)
  summary_stat[,isPeak:= as.integer(region %in% subjectHits(overlap))]
  return(summary_stat)  
}

summary_stats = mcmapply(label_peaks_in_summary,overlaps,summary_stats,SIMPLIFY=FALSE,mc.cores = mcores)

summary_stats_all = do.call(rbind,summary_stats)

peak_stats = mcmapply(extract_summary_stats,
  overlaps,peaks_list,summary_stats,SIMPLIFY=FALSE,mc.cores = mcores)

peaks_all = do.call(rbind,peak_stats)

## lambda = 0
## -log10(1 - phyper(nrow(summary_stats_all[isPeak == 1 & npos > lambda]),
##        nrow(summary_stats_all[npos > lambda]),
##        nrow(summary_stats_all[npos <= lambda]),
##        nrow(summary_stats_all[isPeak == 1])))
         

## library(rpart)
## fit = rpart(isPeak ~ npos  +  prob + M + A,method = "class",data = as.data.frame(summary_stats))


## lowerBounds = c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,100,150,200,250,500,750)
## message("Filtering summary statistics for plots")
## filtered_summary = lapply(lowerBounds,function(x,summary_tables){
##   print(x)
##   filtered = mclapply(summary_tables,
##     function(summary_table,lower){
##     return(.fn_filter(lower,summary_table,"depth"))}         
##     ,x,mc.cores = mcores)
##   return(do.call(rbind,filtered))},summary_stats)

## positions_reads_map(peaks_all,150)

## p = filter_regions_plot(lowerBounds,filtered_summary,"prob","fwd/(fwd + bwd)",mcores)

## filter_regions_plot(lowerBounds,filtered_summary,"dw_ratio","depth/width",mcores,TRUE)

## filter_label_plot(lowerBounds,filtered_summary,mcores)

## filter_regions_plot(lowerBounds,filtered_summary,"pbc","npos/depth",mcores)
## filter_scatter_plot(lowerBounds,"A","M",filtered_summary,mcores) + geom_abline(slope=0,intercept = 0,linetype =2)

## filter_scatter_plot(lowerBounds,"dw_ratio","pbc",filtered_summary,mcores) + scale_x_continuous(limits = c(-.1,10))+ylab("npos/depth")+xlab("depth/width")

histograms = list()
histograms[[1]] = ggplot(summary_stats_all,aes(prob))+geom_histogram(aes(y=..density..))+
  facet_grid(isPeak~.)+xlab("fwd strand ratio")
histograms[[2]] = ggplot(summary_stats_all,aes(npos))+geom_histogram(aes(y=..density..))+
  facet_grid(isPeak~.)+xlab("number of unique positions")+
  scale_x_continuous(limits = c(0,1e3))
histograms[[3]] = ggplot(summary_stats_all,aes(depth))+geom_histogram(aes(y=..density..))+
  facet_grid(isPeak~.)+scale_x_continuous(limits = c(0,1000))
histograms[[4]] = ggplot(summary_stats_all,aes(dw_ratio))+geom_histogram(aes(y=..density..))+
  facet_grid(isPeak~.)+scale_x_continuous(limits = c(0,5))+
  xlab("depth/width")
histograms[[5]] = ggplot(summary_stats_all,aes(pbc))+geom_histogram(aes(y=..density..))+
  facet_grid(isPeak~.)+scale_x_continuous(limits = c(0,1.2))+
  xlab("npos/depth")


pdf(file = file.path(figsdir,paste0(filecodename,"_peaks_comp_hist.pdf")))
u = lapply(histograms,print)
dev.off()

lowerBounds = c(0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,100,150,200,250,500,750)
rf = colorRampPalette(rev(brewer.pal(11,"Spectral")))
r = rf(16)

# MA-plots
maplots = mclapply(lowerBounds,
  function(x){
    ggplot(summary_stats_all[label == "both" & depth > x],aes(A,M))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+facet_grid(isPeak~.)+      
      ggtitle(paste0("depth > ",x))},mc.cores =mcores)

# dw_ratio - pbc
plots = mclapply(lowerBounds,
  function(x){
    ggplot(summary_stats_all[label == "both" & depth > x],aes(dw_ratio,pbc))+
      stat_binhex(bins = 70)+ scale_fill_gradientn(colours = r,trans = 'log10')+
        facet_grid(isPeak~.)+ scale_x_continuous(limits = c(0,10))+
        scale_y_continuous(limits = c(0,1))+ggtitle(paste0("depth > ",x))+
        xlab("depth/width")+ylab("npos/depth")},mc.cores = mcores)


pdf(file = file.path(figsdir,paste0(filecodename,"_peaks_MA_plots.pdf")))
u = lapply(maplots,print)
dev.off()

pdf(file = file.path(figsdir,paste0(filecodename,"_peaks_depth_width_npos_plots.pdf")))
u = lapply(plots,print)
dev.off()

