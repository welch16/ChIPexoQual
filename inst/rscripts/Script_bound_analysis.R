#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_bound_analysis.R - Given a list of fixed lower bounds:

lowerBounds = c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,75,100,125,150,200,250,500,750)

    this script, calculates the following summary statistics:

    fwd_strand_ratio = f / (f + r)
    depth_width_ratio = (f+r) / width
    M = log(f*r/ width^2)
    A = log(f/r)
    pbc = npos / (f + r)

where:
   f = Number of forward reads that overlap a given region
   r = Number of reverse reads that overlap a given region
   width = Width of the region
   npos = Number of unique positions in a given region

  and plots the variation of this statistics for all islands respect different lower bounds.

  Arguments:

    -- exofile

      This is the ChIPexo file of reads

    -- outdir

      This is the directory where all the output files are going to be saved

    -- filecodename

      This is the code name used to save all the files that correspond to this run

    -- mcores

      This is the number of cores to be used with parallel programming
 
  Author:

    Rene Welch, Department of Statistics, University of Wisconsin - Madison

");q()}


stopifnot(length(args)==4)

exofile = args[1]
outdir = args[2]
filecodename = args[3]
mcores = args[4]

codedir = "/p/keles/ChIPexo/volume4/ChIPexoQual"
library(devtools)
library(scales)
library(hexbin)
library(RColorBrewer)

load_all(codedir)

## exofile = "/p/keles/ChIPexo/volume3/CarrollData/human/ERR336933.bam"
## outdir = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human/ER-Rep1new"
## filecodename = "ER-Rep1"
## mcores = 8

message("File: ",exofile)
message("Dir to save: ",outdir)
message("Code name: ",filecodename)
message("Number of cores: ",mcores)


message("Performing analysis")
results = bound_analysis(exofile,mc=as.numeric(mcores))

# Create folders
message("Creating directories")
figsdir = file.path(outdir,"figs")
datadir = file.path(outdir,"data")

if(!file.exists(outdir))dir.create(outdir)
if(!file.exists(figsdir)) dir.create(figsdir)
if(!file.exists(datadir))dir.create(datadir)

# Make plots
plots = results$plots

pdf(file = file.path(figsdir,paste0(filecodename,"_bound_VS_fwd_strand_ratio.pdf")),height = 8,width = 12)
u = print(plots[[1]])
dev.off()

pdf(file = file.path(figsdir,paste0(filecodename,"_bound_VS_depth_width__ratio.pdf")),height = 8,width = 12)
u = print(plots[[2]])
dev.off()

pdf(file = file.path(figsdir,paste0(filecodename,"_bound_VS_label_area.pdf")),height = 8,width = 12)
u = print(plots[[3]])
dev.off()

pdf(file = file.path(figsdir,paste0(filecodename,"_bound_VS_npos_depth_ratio.pdf")),height = 8,width = 12)
u = print(plots[[4]])
dev.off()


pdf(file = file.path(figsdir,paste0(filecodename,"_MA_plots.pdf")),height = 12,width = 12)
u = print(plots[[5]])
dev.off()


pdf(file = file.path(figsdir,paste0(filecodename,"_npos_VS_depth_map.pdf")),height = 8,width = 14)
u = print(plots[[6]])
dev.off()

pdf(file = file.path(figsdir,paste0(filecodename,"_dwRatio_VS_nposDepthRatio_plot.pdf")),height = 8,width = 8)
u = print(plots[[7]])
dev.off()

ggsave(
  file.path(figsdir,paste0(filecodename,"_bound_VS_fwd_strand_ratio.png")),
  plots[[1]],height = 4,width = 6,dpi=600)

ggsave(
  file = file.path(figsdir,paste0(filecodename,"_bound_VS_depth_width__ratio.png")),
  plots[[2]],height = 4,width = 6,dpi=600)

ggsave(
  file = file.path(figsdir,paste0(filecodename,"_bound_VS_label_area.png")),
  plots[[3]],height = 4,width = 6,dpi=600)

ggsave(
  file = file.path(figsdir,paste0(filecodename,"_bound_VS_npos_depth_ratio.png")),
  plots[[4]],height = 4,width = 6,dpi=600)

ggsave(
  file = file.path(figsdir,paste0(filecodename,"_MA_plots.png")),
  plots[[5]],height = 7,width = 7,dpi=1200)

ggsave(
  file = file.path(figsdir,paste0(filecodename,"_npos_VS_depth_map.pdf")),
  plots[[6]],height = 8,width = 14,dpi=1000)

ggsave(
  file = file.path(figsdir,paste0(filecodename,"_dwRatio_VS_nposDepthRatio_plot.pdf")),
  plots[[7]],height = 8,width = 8,dpi = 1000)


# Save data
regions = results$regions
depth = results$depth
boundRegionsTable = results$boundRegions
reads_table = results$subset_reads
summary_stats = results$summary_stats

save(list = "regions",file = file.path(datadir,paste0(filecodename,"_regions.RData")))
save(list = "depth",file = file.path(datadir,paste0(filecodename,"_depth.RData")))
save(list = "boundRegionsTable",file = file.path(datadir,paste0(filecodename,"_boundRegionsTable.RData")))
save(list = "reads_table",file = file.path(datadir,paste0(filecodename,"_reads_by_region.RData")))
save(list = "summary_stats",file = file.path(datadir,paste0(filecodename,"_summary_statistics.RData")))
