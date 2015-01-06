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

where:
   f = Number of forward reads that overlap a given region
   r = Number of reverse reads that overlap a given region
   width = Width of the region

  and plots the variation of this statistics for all islands respect different lower bounds.

  Arguments:

    -- exofile

      This is the ChIPexo file of reads

    -- outdir

      This is the directory where all the output files are going to be saved

    -- mcores

      This is the number of cores to be used with parallel programming
 
  Author:

    Rene Welch, Department of Statistics, University of Wisconsin - Madison

");q()}


library(devtools)
library(scales)
library(RColorBrewer)
library(data.table)

codeDr = "/p/keles/ChIPexo/volume4/ChIPexoQual"
load_all(codeDr)


## dataDr = "/p/keles/ChIPexo/volume3/SerrandourData/human"
## dataDr = "/p/keles/ChIPexo/volume3/RenData/BAMfiles"
## dataDr = "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"


exofiles = "edsn1310_042814_qc.sorted.bam" 


# Cell line: MCF7
# TF : ER1

## #param = NULL
## param = ScanBamParam( what = "mapq")
## mc =  8


## exofiles = paste0("ERR3369",c(33,50,58),".bam")
## chipseqfiles = paste0("ERR3369",c(41,36,52,54,59,48,53),".bam")
## names(exofiles) = paste0("rep",c("A","B","C"))
## names(chipseqfiles) = c(paste0("rep", rep(c("A","B","C"),each=2),"_",rep(1:2,3)),"Input")

## exofiles = paste0("AY",552:554,".R1.sort.bam")
## exofiles = "edsn1310_042814_qc.sorted.bam" 


## reads = readGAlignmentsFromBam(file.path(dataDr,exofiles[1]),param = param)
## depth = length(reads)

## ## Ren H3k27ac depth 29,599,796
## ## Ren H3k4me3  28,794,319
## ##       rep2 31,818,368

## ## Serrandour repA 9,289,835
## ##            repB 11,041,833 

## ## Landick data edsn1310 3,909,669


## gr = as(reads,"GRanges")
## grF = subset(gr,strand == "+")
## grR = subset(gr,strand == "-")

## # experiment to define the optimal region

## regions = create_regions(gr,1)

## #regions = mclapply(lowerBounds,function(x)create_regions(gr,x),mc.cores = mc)

## # the idea is to study with basic properties of the islands which is the best
## # strategy to build them

## ## nRegions = data.table(bounds = lowerBounds,
## ##   nr = sapply(regions,function(x)sum(sapply(x,length))))

## grF = sort(grF)
## grR = sort(grR)
## grF = split(grF,seqnames(grF))
## grR = split(grR,seqnames(grR))


## all_chr = names(seqlengths(gr))
## fwd_overlaps = mcmapply(reads_overlaps,all_chr,regions,grF,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
## bwd_overlaps = mcmapply(reads_overlaps,all_chr,regions,grR,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

## fwd_reads = mcmapply(extract_reads,all_chr,grF,fwd_overlaps,regions, MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
## bwd_reads = mcmapply(extract_reads,all_chr,grR,bwd_overlaps,regions, MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
## subset_reads = mcmapply(function(x,y)rbind(x,y),fwd_reads,bwd_reads,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

## chr_Nregions = mclapply(regions,length,mc.cores = mc)
## fwd_depths = mcmapply(depth_from_reads,fwd_reads,chr_Nregions,SIMPLIFY=FALSE,mc.cores = mc)
## bwd_depths = mcmapply(depth_from_reads,bwd_reads,chr_Nregions,SIMPLIFY=FALSE,mc.cores = mc)
## chr_depths = mcmapply("+",fwd_depths,bwd_depths,SIMPLIFY=FALSE,mc.cores= mc)

## fwd_strand_ratio = mcmapply("/",fwd_depths,chr_depths,SIMPLIFY=FALSE,mc.cores = mc)
## labels = mcmapply(function(x,y)ifelse(x >0,ifelse(y > 0,"both","fwd"),"bwd"),fwd_depths,bwd_depths,SIMPLIFY = FALSE,mc.cores = mc)
## depth_width_ratio = mcmapply(function(chr_regions,chr_depth){
##   return(chr_depth/width(chr_regions))}
##   ,regions,chr_depths,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)

## M_values = mcmapply(function(chr_regions,fwd_depth,bwd_depth){
##   return(fwd_depth * bwd_depth / width(chr_regions)^2)},regions,fwd_depths,bwd_depths,
##   MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)
## A_values = mcmapply("/",fwd_depths,bwd_depths,MoreArgs= list(),
##   SIMPLIFY=FALSE,mc.cores =mc)

## M_values = mcmapply(log10,M_values,SIMPLIFY=FALSE,mc.cores = mc)
## A_values = mcmapply(log10,A_values,SIMPLIFY=FALSE,mc.cores=mc)


## plots = list()
## plots[[1]] = filter_regions_plot(lowerBounds,chr_depths,fwd_strand_ratio,"fwd/(fwd + bwd)")
## plots[[2]] = filter_regions_plot(lowerBounds,chr_depths,depth_width_ratio,"depth/width",TRUE)
## plots[[3]] = filter_label_plot(lowerBounds,chr_depths,labels)
## plots[[4]] = filter_MA_plot(lowerBounds,chr_depths,M_values,A_values)















