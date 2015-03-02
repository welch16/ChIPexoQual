#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_plot_peaks.R - Given a directory and a family code name,
    this script compares across all data sets in the family, the
    average read coverage vs unique read coverage rate plot.

    It builds a gold standard set of regions by intersecting the set
    of islands for all samples in the family and visualizing the
    behaviour of the signal to noise plot for the gold standard
    set of regions

  Arguments:

    -- indir

      This is the directory where all the output files are going to be saved

    -- familycodename

      This is the code name used to select all samples of the same
      experiment. The availabel families are:

      - Carroll_human
      - Carroll_mouse
      - Ren_H3k4me1
      - Landick_Sig70_exp
      - Landick_Sig70_stat
      - Landick_SigmaS_exp
      - Landick_SigmaS_stat
      - Landick_Sig70_rif0
      - Landick_Sig70_rif20
      - Landick_Beta_rif0
      - Landick_Beta_rif20
      - Landick_BetaPrimeFlag_rif0
      - Landick_BetaPrimeFlag_rif20

    -- mcores

      This is the number of cores to be used with parallel programming
 
  Author:

    Rene Welch, Department of Statistics, University of Wisconsin - Madison

");q()}


stopifnot(length(args)==3)

indir = args[1]
familycodename = args[2]
mcores = args[3]

codedir = "/p/keles/ChIPexo/volume4/ChIPexoQual"
library(devtools)
library(parallel)
library(data.table)
library(RColorBrewer)

 ## indir = "/p/keles/ChIPexo/volume3/Analysis/Ren"
 ## familycodename = "ren_h3k4me1"
 ## mcores = 16

load_all(codedir)

stopifnot(tolower(familycodename) %in%
  c("carroll_human","carroll_mouse",
    "ren_h3k4me1",
    "landick_sig70_exp","landick_sig70_stat",
    "landick_sigmas_exp","landick_sigmas_stat",
    "landick_sig70_rif0","landick_sig70_rif20",
    "landick_beta_rif0","landick_beta_rif20",
    "landick_betaprimeflag_rif0","landick_betaprimeflag_rif20"))

filecodenames <- function(familycode)
{
  switch(tolower(familycode),
    carroll_human={
      out <- paste0("ER-rep",1:3)
    },
    carroll_mouse={
      out <- paste0("FoxA1-rep",1:3)
    },
    ren_h3k4me1={
      out <- paste0("H3k4me1-rep",1:2)
    },
    landick_sig70_exp={
      out <- paste0("Sig70-exp-rep",1:2)
    },
    landick_sig70_stat={
      out <- paste0("Sig70-stat-rep",1:2)
    },
    landick_sigmas_exp={
      out <- paste0("SigmaS-exp-rep",1:2)
    },
    landick_sigmas_stat={
      out <- paste0("SigmaS-stat-rep",1:2)
    },
    landick_sig70_rif0={
      out <- paste0("Sig70-rif0-rep",1:2)
    },
    landick_sig70_rif20={
      out <- paste0("Sig70-rif20-rep",1:2)
    },
    landick_beta_rif0={
      out <- paste0("Beta-rif0-rep",1:2)
    },
    landick_beta_rif20={
      out <- paste0("Beta-rif20-rep",1:2)
    },
    landick_betaprimeflag_rif0={
      out <- paste0("BetaPrimeFlag-rif0-rep",1:2)
    },
    landick_betaprimeflag_rif20={
      out <- paste0("BetaPrimeFlag-rif20-rep",1:2)
    } )
  return(out)
}

familynames = filecodenames(familycodename)

# load_all_regions

message("Loading regions")
regions = mclapply(familynames,function(x,indir){
  load(file.path(indir,x,"data",paste0(x,"_regions.RData")))
  return(regions)},indir,mc.cores = mcores)

message("Loading reads")
reads = mclapply(familynames,function(x,indir){
  load(file.path(indir,x,"data",paste0(x,"_reads_by_region.RData")))
  return(reads_table)},indir,mc.cores = mcores)

filter_strand <- function(set_reads,st,mcores)
{
  out = mclapply(set_reads,function(x,st){
    setkey(x,strand)
    chr = x[,(seqnames)][1]
    return(GRanges(seqnames = chr,ranges = reads2IRanges(x[st]),
      strand = st))},st,mc.cores = mcores)
  return(out)  
}

readsF = lapply(reads,filter_strand,"+",mcores)
readsR = lapply(reads,filter_strand,"-",mcores)


readLengthF = sapply(readsF,function(x,mcores){
  out = mclapply(x,function(y)mean(width(y)),mc.cores = mcores)
  return(out)},mcores)
readLengthR = sapply(readsR,function(x,mcores){
  out = mclapply(x,function(y)mean(width(y)),mc.cores = mcores)
  return(out)},mcores)

readLength = round(mean(c(unlist(readLengthF),unlist(readLengthR))),0)

message("Creating gold standard regions by intersection")
message("Filtering the regions with width <= ",readLength)
if(length(regions) == 2){
  gold_regions = mcmapply(function(...,K){
    chr_regions = list(...)
    gold = reduce(do.call(intersect,chr_regions))
    gold = gold[width(gold) > K]
    return(gold)
  },regions[[1]],regions[[2]],
  MoreArgs = list(K = readLength),SIMPLIFY=FALSE,mc.cores = mcores)
}else if(length(regions) == 3){
  gold_regions = mcmapply(function(...,K){
    chr_regions = list(...)
    gold = reduce(do.call(intersect,chr_regions))
    gold = gold[width(gold) > K]
    return(gold)
  },regions[[1]],regions[[2]],regions[[3]],
  MoreArgs = list(K = readLength),SIMPLIFY=FALSE,mc.cores = mcores)
}

chr = names(gold_regions)


message("Calculating overlaps between reads and regions")
overlapsF = lapply(readsF,function(x,chr,gold_regions,mcores){
  out = mcmapply(reads_overlaps,chr,gold_regions,x,SIMPLIFY=FALSE,mc.cores = mcores)
  return(out)},chr,gold_regions,mcores)
overlapsR = lapply(readsR,function(x,chr,gold_regions,mcores){
  out = mcmapply(reads_overlaps,chr,gold_regions,x,SIMPLIFY=FALSE,mc.cores = mcores)
  return(out)},chr,gold_regions,mcores)

message("Extracting reads per region")
extract_readsF = mapply(function(samp_reads,samp_overlaps,chr,gold_regions,mcores){
  out = mcmapply(extract_reads,chr,samp_reads,samp_overlaps,gold_regions,
    SIMPLIFY=FALSE,mc.cores = mcores)
  return(out)}, readsF,overlapsF,MoreArgs = list(chr,gold_regions,mcores),SIMPLIFY=FALSE)
extract_readsR = mapply(function(samp_reads,samp_overlaps,chr,gold_regions,mcores){
  out = mcmapply(extract_reads,chr,samp_reads,samp_overlaps,gold_regions,
    SIMPLIFY=FALSE,mc.cores = mcores)
  return(out)}, readsR,overlapsR,MoreArgs = list(chr,gold_regions,mcores),SIMPLIFY=FALSE)

gold_nregions = mclapply(gold_regions,length,mc.cores = mcores)

message("Calculating base statistics")
depthF = lapply(extract_readsF,function(x,mcores)depth_by_region(x,mcores),mcores)
depthR = lapply(extract_readsR,function(x,mcores)depth_by_region(x,mcores),mcores)
nposF = lapply(extract_readsF,function(x,mcores)npos_by_region(x,"+",mcores),mcores)
nposR = lapply(extract_readsR,function(x,mcores)npos_by_region(x,"-",mcores),mcores)

widths = mclapply(gold_regions,function(region){
    data.table(region = 1:length(region),V1 = width(region))},
    mc.cores = mcores)

message("Forming summary statistics tables")
base = mclapply(gold_nregions,.create_summary,mc.cores =mcores)
base = mcmapply(.add_variable,base,widths,MoreArgs = list(name="width"),
  SIMPLIFY=FALSE,mc.cores =mcores)
names = c("fwd_depth","bwd_depth","fwd_npos","bwd_npos")

summaries <- mapply(function(...,base,names,mcores){
  variables = list(...)
  names(variables) = names
  out = base
  for(i in names){
    out = mcmapply(.add_variable,out,variables[[i]],MoreArgs = list(name=i),
      SIMPLIFY=FALSE,mc.cores =mcores)
  }
  return(out)
},depthF,depthR,nposF,nposR,
  MoreArgs = list(base=base,names=names,mcores=mcores),SIMPLIFY=FALSE)

                 
summaries = lapply(summaries,function(x,mcores){
  out = mclapply(x,summary_statistics,mc.cores=mcores)
  return(out)},mcores)
names(summaries) = familynames

## stack
stack_summaries <- function(summary)
{
 if(is.null(names(summary))){
   names = 1:length(summary)
 }else{
   names = names(summary)
 }
 summaries = mapply(function(x,nm){
   x[,sample:=nm]
   return(x)},summary,names,SIMPLIFY=FALSE)
 return(do.call(rbind,summaries))  
}

summaries = lapply(summaries,function(x)do.call(rbind,x))
gold_set = stack_summaries(summaries)

rf = colorRampPalette(rev(brewer.pal(11,"Spectral")))
r = rf(16)

lowerBounds = c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,100,150)


plots = mclapply(lowerBounds,function(x,gold_set){
  p  = ggplot(gold_set[npos > x],aes(dw_ratio,pbc))+stat_binhex(bins=50)+
    scale_fill_gradientn(colours =r,trans = 'log10')+facet_grid(.~sample)+
    scale_x_continuous(limits = c(-.1,10))+theme(legend.position ="top")+
    xlab("average read coverage")+ylab("unique read coverage rate")+
    ggtitle(paste0("npos > ",x))
  return(p)},gold_set,mc.cores = mcores)

pdf(file = file.path(indir,paste0(familycodename,"_gold_std_comp.pdf")),height = 5.5,
    width = 4*length(familynames))
u = lapply(plots,print)
dev.off()


plots = mclapply(lowerBounds,function(x,gold_set){
  p  = ggplot(gold_set[npos > x],aes(dw_ratio,pbc,colour = prob))+
    geom_point(alpha = I(9/10),shape = 20,size=2)+
    facet_grid(.~sample)+
    scale_color_gradient2(midpoint = .5,mid = "navyblue",low = "goldenrod",
      high = "goldenrod",guide= "colourbar",breaks = seq(-.25,1.25,by=.25))+
    scale_x_continuous(limits = c(-.1,10))+theme(legend.position ="top")+
    xlab("average read coverage")+ylab("unique read coverage rate")+
    ggtitle(paste0("npos > ",x))
  return(p)},gold_set,mc.cores = mcores)

pdf(file = file.path(indir,paste0(familycodename,"_gold_std_comp_fsr.pdf")),height = 5.5,
    width = 4*length(familynames))
u = lapply(plots,print)
dev.off()
