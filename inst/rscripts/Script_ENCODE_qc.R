#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_ENCODE_qc.R - Given a directory and a filecodename this
    script calculates both the strand cross-correlation curve, the
    normalized strand cross correlation score and the PCR botteneck
    coefficient.

    And it saves everythin in

        {indir}/{filecodename}/data


  Arguments:

    -- indir

      This is the directory where all the output files are going to be saved

    -- filecodename

      This is the code name used to save all the files that correspond to this run

    -- mcores

      This is the number of cores to be used with parallel programming
 
  Author:

    Rene Welch, Department of Statistics, University of Wisconsin - Madison

");q()}

stopifnot(length(args)==3)

indir = args[1]
filecodename = args[2]
mcores = args[3]

codedir = "/p/keles/ChIPexo/volume4/ChIPexoQual"
library(devtools)
library(grid)
library(data.table)
library(S4Vectors)
library(parallel)
library(ggplot2)
library(GenomicAlignments)

load_all(codedir)

##  indir = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human"
##  filecodename = "ER-rep1"
##  mcores = 24

indir = file.path(indir,filecodename)

message("Loading directory: ",indir)
message("Code name: ",filecodename)
message("Number of cores: ",mcores)

datadir = file.path(indir,"data")
figsdir = file.path(indir,"figs")
message("Loading statistics, regions and reads")
load(file.path(datadir,paste0(filecodename,"_regions.RData")))
load(file.path(datadir,paste0(filecodename,"_reads_by_region.RData")))


# Calculate pcr bottleneck coefficient
message("Calculating PBC")
pcr_coeff= pbc(reads_table,mcores)


# Calculate nsc
message("Calculating strand cross-correlation")
weights = do.call(c,mclapply(reads_table,nrow,mc.cores = mcores))
weights = weights / sum(weights)


fwd_cover = mclapply(reads_table,function(x){
  coverage(reads2IRanges(x[strand == "+"]))},mc.cores = mcores)
bwd_cover = mclapply(reads_table,function(x){
  coverage(reads2IRanges(x[strand == "-"]))},mc.cores = mcores)

cross_corr_curves = mcmapply(function(fwd,bwd,w,maxShift){
  start = max(runLength(fwd)[1],runLength(bwd)[1])
  end = min(tail(cumsum(runLength(fwd)),1),
    tail(cumsum(runLength(bwd)),1))
  fwd = window(fwd,start,end)
  bwd = window(bwd,start,end)
  cc1 = shiftApply(0:maxShift,fwd,bwd,cor)
  cc2 = shiftApply(1:maxShift,bwd,fwd,cor)
  dt = rbind(data.table(shift = 0:maxShift,cross.corr = cc1),
    data.table(shift = seq(-1,-maxShift,by=-1),cross.corr = cc2))
  dt$cross.corr = w * dt$cross.corr
  return(dt)  
},fwd_cover,bwd_cover,weights,MoreArgs = list(maxShift=300),
  SIMPLIFY=FALSE,mc.cores= mcores)

cc = do.call(rbind,cross_corr_curves)
cc = cc[,sum(cross.corr),by =shift]
setnames(cc,names(cc),c("shift","cross.corr"))

ccplot = ggplot(cc,aes(shift,cross.corr))+geom_line()+scale_y_continuous(limits = c(-1,1))+geom_vline(xintercept = 0,linetype= 2)+geom_abline(slope = 0,intercept = 0,linetype = 2)


# Make plots
pdf(file = file.path(figsdir,paste0(filecodename, "_cross.corr.pdf")))
print(ccplot)
dev.off()


# Save output
save(pcr_coeff,file  = file.path(datadir,paste0(filecodename, "_PBC.RData")))
save(cc,file = file.path(datadir,paste0(filecodename,"_cross.corr.RData")))
