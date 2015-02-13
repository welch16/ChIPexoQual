#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    Script_calculate_cross_corr.R - Given a directory and a filecodename
    this script calculates local strand cross-correlation for the regions

    Under this setting and following the filtering of Script_plot_peaks.R, i.e.
      - summary_stats[ npos >2 ]
      - summary_stats[depth > median(depth)]
      - summary_stats[label == 'both']

    Right now, it updates the summary_statistics by adding (whenever possible) the
    following summary statistics:

    - number of local maxima
    - number of local maxima in window ( -readLength , readLength)
    - shift of global maxima
    - voting rule of local maxima in window (consider the sign of the shift by number of
      votes, if there is a tie, use global maxima's sign, zeroes goes to majority, if there
      are no elements in window then there is no sig)n

    Finally saves a different data.table with all the local maxima of a given region

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
library(grid)
library(data.table)
library(S4Vectors)
library(parallel)
library(ggplot2)
library(GenomicAlignments)

load_all(codedir)

##  indir = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human"
##  filecodename = "ER-rep1"
##  chr = "chr1"
## # readLength = 36
##  mcores = 24

indir = file.path(indir,filecodename)

message("Loading directory: ",indir)
message("Code name: ",filecodename)
message("Number of cores: ",mcores)

figsdir = file.path(indir,"figs")
datadir = file.path(indir,"data")

message("Loading summary statistics and reads")
statsfile = file.path(datadir,paste0(filecodename,"_summary_statistics.RData"))
load(statsfile)
load(file.path(datadir,paste0(filecodename,"_reads_by_region.RData")))

stats = summary_stats[[chr]]
reads = reads_table[[chr]]

readLength = floor(median(reads$end - reads$start + 1))

# Apply same filter as in Script_plot_peaks.R
stats = stats[npos > 2]
stats = stats[depth > median(depth)]
stats = stats[label == "both"]

regions = stats[,.(region)][[1]]
setkey(reads,region)
reads = reads[region == regions]

message("Calculating forward coverages")
fwd_cover = mclapply(regions,function(x){
  coverage(reads2IRanges(reads[region == x & strand == "+"]))},
  mc.cores = mcores)

message("Calculating backward coverages")
bwd_cover = mclapply(regions,function(x){
  coverage(reads2IRanges(reads[region == x & strand == "-"]))},
  mc.cores = mcores)

message("Calculating strand-cross correlation limits")
mins = mcmapply(function(X,Y)max(runLength(X)[1],runLength(Y)[1]),fwd_cover,bwd_cover,
  SIMPLIFY = FALSE,mc.cores = mcores)
maxs = mcmapply(function(X,Y)min(max(cumsum(runLength(X))),max(cumsum(runLength(Y)))),
  fwd_cover,bwd_cover,SIMPLIFY = FALSE,mc.cores = mcores)
fwd_cover= mcmapply(function(cover,m,M)window(cover,start = m,end = M),
  fwd_cover,mins,maxs,SIMPLIFY = FALSE,mc.cores = mcores)
bwd_cover = mcmapply(function(cover,m,M)window(cover,start = m,end = M),
  bwd_cover,mins,maxs,SIMPLIFY = FALSE,mc.cores = mcores)

widths = stats[,.(width)][[1]]

cross_corr_data.table <- function(X,Y,width,rr,mcores)
{
  message(rr)
  shift = c(0,seq_len(.5*width))
  shift = as.integer(shift[shift < length(X)])
  cc1 = mcshiftapply(shift,cor,X,Y,mcores)
  cc2 = mcshiftapply(shift,cor,Y,X,mcores)

  dt = rbind(
    data.table(shift=shift,cross.cor = cc1),
    data.table(shift=-shift[-1],cross.cor=cc2[-1]))
  dt = dt[!is.nan(cross.cor)]
  dt[,region:= rr]
  
  return(dt)
}

message("Calculating cross-corr curves")
cross.corr =  mapply(cross_corr_data.table,
  fwd_cover,bwd_cover,widths,regions,MoreArgs = list(mcores),SIMPLIFY= FALSE)#,mc.cores = mcores)
            
local_maxima_data.table <- function(dt,readLength = NULL)
{
  if(nrow(dt)>0){
    maxshift = max(abs(dt$shift))
    if(nrow(dt) > 1){
      dt = dt[localMaxima(cross.cor)]
    }      
    dt = dt[cross.cor > 0]
    dt = dt[abs(shift) < maxshift]
    dt[,bounds:="no"]
    if(!is.null(readLength)){
      dt[data.table::between(shift,-readLength-1,readLength+1),bounds:="yes"]
    }
  }else{
    dt[,bounds:= character()]
  }
  return(dt)
}

message("Extract local non-negative local maxima")
local.max = mclapply(cross.corr,local_maxima_data.table,readLength,mc.cores = mcores)

# Calculate new summary statistics
message("Calculating cross.corr summary statistics")
new_stats = data.table(region = stats[,(region)])
new_stats[,nr_loc_max := do.call(c,mclapply(local.max,nrow,mc.cores = mcores))]
new_stats[,nr_loc_max_wind := do.call(c,
  mclapply(local.max,function(x)nrow(x[bounds=="yes"]),mc.cores = mcores))]
new_stats[,global_shift:=do.call(c,
  mclapply(local.max,function(x){
    if(nrow(x) > 0){
      out = x[which.max(cross.cor),(shift)]
    }else{
      out = NA
    }
    return(out)
  },mc.cores = mcores))
]
new_stats[,digestion:= do.call(c,
  mclapply(local.max,function(x){
    out = NA
    if(nrow(x) > 0){
      x = x[bounds == "yes"]
      if(nrow(x) > 0){
        zeros = nrow(x[shift == 0])
        pos = nrow(x[shift > 0])
        neg = nrow(x[shift < 0])
        if(pos == neg){
          out = sign(x[which.max(cross.cor),(shift)])
        }else{
          out = 1
          if(pos < neg)out = -1
        }
      }
    }
    return(out)
  },mc.cores = mcores))]
setkey(new_stats,region)

## new_stats = data.table(left_join(summary_stats[[chr]],new_stats,by="region"))
## setkey(new_stats,region)

cross_corr_plot <- function(cross.cor,loc.max,readLength=NULL)
{
  p = ggplot(cross.cor,aes(shift,cross.cor))+geom_line()+
    geom_abline(slope=0,intercept=0,linetype=2)
    if(!is.null(readLength)){
      p = p +
        geom_vline(xintercept = loc.max[bounds=="no",(shift)],colour = "red",linetype = 2)+
        geom_vline(xintercept = loc.max[bounds=="yes",(shift)],colour = "blue",linetype = 2)          
    }else{
      p = p + geom_vline(xintercept = loc.max[,(shift)],colour = "red",linetype=2)    
    }
  annot = loc.max[,(shift)]
  annot = paste0(annot,",")
  annot = do.call(paste0,lapply(annot,function(x)x))
  p = p + ggtitle(annot)+scale_y_continuous(limits = c(-1,1))
  
  return(p)
}

cc_plots = mcmapply(cross_corr_plot,cross.corr,local.max,readLength,SIMPLIFY = FALSE,
  mc.cores = mcores)

# This sorts the peaks indexed with the decreasing npos (unique number of positions)
idd = sort(stats$npos,decreasing =TRUE,index.return=TRUE)

message("Plotting peaks")
pdf(file = file.path(figsdir,paste0(filecodename,"-",chr,"_cross_corr.pdf")))
u = lapply(idd$ix,function(j)print(cc_plots[[j]]))
dev.off()


statsfile = file.path(datadir,paste0(filecodename,"-",chr,"_cross_corr_statistics.RData"))
cross_corr_stats = new_stats
save(cross_corr_stats,file = statsfile)
cross_corr = do.call(rbind,cross.corr)
save(cross_corr,file = file.path(datadir,paste0(filecodename,"-",chr,"_cross_corr.RData")))
                  
# criteria for local maxima of strand cross correlation
#
# 1 - The local maxima must be a global
# 2 - The position where the local maxima must be close to zero
# 3 - The local maxima must be non-negative
# 4 - Can remove local maximae that are close to each other. in that case, we can change the collection of close local maximas by their average
# 5 - There is a pattern, where there is some simmetry with the local maximae, that may suggests that there are more than one binding events in  the region, which may suggest:
#   a - nmaxima may be an indicator of > 1 possible binding events in region
#   b - this suggest the need to improve the partition in that region (may be an overkill)
# 6 - Must remove a local maxima which is on an extreme point of the shift
# 7 - When the SHIFT vector is very long, there may be cases where:
#   a - The behaviour at the extremes it looks kind off noisy
#   b - Several local maximae in those cases
#   c - This may be avoided if we consider a max. shift of a * width(region) where
#       0 <= a <= 1, right now we are using a = .75 but may be worth checking for smaller a's


# there are two effects to take into account:
# 1 - as the depth (or npos) of the region decreases, then the cross - correlation function
#   becomes noisier, therefore it has more local - maxima and therefore it may be harder to
#   find the correct shift. Under this case, the local maxima may be closer to each other
# 2 - when there are more than one binding events, there are more local maxima
