##' @importFrom S4Vectors subjectHits
NULL

##' ExoDataBlacklist
##'
##' \code{ExoDataBlacklist} separates the regions in an \code{ExoData} object by overlapping them
##' with a set of blacklisted regions and calculates the quality parameters in both collections of 
##' islands.
##' 
##' @param exo a \code{ExoData} object.
##' @param blacklist a \code{GRanges} object with the blacklisted regions or a \code{character} indicating 
##' which of the blacklist included in \code{ChIPexoQual} to use.
##' @param which.param a character value with either \code{"beta1"} or 
##' \code{"beta2"} that determines which parameters in the model 
##' depth_i ~ uniquePos_i + width_i to plot. The default value is 
##' \code{"beta1"}.
##' @param nregions a numeric value indicating the number of regions sampled to 
##' estimate the quality parameter distributions. The default value is extracted from \code{exo}.
##' @param ntimes a numeric value indicating the number of times that regions are 
##' sampled to estimate the quality parameter distributions. The default value
##' is extracted from \code{object}.
##' @return A \code{ggplot} object with a boxplot that compares the quality scores distribution when the regions
##' overlap a pre-defined collection of blacklists.
##' 
##' @examples 
##' data(exoExample)
##' data(blacklists)
##' ExoDataBlacklist(exoExample,blacklists[["mm9"]],ntimes = 10,nregions = 500)
##'
##'
##' @rdname ExoDataBlacklist
##' @export
ExoDataBlacklist <- function(exo, blacklist , which.param = "beta1",
                             nregions = NULL,ntimes = NULL)
{
  term <- overlap <- estimate <- depth <- uniquePos <- NULL
  
  stopifnot(class(exo) == "ExoData")
  stopifnot(class(blacklist) %in% c("GRanges"))
  stopifnot(which.param %in% c("beta1","beta2"))

  if(is.null(nregions)){
    nregions <- metadata(exo)$nregions
  }
  
  if(is.null(ntimes)){
    ntimes <- metadata(exo)$ntimeskin
  }
  
  stopifnot(length(exo) > nregions)
  
  overlaps <- findOverlaps(exo,blacklist)
  
  inBlacklist <- formatRegions(exo[subjectHits(overlaps)])
  outBlacklist <- formatRegions(exo[-subjectHits(overlaps)])
  
  inBlacklistParam <- lapply(seq_len(ntimes),calculateParamDist,inBlacklist,nregions)
  outBlacklistParam <- lapply(seq_len(ntimes),calculateParamDist,outBlacklist,nregions)
  allParam <- lapply(seq_len(ntimes),calculateParamDist,formatRegions(exo),nregions)
 
  inBlacklistParam <- rbindlist(inBlacklistParam)
  outBlacklistParam <- rbindlist(outBlacklistParam)
  allParam <- rbindlist(allParam)
  
  
  if(which.param == "beta1"){
    DT <- rbindlist(list(
      allParam[term == "uniquePos"][,overlap := "All regions"],
      inBlacklistParam[term == "uniquePos"][,overlap := "Overlap blacklists"],
      outBlacklistParam[term == "uniquePos"][,overlap := "Don't overlap blacklists"]
    ))
  }else{
    DT <- rbindlist(list(
      allParam[term != "uniquePos"][,overlap := "All regions"],
      inBlacklistParam[term != "uniquePos"][,overlap := "Overlap blacklists"],
      outBlacklistParam[term != "uniquePos"][,overlap := "Don't overlap blacklists"]
    ))
    DT <- DT[,estimate := - estimate]
  }
  out <- ggplot(DT,aes_string("overlap","estimate",colour = "overlap"))+geom_boxplot()+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          legend.position = "top"
          )+
    scale_color_brewer(palette = "Set1",name = "")
 
  if(which.param == "beta1"){
    out <- out + ylab(expression(beta[1]))
  }else{
    out <- out + ylab(expression(beta[2]))
  }
  out

}


# ##' availableBlacklists
# ##'
# ##' \code{availableBlacklists} returns the blacklists that are included in the \code{ChIPexoQual} package.
# ##' 
# ##' @return A \code{character} vectors with the blacklists included in \code{ChIPexoQual}
# ##'
# ##' @examples
# ##' availableBlacklists()
# ##' 
# ##' @rdname availableBlacklists
# ##' @export
# availableBlacklists <- function(){
#   blacklists <- NULL
#   data(blacklists)
#   names(blacklists)
# }
