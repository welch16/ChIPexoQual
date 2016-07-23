##' @import S4Vectors
NULL

##' @rdname beta1-methods
##' @aliases beta1
##' @docType methods
setMethod("beta1",
          signature = signature(object = "ExoData"),
          definition = function(object){
              
              object@param_dist[["beta1"]]
              
          })

##' @rdname beta2-methods
##' @aliases beta2
##' @docType methods
setMethod("beta2",
          signature = signature(object = "ExoData"),
          definition = function(object){
              
              object@param_dist[["beta2"]]
              
          })

##' @rdname param_dist-methods
##' @aliases param_dist
##' @docType methods
setMethod("param_dist",
          signature = signature(object = "ExoData"),
          definition = function(object = "ExoData"){
              
              DataFrame(beta1 = object@param_dist[["beta1"]],
                        beta2 = object@param_dist[["beta2"]])
              
          })

##' @rdname nreads-methods
##' @aliases nreads
##' @docType methods
##' @exportMethod nreads
setMethod("nreads",
          signature = signature(object = "ExoData"),
          definition = function(object){
              return(metadata(object)$nreads)
          }
)

##' @rdname .MA_DF-methods
##' @aliases .MA_DF
##' @docType methods
setMethod(".MA_DF",
          signature = signature(object = "ExoData"),
          definition = function(object){
              
              A = NULL
              
              DF = mcols(object)[,c("M","A")]
              DF = subset(DF,!is.infinite(A))
              DF
             
          })

##' @rdname .ARC_URC_DF-methods
##' @aliases .ARC_URC_DF
##' @docType methods
setMethod(".ARC_URC_DF",
          signature = signature(object = "ExoData",
                                both_strand = "logical"),
          definition = function(object,both_strand){
              
              f = NULL; r = NULL
              
              DF = mcols(object)[,c("ARC","URC","f","r")]
              if(both_strand){
                  DF = subset(DF,f > 0 & r > 0)
              }
              DF[,c("ARC","URC")]
              
          })

##' @rdname .FSR_dist_DF-methods
##' @aliases .FSR_dist_DF
##' @docType methods
setMethod(".FSR_dist_DF",
          signature = signature(object = "ExoData"),
          definition = function(object,quantiles,
                                depth_values,
                                both_strand){
              
              f = NULL; r = NULL
              
              base_DF = mcols(object)[,c("f","r","d","FSR")]

              if(both_strand){
                  base_DF = subset(base_DF,f > 0 & r > 0)
              }
              
              DF_list = lapply(depth_values,.filter_quantiles,
                               base_DF,quantiles)
              do.call(rbind,DF_list)
          })

##' @rdname .region_comp_DF-methods
##' @aliases .region_comp_DF
##' @docType methods
setMethod(".region_comp_DF",
          signature = signature(object = "ExoData"),
          definition = function(object,depth_values){
              
              base_DF = mcols(object)[,c("f","r","d")]
              DF_list = lapply(depth_values,.filter_region_comp,
                               base_DF)
              do.call(rbind,DF_list)

          }) 
