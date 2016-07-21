##' @importFrom data.table as.data.table
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
              
              as.data.table(object@param_dist)
              
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

##' @rdname .MA_DT-methods
##' @aliases .MA_DT
##' @docType methods
setMethod(".MA_DT",
          signature = signature(object = "ExoData"),
          definition = function(object){
              
              DT = data.table(M = object$M,A = object$A)
              DT = DT[!is.infinite(A)]
              DT
             
          })

##' @rdname .ARC_URC_DT-methods
##' @aliases .ARC_URC_DT
##' @docType methods
setMethod(".ARC_URC_DT",
          signature = signature(object = "ExoData"),
          definition = function(object,both_strand = FALSE){
              
              DT = data.table(ARC = object$ARC , URC = object$URC)
              if(both_strand){
                  DT = DT[ as.character(object$label) == "both"]
              }
              DT
              
          })

##' @rdname .FSR_dist-methods
##' @aliases .FSR_dist
##' @docType methods
