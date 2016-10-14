##' @import S4Vectors
NULL

##' @rdname beta1-methods
##' @aliases beta1
##' @docType methods
setMethod("beta1",
          signature = signature(object = "ExoData"),
          definition = function(object){
              
              object@paramDist[["beta1"]]
              
          })

##' @rdname beta2-methods
##' @aliases beta2
##' @docType methods
setMethod("beta2",
          signature = signature(object = "ExoData"),
          definition = function(object){
              
              object@paramDist[["beta2"]]
              
          })

##' @rdname paramDist-methods
##' @aliases paramDist
##' @docType methods
setMethod("paramDist",
          signature = signature(object = "ExoData"),
          definition = function(object = "ExoData"){
              
              DataFrame(beta1 = object@paramDist[["beta1"]],
                        beta2 = object@paramDist[["beta2"]])
              
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

