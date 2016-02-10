
##' @rdname reads-methods
##' @aliases reads
##' @docType method
##' @exportMethod reads
setMethod("reads",
  signature = signature(object = "ChIPexo_experiment"),
  definition = function(object)object@reads
)

##' @rdname height-methods
##' @aliases height
##' @docType method
##' @exportMethod height
setMethod("height",
  signature = signature(object = "ChIPexo_experiment"),
  definition = function(object)object@height
)         

##' @rdname depth-methods
##' @aliases depth
##' @docType method
##' @exportMethod depth
setMethod("depth",
  signature = signature(object = "ChIPexo_experiment"),
  definition = function(object)object@depth
)         
