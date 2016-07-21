
##' nreads methods
##' 
##' \code{nreads} returns the number of reads in the object.
##' 
##' @param object A \code{ExoData} object.
##' 
##' @return The number of reads in the \code{ExoData} object.
##' @export
##' @docType methods
##' @rdname nreads-methods
##' @examples
##' a = 1
setGeneric("nreads",
           function(object)
               standardGeneric("nreads")
)

##' .MA_DT methods
##' 
##' \code{.MA_DT} returns a \code{data.table} with the info. to generate a MA
##' plot to analyze the strand imbalance in ChIP-exo data.
##' 
##' @param object A \code{ExoData} object.
##' 
##' @return A \code{data.table} with two columns: M and A.
##' @docType methods
##' @rdname .MA_DT-methods
##' @examples 
##' a =1 
setGeneric(".MA_DT",
            function(object)
                standardGeneric(".MA_DT")
)

##' .ARC_URC_DT methods
##' 
##' \code{.ARC_URC_DT} returns a \code{data.table} with the info. to generate
##' a ARC vs URC plot to analyze enrichment and library complexity in 
##' ChIP-exo data.
##' 
##' @param object A \code{ExoData} object.
##' @param both_strand A logical value indicating if the \code{data.table} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{data.table} with two columns: ARC and URC
##' @docType methods
##' @rdname .ARC_URC_DT-methods
##' @examples 
##' a =1 
setGeneric(".ARC_URC_DT",
           function(object,...)
               standardGeneric(".ARC_URC_DT")
)

##' beta1 methods
##' 
##' \code{beta1} returns a vector with all the estimated values of the 
##' \eqn{d_i  = \beta_1 u_i + \beta_2 w_i + \epsilon_i} models fitted by 
##' \code{ChIPexoQual}
##' 
##' @param object a \code{ExoData} object.
##' @return A numeric vector with estimated values for \eqn{\beta_1}.
##' 
##' @docType methods
##' @rdname beta1-methods
##' @examples 
##' a = 1
setGeneric("beta1",
           function(object)
               standardGeneric("beta1"))

##' beta2 methods
##' 
##' \code{beta2} returns a vector with all the estimated values of the 
##' \eqn{d_i  = \beta1 u_i + \beta2 w_i + \epsilon_i} models fitted by 
##' \code{ChIPexoQual}
##' 
##' @param object a \code{ExoData} object.
##' @return A numeric vector with estimated values for \eqn{\beta_2}.
##' 
##' @docType methods
##' @rdname beta2-methods
##' @examples
##' a = 1
setGeneric("beta2",
           function(object)
               standardGeneric("beta2"))

##' param_dist methods
##' 
##' \code{param_dist} returns a \code{data.table} with all the estimated 
##' coefficients in the \eqn{d_i  = \beta_1 u_i + \beta_2 w_i + \epsilon_i} models 
##' fitted by \code{ChIPexoQual}
##' 
##' @param object a \code{ExoData} object.
##' @return A \code{data.table} with the fitted values of \eqn{\beta_1} and 
##' \eqn{\beta_2}.
##' 
##' @docType methods
##' @rdname param_dist-methods
##' @examples 
##' a = 1
setGeneric("param_dist",
           function(object)
               standardGeneric("param_dist"))







