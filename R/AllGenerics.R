
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
##' data(exoExample)
##' nreads(exoExample)
setGeneric("nreads",
           function(object)
               standardGeneric("nreads")
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
##' @export
##' @examples 
##' data(exoExample)
##' beta1(exoExample)
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
##' @export
##' @examples
##' data(exoExample)
##' beta2(exoExample)
setGeneric("beta2",
           function(object)
               standardGeneric("beta2"))

##' paramDist methods
##' 
##' \code{paramDist} returns a \code{DataFrame} with all the estimated 
##' coefficients in the \eqn{d_i  = \beta_1 u_i + \beta_2 w_i + \epsilon_i} models 
##' fitted by \code{ChIPexoQual}
##' 
##' @param object a \code{ExoData} object.
##' @return A \code{DataFrame} with the fitted values of \eqn{\beta_1} and 
##' \eqn{\beta_2}.
##' 
##' @docType methods
##' @rdname paramDist-methods
##' @export
##' @examples 
##' data(exoExample)
##' paramDist(exoExample)
setGeneric("paramDist",
           function(object)
               standardGeneric("paramDist"))

