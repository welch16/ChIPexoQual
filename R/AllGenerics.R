
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
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' nreads(exampleExoData[[1]])
setGeneric("nreads",
           function(object)
               standardGeneric("nreads")
)

##' .MA_DF methods
##' 
##' \code{.MA_DF} returns a \code{DataFrame} with the info. to generate a MA
##' plot to analyze the strand imbalance in ChIP-exo data.
##' 
##' @param object A \code{ExoData} object.
##' 
##' @return A \code{DataFrame} with two columns: M and A.
##' @docType methods
##' @rdname .MA_DF-methods
setGeneric(".MA_DF",
            function(object)
                standardGeneric(".MA_DF")
)

##' .FSR_dist_DF methods
##' 
##' \code{.FSR_dist_DF} return a \code{DataFrame} with the info. necessary 
##' to generate the FSR distribution plot to analyze strand imbalance in 
##' ChIP-exo data.
##' 
##' @param object a \code{ExoData} object.
##' @param quantiles a numeric vector with the quantiles used to estimate the
##' FSR distribution at a given depth. The default value is 
##' \code(c(0,.25,.5,.75,1)).
##' @param depth_values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' @param both_strand a logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' @return A \code{DataFrame} with three columns: depth, quantiles and FSR.
##' @rdname .FSR_dist_DF-methods
##' @docType methods
setGeneric(".FSR_dist_DF",
           function(object,...)
               standardGeneric(".FSR_dist_DF"))

##' .region_comp_DF methods
##' \code{.region_comp_DF} returns a \code{DataFrame} with the info. 
##' necessary to generate the Region Composition plot to analyze strand
##' imbalance in ChIP-exo data.
##' 
##' @param object a \code{ExoData} object.
##' @param depth_values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' @param both_strand a logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' @return A \code{DataFrame} with three columns: depth, label and propotion
##' of regions at a given depth level with the respective label.
##' @rdname .region_comp_DF-methods
##' @docType methods
setGeneric(".region_comp_DF",
           function(object,...)
               standardGeneric(".region_comp_DF"))


##' .ARC_URC_DF methods
##' 
##' \code{.ARC_URC_DF} returns a \code{DataFrame} with the info. to generate
##' a ARC vs URC plot to analyze enrichment and library complexity in 
##' ChIP-exo data.
##' 
##' @param object A \code{ExoData} object.
##' @param both_strand A logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{DataFrame} with two columns: ARC and URC
##' @docType methods
##' @rdname .ARC_URC_DF-methods
setGeneric(".ARC_URC_DF",
           function(object,...)
               standardGeneric(".ARC_URC_DF")
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
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' beta1(exampleExoData[[1]])
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
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' beta2(exampleExoData[[1]])
setGeneric("beta2",
           function(object)
               standardGeneric("beta2"))

##' param_dist methods
##' 
##' \code{param_dist} returns a \code{DataFrame} with all the estimated 
##' coefficients in the \eqn{d_i  = \beta_1 u_i + \beta_2 w_i + \epsilon_i} models 
##' fitted by \code{ChIPexoQual}
##' 
##' @param object a \code{ExoData} object.
##' @return A \code{DataFrame} with the fitted values of \eqn{\beta_1} and 
##' \eqn{\beta_2}.
##' 
##' @docType methods
##' @rdname param_dist-methods
##' @export
##' @examples 
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' param_dist(exampleExoData[[1]])
setGeneric("param_dist",
           function(object)
               standardGeneric("param_dist"))

