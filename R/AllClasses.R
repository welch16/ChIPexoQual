##' ChIP-exo experiment class description.
##' 
##' Contains the reads of a ChIP exo experiment.
##' 
##' @slot file Character vector wih the name of the bam file of a ChIP-exo experiment.
##' @slot reads A \code{GRanges} object with the reads of the experiment.
##' @slot depth A numeric value with the experiment's depth.
##' @slot height Tunning parameter to partition the experiment into regions.
##' @slot summary_stats A \code{data.table} with a collection of summary statistics of the ChIP-exo experiment.
##' 
##' @name ChIPexo_experiment-class
##' @rdname ChIPexo_experiment-class
##' @exportClass ChIPexo_experiment
setClass("ChIPexo_experiment",
    representation = representation(
      file = "character",
      reads = "GRanges",
      depth = "numeric",
      height = "numeric",
      summary_stats = "data.table"
    ),
    prototype = prototype(
      file = "",
      reads = GRanges(),
      depth = 0,
      height = 1,
      summary_stats = data.table()
    )
)