
#' @import ggplot2
NULL

.fn_filter <- function(chr_measurement,chr_depth,lower)return(chr_measurement[chr_depth > lower])

#' Returns a ggplot object with a boxplot of lowerBounds vs chr_measurement. It filters the measurement by the islands such that the depth is greater than lowerBound
#'
#' @param lowerBounds An integer vector with the lower bounds used to filter
#'
#' @param chr_depths A list of vectors with the islands depths per chromosome
#'
#' @param chr_measurement A list of vector with the islands given measurements per chromosome
#'
#' @param measurement_label A string with the ylab of the plot
#'
#' @param log A boolean variable indicating if the y-axis is going to be rescaled with log10
#'
#' @param mc Number of cores used
#'
#' @export
filter_regions_plot <- function(lowerBounds,chr_depths,chr_measurement,measurement_label,mc,log = FALSE)
{
  filter_regions = lapply(lowerBounds,function(x){
    mcmapply(.fn_filter,chr_measurement,chr_depths,MoreArgs = list(lower = x),
      SIMPLIFY = FALSE,mc.cores = mc)})           
  names(filter_regions) = lowerBounds
 
  filter_regions = mcmapply(function(x,y){
    data.table(bound = x,measurement = do.call(c,y))},lowerBounds,filter_regions,MoreArgs = list(),
    SIMPLIFY = FALSE,mc.cores = mc)

  filter_regions = do.call(rbind,filter_regions)
  filter_regions$bound = factor(filter_regions$bound)

  p = ggplot(filter_regions,aes(bound,measurement,colour=bound))+
    geom_boxplot(outlier.colour = alpha("black",1/50))+
    theme(legend.position = "none")+
    ylab(measurement_label)

  if(log){
    p = p + scale_y_log10(  labels=trans_format('log10', math_format(10^.x)))
  }
  return(p)
}


#' Returns a ggplot object with an area plot of the label distribution per lower bound
#'
#' @param lowerBounds An integer vector with the lower bounds used to filter
#'
#' @param chr_depths A list of vectors with the islands depths per chromosome
#'
#' @param labels A list of character vectors with the labels for each island
#'
#' @param mc Number of cores used
#'
#' @export

filter_label_plot <- function(lowerBounds,chr_depths,labels,mc)
{
  filter_labels= lapply(lowerBounds,function(x){
    mcmapply(.fn_filter,labels,chr_depths,MoreArgs = list(lower = x),
    SIMPLIFY = FALSE,mc.cores = mc)})           
  names(filter_labels) = lowerBounds

  filter_labels = mcmapply(function(x,y){
    data.table(bound = x,measurement = do.call(c,y))},lowerBounds,filter_labels,MoreArgs = list(),
    SIMPLIFY = FALSE,mc.cores = mc)

  filter_labels = do.call(rbind,filter_labels)
  filter_labels$bound = factor(filter_labels$bound)
  filter_labels$measurement = factor(filter_labels$measurement)
 
  p = ggplot(filter_labels,aes(bound,fill = measurement))+
    geom_bar(position = "fill")+
    scale_fill_brewer(palette = "Set1")+
    theme(legend.position = "bottom")+
    ylab("proportion by label")
  
  return(p) 
}


#' Returns a ggplot object with a collection of MA plots faceted by lower bound
#'
#' @param lowerBounds An integer vector with the lower bounds used to filter
#'
#' @param chr_depths A list of vectors with the islands depths per chromosome
#'
#' @param chr_M_values A list of vectors with the M coordinate of an MA plot
#'
#' @param chr_A_values A list of vectors with the A coordinate of an MA plot
#'
#' @param mc Number of cores used
#'
#' @export
filter_MA_plot <- function(lowerBounds,chr_depths,chr_M_values,chr_A_values,mc,smooth=FALSE)
{

  filter_M = lapply(lowerBounds,function(x){
    mcmapply(.fn_filter,chr_M_values,chr_depths,MoreArgs = list(lower = x),
      SIMPLIFY = FALSE,mc.cores = mc)})           

  filter_A = lapply(lowerBounds,function(x){
    mcmapply(.fn_filter,chr_A_values,chr_depths,MoreArgs = list(lower = x),
      SIMPLIFY = FALSE,mc.cores = mc)})           
    
  filter_regions = mcmapply(function(x,Mval,Aval){
    data.table(bound = x,A = do.call(c,Mval),M = do.call(c,Aval))},lowerBounds,filter_M,filter_A,MoreArgs = list(),
    SIMPLIFY = FALSE,mc.cores = mc)

  filter_regions = do.call(rbind,filter_regions)
  filter_regions$bound = factor(filter_regions$bound)

  filter_regions = filter_regions[ !is.infinite(M) & !is.infinite(A),]

  rf = colorRampPalette(rev(brewer.pal(11,"Spectral")))
  r = rf(16)
  
  p = ggplot(filter_regions,aes(A,M))+stat_binhex(bins = 70)+
    facet_wrap(~bound,ncol =4)+
    scale_fill_gradientn(colours =r,trans='log10')+
    geom_abline(slope=0,intercept = 0,linetype =2)+
    theme(legend.position = "top")

  if(smooth){
      p = p + geom_smooth(method = "loess",na.rm=TRUE,se=FALSE,colour = I("red"))
  }
  
  return(p)
}

