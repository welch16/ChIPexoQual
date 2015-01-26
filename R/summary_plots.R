
#' @import ggplot2
NULL

.fn_filter <- function(lower,summary_table,byVar)
{
  setkeyv(summary_table,cols = byVar)
  return(summary_table[ eval(parse(text=byVar)) > lower,])
}

#' Returns a ggplot object with a boxplot of lowerBounds vs chr_measurement. It filters the measurement by the islands such that the depth is greater than lowerBound
#'
#' @param lowerBounds An integer vector with the lower bounds used to filter
#'
#' @param filtered_summary A list of data.tables with the already filtered summary statistics
#'
#' @param measure_var A string with the name of the variable to be drawn
#'
#' @param measure_label A string with the y-axis label
#'
#' @param log A boolean variable indicating if the y-axis is going to be rescaled with log10
#'
#' @param mc Number of cores used
#'
#' @export
filter_regions_plot <- function(lowerBounds,filtered_summary,
  measure_var,measure_label,mc,log=FALSE)
{

  plot_data = mcmapply(function(lower,filtered){
    dt = data.table(bound=lower,
      measurement= filtered[[measure_var]])
    return(dt)},lowerBounds,filtered_summary,MoreArgs = list(),
  SIMPLIFY=FALSE,mc.cores = mc)
  
  plot_data = do.call(rbind,plot_data)
  plot_data$bound = factor(plot_data$bound)
  
  p = ggplot(plot_data,aes(bound,measurement,colour=bound))+
    geom_boxplot(outlier.colour = alpha("black",1/50))+
    theme(legend.position = "none")+
    ylab(measure_label)

  if(log){
    p = p + scale_y_log10(  labels=trans_format('log10',
                              math_format(10^.x)))
  }
  return(p)
}



#' Returns a ggplot object with an area plot of the label distribution per lower bound. It filters the labels by the islands such that the depth is greater than lowerBound
#'
#' @param lowerBounds An integer vector with the lower bounds used to filter
#'
#' @param filtered_summary A list of data.tables with the already filtered summary statistics
#'
#' @param mc Number of cores used
#'
#' @export

filter_label_plot <- function(lowerBounds,filtered_summary,mc)
{
  plot_data = mcmapply(function(lower,filtered){
    dt = data.table(bound=lower,
      measurement= filtered[["label"]])
    return(dt)},lowerBounds,filtered_summary,MoreArgs = list(),
  SIMPLIFY=FALSE,mc.cores = mc)
  
  plot_data = do.call(rbind,plot_data)
  plot_data$bound = factor(plot_data$bound)
  
  p = ggplot(plot_data,aes(bound,fill = measurement))+
    geom_bar(position = "fill")+
    scale_fill_brewer(palette = "Set1")+
    theme(legend.position = "bottom")+
    ylab("proportion by label")
  
  return(p) 
}


#' Returns a ggplot object with a collection of scatter plots faceted by lower bound. It filters the labels by the islands such that the depth is greater than lowerBound
#'
#' @param lowerBounds An integer vector with the lower bounds used to filter
#'
#' @param var1 String with the name of which variable is going to be plotted in the x-axis
#'
#' @param var2 String with the name of which variable is going to be plotted in the y-axis
#'
#' @param filtered_summary A list of data.tables with the already filtered summary statistics
#'
#' @param mc Number of cores used
#'
#' @param smooth A boolean variable indicating if a loess regression line is going to be added to the plot
#'
#' @export

filter_scatter_plot <- function(lowerBounds,var1,var2,filtered_summary,mc,smooth=FALSE)
{
  plot_data = mcmapply(function(lower,filtered){
    dt = data.table(bound=lower,
      A= filtered[[var1]],M = filtered[[var2]])
    return(dt)},lowerBounds,filtered_summary,MoreArgs = list(),
  SIMPLIFY=FALSE,mc.cores = mc)
  
  plot_data = do.call(rbind,plot_data)
  plot_data$bound = factor(plot_data$bound)

  rf = colorRampPalette(rev(brewer.pal(11,"Spectral")))
  r = rf(16)

  p = ggplot(plot_data[!is.infinite(A) & !is.infinite(M)],aes(A,M))+stat_binhex(bins = 70)+
    facet_wrap(~bound,ncol =4)+
    scale_fill_gradientn(colours =r,trans='log10')+
    theme(legend.position = "top")

  if(smooth){
      p = p + geom_smooth(method = "loess",na.rm=TRUE,se=FALSE,colour = I("red"))
  }
  
  return(p)
}



#' Returns a ggplot object with a map of the number of positions vs the number of reads, where the color denotes how many regions agree with that specific combination
#'
#' @param summary_table A data.table with the summary statistics per region containing at least the npos and depth fields
#'
#' @param mp Max number of positions
#'
#' @param md Max depth
#'
#' @export
positions_reads_map <- function(summary_table,mp=Inf,md=Inf)
{
  df = summary_table
  setkey(df,depth,npos)
  mat = df[,length(chrID),by=list(depth,npos)]
  setnames(mat,names(mat),c(names(mat)[1:2],"nr_islands"))
  rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r = rf(16)
  max_pos = mp
  max_depth = md
  maxdepth = max(mat[npos <= max_pos & depth <= max_depth,list(depth)])
  p = ggplot(mat[npos <= max_pos],aes(as.factor(npos),as.factor(depth),fill = nr_islands))+
   geom_tile()+scale_fill_gradientn(name = "number of islands",colours = r,trans="log10",labels = trans_format('log10',math_format(10^.x)))+
   scale_y_discrete(breaks = c(seq(0,maxdepth,by=20),maxdepth) )+
   theme(legend.position = "top",axis.text.x = element_text(angle = 90))+xlab("number of positions per island")+
    ylab("number of reads per island")
  return(p)
}


