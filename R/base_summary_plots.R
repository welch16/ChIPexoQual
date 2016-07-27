##' @importFrom stats quantile
##' @import S4Vectors
##' @import ggplot2
##' @import GenomicRanges
##' @import scales
##' @importFrom data.table data.table
##' @importFrom viridis viridis
##' @importFrom RColorBrewer brewer.pal
##' @importFrom dplyr desc arrange
##' @import hexbin
NULL

##' Calculates the FSR distribution of all regions with depth > value
##'  
##' Calculates the FSR distribution of all regions with depth > value
##'  
##' @param value a numeric value with the depth lower bound.
##' @param DF a \code{DataFrame} with depth and FSR.
##' @param quantiles a numeric vector with the quantiles used to summarize
##' the FSR distribution.
##' 
##' @return a \code{DataFrame} with three columns depth, quantiles and FSR.
##' 
##' @rdname .filter_quantiles
##' @name .filter_quantiles
##' 
.filter_quantiles = function(value,DF,quantiles)
{
    d = NULL
    dist = S4Vectors::subset(DF,d > value)
    quant = quantile(dist$FSR,probs = quantiles)
    names(quant) = NULL
    dist = DataFrame(d = value ,quantiles , FSR = quant)
    dist
}    

##' Calculates the region composition probabilities for all 
##' regions with depth > value
##'  
##' Calculates the region composition probabilities for all 
##' regions with depth > value
##'  
##' @param value a numeric value with the depth lower bound.
##' @param DF a \code{DataFrame} with depth and FSR.
##' 
##' @return a \code{DataFrame} with three columns d, label and prob.
##' 
##' @rdname .filter_region_comp
##' @name .filter_region_comp
##' 
.filter_region_comp = function(value,DF)
{
    d = NULL
    DF = S4Vectors::subset(DF,d > value)
    lab = c("both","fwd","bwd")
    tabl = table(fwd = factor(DF$f > 0,levels = c(TRUE,FALSE)) ,
                 bwd = factor(DF$r > 0,levels = c(TRUE,FALSE)))
    counts = c(tabl[1,1],tabl[1,2],tabl[2,1])
    props = counts / sum(counts)
    DataFrame(d = value,lab,prob = props )
}

.generate_names = function(names_input,nms, length){

    if(is.null(names_input)){
        if(is.null(nms)){
            nms = paste0("Sample: ",seq_len(length))
        }
    }else{
        nms = names_input
    }
    nms
}

.name_and_join = function(my_list,nms){
    
    DF = do.call(rbind,my_list)
    nms = mapply(rep,nms,each = vapply(my_list,nrow,1L),SIMPLIFY = FALSE)
    nms = do.call(c,nms)
    DF$sample = nms
    data.table(as.data.frame(DF))
    
}

##' MA_plot 
##' 
##' \code{MA_plot} returns a \code{ggplot} object with the MA plot to 
##' analyze the strand imbalance in ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names_input a character vector with the names to use in the
##' plot. If it is empty \code{ARC_URC_plot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' 
##' @return A \code{ggplot2} object with the MA plot.
##' @rdname MA_plot
##' @name MA_plot
##' @export
##' @examples 
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' MA_plot(exampleExoData[[1]])
MA_plot = function(...,names_input = NULL)
{
    M = NULL; A = NULL; .x = NULL
    args = unlist(list(...))
    if(!is.null(names_input)){
        stopifnot(length(names_input) == length(args))
    }
    MA_list = lapply(args,.MA_DF)
    nsamples = length(MA_list)
    nms = .generate_names(names_input,names(MA_list), 
                          nsamples)
    names(MA_list) = NULL
    MA_DF = .name_and_join(MA_list,nms)
    r = viridis::viridis(1e3, option = "D")
    p = ggplot(MA_DF,aes(M,A))+stat_binhex(bins = 70)+
        scale_fill_gradientn(colours = r,trans = 'log10',
            labels=trans_format('log10',math_format(10^.x)) )+
        theme_bw()+theme(legend.position = "top")
    if(nsamples <= 3){
        p = p + facet_wrap( ~ sample,nrow = 1)
    }else{
        p = p + facet_wrap( ~ sample,ncol = 4)
    }
    p
}

##' ARC_URC_plot 
##' 
##' \code{ARC_URC plot} returns a \code{ggplot} object with the ARC vs
##' URC plot to analyze enrichment and library complexity in ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names_input a character vector with the names to use in the
##' plot. If it is empty \code{ARC_URC_plot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' @param both_strand A logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{ggplot2} object with the ARC vs URC plot.
##' @rdname ARC_URC_plot
##' @name ARC_URC_plot
##' @export
##' @examples 
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' ARC_URC_plot(exampleExoData[[1]])
ARC_URC_plot = function(...,names_input = NULL,both_strand = FALSE)
{
    ARC = NULL; URC = NULL; .x = NULL
    
    args = unlist(list(...))
    if(!is.null(names_input)){
        stopifnot(length(names_input) == length(args))
    }
    ARC_URC_list = lapply(args,.ARC_URC_DF,both_strand)
    nsamples = length(ARC_URC_list)
    nms = .generate_names(names_input,names(ARC_URC_list), 
                          nsamples)
    names(ARC_URC_list) = NULL
    ARC_URC_DF = .name_and_join(ARC_URC_list,nms)
    r = viridis(1e3, option = "D")
    p = ggplot(ARC_URC_DF,aes(ARC,URC))+stat_binhex(bins = 50)+
        scale_fill_gradientn(colours = r,trans = 'log10',
                             labels=trans_format('log10',math_format(10^.x)) )+
        theme_bw()+theme(legend.position = "top")
    if(nsamples <= 3){
        p = p + facet_wrap( ~ sample,nrow = 1)
    }else{
        p = p + facet_wrap( ~ sample,ncol = 4)
    }
    p = p + xlim(0,3)+ylim(0,1)
    p
}

##' FSR_dist_plot 
##' 
##' \code{FSR_dist_plot} returns a \code{ggplot} object with the Forward 
##' Strand Ratio distribution plot to analyze strand imbalance in 
##' ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names_input a character vector with the names to use in the
##' plot. If it is empty \code{FSR_dist_plot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' @param quantiles a numeric vector with the quantiles used to estimate the
##' FSR distribution at a given depth. The default value is 
##' \code{c(0,.25,.5,.75,1))}
##' @param depth_values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' @param both_strand a logical value indicating if the \code{DataFrame} 
##' contains only regions with reads aligned to both strand or all. The default
##' value is \code{FALSE}.
##' 
##' @return A \code{ggplot2} object with the FSR distribution plot.
##' @rdname FSR_dist_plot
##' @name FSR_dist_plot
##' @export
##' @examples 
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' FSR_dist_plot(exampleExoData[[1]],exampleExoData[[2]])
FSR_dist_plot = function(...,names_input = NULL,
                         quantiles = c(0,.25,.5,.75,1),
                         depth_values = seq_len(30),
                         both_strand = FALSE)
{
    d = NULL; FSR = NULL
    
    args = unlist(list(...))
    if(!is.null(names_input)){
        stopifnot(length(names_input) == length(args))
    }
    FSR_list = lapply(args,.FSR_dist_DF,quantiles,
                      depth_values,
                      both_strand)
    nsamples = length(FSR_list)
    nms = .generate_names(names_input,names(FSR_list), 
                          nsamples)
    names(FSR_list) = NULL
    FSR_DF = .name_and_join(FSR_list,nms)
    p <- ggplot(FSR_DF,aes(d,FSR,colour = as.factor(quantiles)))+
        geom_line(size = 1)+
        theme_bw()+theme(legend.position = "top")+facet_grid(sample ~ .)+
        scale_color_brewer(palette = "Dark2",name = "")+
        xlab("Minimum number of reads")+ylab("Forward Strand Ratio \n (FSR)")+
        ylim(0,1)
    p
}

##' region_comp_plot 
##' 
##' \code{region_comp_plot} returns a \code{ggplot} object with the 
##' Region Composition plot to analyze strand imbalance in ChIP-exo data.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param names_input a character vector with the names to use in the
##' plot. If it is empty \code{region_comp_plot} is going to create the names
##' as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' @param depth_values a numeric vector indicating the regions with depth 
##' less or equal to, that are going to be filtered out. The defaulta values 
##' are \code{seq_len(50)}.
##' 
##' @return A \code{ggplot2} object with the Region Composition plot.
##' @rdname region_comp_plot
##' @name region_comp_plot
##' @export
##' @examples 
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' region_comp_plot(exampleExoData[[1]],exampleExoData[[2]])
region_comp_plot = function(...,names_input = NULL,
                            depth_values = seq_len(15))
{
    lab = NULL; d = NULL; prob = NULL
    
    args = unlist(list(...))
    if(!is.null(names_input)){
        stopifnot(length(names_input) == length(args))
    }
    region_list = lapply(args,.region_comp_DF,depth_values)
    nsamples = length(region_list)
    nms = .generate_names(names_input,names(region_list), 
                          nsamples)
    names(region_list) = NULL
    region_DF = .name_and_join(region_list,nms)
    r <- brewer.pal(name = "Set1",3)
    region_DF = region_DF[,lab := 
        factor(lab,levels = c("both","fwd","bwd"))]
    p <- ggplot(arrange(region_DF,lab, desc(prob)),
                  aes(d,prob,fill = lab))+
        geom_bar(stat = "identity")+
        theme_bw()+theme(legend.position = "top")+
        facet_grid(sample ~ .)+
        scale_fill_brewer(palette = "Pastel1",name = "Strand composition")+
        xlab("Minimum number of reads")+ylab("Proportion of regions")
    p
    
}

##' param_dist_boxplot 
##' 
##' \code{param_dist_boxplot} returns a \code{ggplot} object with a 
##' boxplot comparing the \code{ntimes} estimations of the chosen 
##' parameter.
##' 
##' @param ... a \code{list} of \code{ExoData} objects, or several 
##' \code{ExoData} objects by themselves.
##' @param which_param a character value with either \code{"beta1"} or 
##' \code{"beta2"} that determines which paramters in the model 
##' d_i ~ u_i + w_i to plot. The default value is \code{"beta1"}.
##' @param names_input a character vector with the names to use in the
##' plot. If it is empty \code{param_dist_boxplot} is going to create the 
##' names as the names of the list when they are available or is going to 
##' name them as Sample: 1 ,... , Sample: k.
##' 
##' @return A \code{ggplot2} object with the boxplot of the chosen 
##' parameter
##' @rdname param_dist_boxplot
##' @name param_dist_boxplot
##' @export
##' @examples 
##' library(ChIPexoQualExample)
##' data(exampleExoData)
##' param_dist_boxplot(exampleExoData[[1]],exampleExoData[[2]])
param_dist_boxplot = function(...,names_input = NULL,
                            which_param = "beta1")
{
    lab = NULL; d = NULL; prob = NULL
    
    stopifnot(which_param %in% c("beta1","beta2"))
    args = unlist(list(...))
    if(!is.null(names_input)){
        stopifnot(length(names_input) == length(args))
    }
    param_list = lapply(args,param_dist)
    nsamples = length(param_list)
    nms = .generate_names(names_input,names(param_list), 
                          nsamples)
    names(param_list) = NULL
    param_DF = .name_and_join(param_list,nms)
    param_DF$beta2 = - param_DF$beta2
    p <- ggplot(param_DF,
                aes_string(x = "sample",
                           y = which_param))+
        geom_boxplot()+
        theme_bw()+
        theme(legend.position = "top",
              axis.title.y = element_text(angle = 0))
    if(which_param == "beta1"){
        p = p + ylab(expression(beta[1]))+
            geom_abline(slope = 0,intercept = 10,linetype = 2)
    }else{
        p = p + ylab(expression(beta[2]))+
            geom_abline(slope = 0,intercept = 0,linetype = 2)
    }
    p+xlab("Sample")
}

