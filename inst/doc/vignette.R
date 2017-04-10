## ----style, echo = FALSE, results = 'asis'-------------------------------

    library(BiocStyle)
    markdown(css.files = c('custom.css'))


## ----extraload,include = FALSE,echo = FALSE,eval = TRUE------------------
    
    library(dplyr,quietly = TRUE)
    library(BiocParallel)
    library(gridExtra)
    library(ChIPexoQual)
    library(ChIPexoQualExample)
    library(knitr)
    opts_chunk$set(fig.align = "center")


## ----load,include=TRUE,echo=TRUE,eval=FALSE------------------------------
#  
#      library(ChIPexoQual)
#      library(ChIPexoQualExample)
#  

## ----ex1,include=TRUE,echo=TRUE,eval=TRUE---------------------------------------------------
    
    files = list.files(system.file("extdata",
        package = "ChIPexoQualExample"),full.names = TRUE)
    basename(files[1])
    ex1 = ExoData(file = files[1],mc.cores = 2L,verbose = FALSE)
    ex1
    reads = readGAlignments(files[1],param = NULL)
    ex2 = ExoData(reads = reads,mc.cores = 2L,verbose = FALSE)
    identical(GRanges(ex1),GRanges(ex2))


## ----pre_create,include=FALSE,echo=FALSE,eval = TRUE----------------------------------------
    
    rm(reads,ex1,ex2)
    

## ----create_exdata,include=TRUE,echo=TRUE,eval=TRUE-----------------------------------------
    
    files = files[grep("bai",files,invert = TRUE)] ## ignore index files
    exampleExoData = lapply(files,ExoData,mc.cores = 2L,verbose = FALSE)


## ----depth,include=TRUE,echo=TRUE,eval=TRUE-------------------------------------------------

    sapply(exampleExoData,nreads)


## ----ARC1,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,fig.height=4,fig.width=9,results='markup'----
    
    ARCvURCplot(exampleExoData,names.input = paste("Rep",1:3,sep = "-"))
    

## ----strand1,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,results="markup",fig.height=4----

    p1 = regionCompplot(exampleExoData,names.input = paste("Rep",1:3,
        sep = "-"),depth.values = seq_len(50))
    p2 = FSRDistplot(exampleExoData,names.input = paste("Rep",1:3,sep = "-"),
        quantiles = c(.25,.5,.75),depth.values = seq_len(100))
    gridExtra::grid.arrange(p1,p2,nrow = 1)
    

## ----ARC2,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,fig.height=4,fig.width=9,results='markup'----
    
    ARCvURCplot(exampleExoData[[1]],
                 subset(exampleExoData[[1]],uniquePos > 10),
                 subset(exampleExoData[[1]],uniquePos > 20),
                 names.input = c("All", "uniquePos > 10", "uniquePos > 20"))
    

## ----param_dist,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,results="markup",fig.height=2.5----
    
    p1 = paramDistBoxplot(exampleExoData,which.param = "beta1", names.input = paste("Rep",1:3,sep = "-"))
    p2 = paramDistBoxplot(exampleExoData,which.param = "beta2", names.input = paste("Rep",1:3,sep = "-"))
    gridExtra::grid.arrange(p1,p2,nrow = 1)


## ----param_eval,include=TRUE,echo=TRUE,eval=TRUE--------------------------------------------
    
    sapply(exampleExoData,function(x)median(beta1(x)))
    sapply(exampleExoData,function(x)median(-beta2(x)))


## ----subsamping,include=TRUE,echo=TRUE,eval=TRUE--------------------------------------------

    sample.depth = seq(1e5,2e5,5e4)
    exoList = ExoDataSubsampling(file = files[3],sample.depth = sample.depth, verbose=FALSE)


## ----subsamplots,include=TRUE,echo=TRUE,eval=TRUE,warning=FALSE,results="markup",fig.height=2.5----
    
    p1 = paramDistBoxplot(exoList,which.param = "beta1")
    p2 = paramDistBoxplot(exoList,which.param = "beta2")
    gridExtra::grid.arrange(p1,p2,nrow = 1)


## ----info,include=TRUE,echo =TRUE,eval=TRUE-------------------------------------------------
sessionInfo("ChIPexoQual")

