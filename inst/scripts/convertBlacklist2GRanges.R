
library(GRanges)

indir = system.file("extdata",package = "ChIPexoQual")
files = list.files(indir,full.names = TRUE)

blacklists = lapply(files,readr::read_tsv,col_names = FALSE)
names(blacklists) = sapply(strsplit(basename(files),"-"),function(x)x[1])

blacklists = lapply(blacklists,
                function(x)GRanges(seqnames = x$X1,
                                   ranges = IRanges(
                                     start = x$X2,
                                     end = x$X3
                                   )))

save(blacklists,file = "data/blacklists.RData")
