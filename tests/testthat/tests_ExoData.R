
library(ChIPexoQual)
library(ChIPexoQualExample)
library(parallel)

context("ExoData object is generated correctly with reads and files")

files <- list.files(system.file("extdata",package = "ChIPexoQualExample"),
                   full.names = TRUE)
files <- files[grep("bai",files,invert = TRUE)]
files <- files[3]
reads <- readGAlignments(files,param = NULL)
data("exoExample")

depths <- nreads(exoExample)
reads_exo <- ExoData(reads = reads,mc.cores = 2L,verbose = FALSE)

test_that("Error on ExoData non-arguments",
          {
              expect_error(ExoData())
          })

test_that("Error on ExoData double arguments",
          {
              expect_error(ExoData(file = files,reads = reads))
          })

test_that("Error on ExoData with missing file",
          {
              expect_error(ExoData(tempfile()))
          })

test_that("ExoData is s4 class",
          {
              expect_s4_class(reads_exo,"ExoData")
          })

test_that("Non-error with save.reads argument",
          {
              expect_s4_class(
                  ExoData(file = files,save.reads = TRUE),
                  "ExoData"
              )
          })

test_that("At least one position in every island",
          {
              expect_true(
                  all(reads_exo$uniquePos >= 1)
              )
          })

test_that("Number of unique positions in islands is <= depth",
          {
              expect_true(
                  all(reads_exo$uniquePos <= reads_exo$depth)
              )
              
          })

test_that("Writing correct nr. of reads in metadata",
          {
              expect_identical(
                  nreads(reads_exo),depths)
          })

test_that("Reads returns the same MADataFrame output",
          {
              expect_equal(
                  MADataFrame(reads_exo),MADataFrame(exoExample))
              
          })

test_that("Reads returns the same ARCvURCDataFrame output",
          {
              expect_equal(
                  ARCvURCDataFrame(reads_exo,FALSE),
                  ARCvURCDataFrame(exoExample,FALSE))
          })

test_that("Reads returns the same FSRDistDataFrame output",
          {
              quantiles <- c(.25,.5,.75)
              depth_values <- seq_len(25)
              expect_equal(
                  FSRDistDataFrame(reads_exo,quantiles = quantiles,
                               depth.values = depth_values,
                               both.strand = FALSE),
                  FSRDistDataFrame(exoExample,quantiles,
                               depth_values,FALSE))
          })

test_that("Reads returns the same regionCompDataFrame output",
          {
              depth_values <- seq_len(10)
              expect_equal(
                  regionCompDataFrame(reads_exo,depth_values),
                  regionCompDataFrame(exoExample,depth_values)
              )
          })

test_that("Error when the provided sample.depth is greater than available 
          number of reads",
          {
              expect_error(ExoDataSubsampling(reads = reads_exo,
                                              sample.depth = depths + 100))
              
          })

test_that("Error when sample.depth is not numeric",
          {
              expect_error(ExoDataSubsampling(reads = reads_exo,
                                              sample.depth = c("1","2")))
          })

