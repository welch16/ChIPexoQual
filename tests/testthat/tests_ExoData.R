
library(ChIPexoQual)
library(ChIPexoQualExample)
library(parallel)

context("ExoData object is generated correctly with reads and files")

files = list.files(system.file("extdata",package = "ChIPexoQualExample"),
                   full.names = TRUE)
files = files[grep("bai",files,invert = TRUE)]
reads = mclapply(files,readGAlignments,param = NULL,
                 mc.cores = getOption("mc.cores",2L))
data("exampleExoData")

depths = sapply(exampleExoData,nreads)
reads_exo = lapply(reads,function(x)ExoData(reads = x))

test_that("Error on ExoData non-arguments",
          {
              expect_error(ExoData())
          })

test_that("Error on ExoData double arguments",
          {
              expect_error(ExoData(file = files[1],reads = reads[[1]]))
              expect_error(ExoData(file = files[2],reads = reads[[2]]))
              expect_error(ExoData(file = files[3],reads = reads[[3]]))
          })

test_that("Error on ExoData with missing file",
          {
              expect_error(ExoData(tempfile()))
          })

test_that("ExoData builds the same object with reads and file",
          {
              expect_identical(
                  GRanges(exampleExoData[[1]]),
                  GRanges(reads_exo[[1]])
              )
              expect_identical(
                  GRanges(exampleExoData[[2]]),
                  GRanges(reads_exo[[2]])
              )
              expect_identical(
                  GRanges(exampleExoData[[3]]),
                  GRanges(reads_exo[[3]])
              )
          })

test_that("ExoData is s4 class",
          {
              expect_s4_class(reads_exo[[1]],"ExoData")
              expect_s4_class(reads_exo[[2]],"ExoData")
              expect_s4_class(reads_exo[[3]],"ExoData")
          })


test_that("ExoData build same object when reads is GRanges and GAlignments",
          {
              expect_identical(
                  GRanges(exampleExoData[[1]]),
                  GRanges(ExoData(reads = as(reads[[1]],"GRanges")))
              )
          })

test_that("ExoData regions are disjoint",
          {
              expect_identical(
                  length(reads_exo[[1]]),
                  length(reduce(GRanges(reads_exo[[1]])))
              )
              expect_identical(
                  length(reads_exo[[2]]),
                  length(reduce(GRanges(reads_exo[[2]])))
              )
              expect_identical(
                  length(reads_exo[[3]]),
                  length(reduce(GRanges(reads_exo[[3]])))
              )
          })


test_that("At least one position in every island",
          {
              expect_true(
                  all(reads_exo[[1]]$u >= 1)
              )
              expect_true(
                  all(reads_exo[[2]]$u >= 1)
              )
              expect_true(
                  all(reads_exo[[3]]$u >= 1)
              )
          })

test_that("Number of unique positions in islands is <= depth",
          {
              expect_true(
                  all(reads_exo[[1]]$u <= reads_exo[[1]]$d   )
              )
              expect_true(
                  all(reads_exo[[2]]$u <= reads_exo[[2]]$d)
              )
              expect_true(
                  all(reads_exo[[3]]$u <= reads_exo[[3]]$d)
              )
              
          })

test_that("Writing correct nr. of reads in metadata",
          {
              expect_identical(
                  nreads(reads_exo[[1]]),depths[[1]])
              expect_identical(
                  nreads(reads_exo[[2]]),depths[[2]])
              expect_identical(
                  nreads(reads_exo[[3]]),depths[[3]])
          })

test_that("Reads returns the same .MA_DF output",
          {
              expect_identical(
                  .MA_DF(reads_exo[[1]]),.MA_DF(exampleExoData[[1]]))
              expect_identical(
                  .MA_DF(reads_exo[[2]]),.MA_DF(exampleExoData[[2]]))
              expect_identical(
                  .MA_DF(reads_exo[[3]]),.MA_DF(exampleExoData[[3]]))
              
          })

test_that("Reads returns the same .ARC_URC_DF output",
          {
              expect_identical(
                  .ARC_URC_DF(reads_exo[[1]],TRUE),
                  .ARC_URC_DF(exampleExoData[[1]],TRUE))
              expect_identical(
                  .ARC_URC_DF(reads_exo[[2]],TRUE),
                  .ARC_URC_DF(exampleExoData[[2]],TRUE))
              expect_identical(
                  .ARC_URC_DF(reads_exo[[3]],FALSE),
                  .ARC_URC_DF(exampleExoData[[3]],FALSE))
          })

test_that("Reads returns the same .FSR_dist_DF output",
          {
              quantiles = c(.25,.5,.75)
              depth_values = seq_len(25)
              expect_identical(
                  .FSR_dist_DF(reads_exo[[1]],quantiles = quantiles,
                               depth_values = depth_values,
                               both_strand = TRUE),
                  .FSR_dist_DF(exampleExoData[[1]],quantiles,
                                depth_values,TRUE))
              expect_identical(
                  .FSR_dist_DF(reads_exo[[2]],quantiles = quantiles,
                               depth_values = depth_values,
                               both_strand = FALSE),
                  .FSR_dist_DF(exampleExoData[[2]],quantiles,
                               depth_values,FALSE))
          })

test_that("Reads returns the same .region_comp_DF output",
          {
              depth_values = seq_len(10)
              expect_identical(
                  .region_comp_DF(reads_exo[[1]],depth_values),
                  .region_comp_DF(exampleExoData[[1]],depth_values)
              )
              expect_identical(
                  .region_comp_DF(reads_exo[[2]],depth_values),
                  .region_comp_DF(exampleExoData[[2]],depth_values)
              )
              expect_identical(
                  .region_comp_DF(reads_exo[[3]],depth_values),
                  .region_comp_DF(exampleExoData[[3]],depth_values)
              )
          })

