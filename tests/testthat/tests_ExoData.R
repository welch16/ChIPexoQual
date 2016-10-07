
library(ChIPexoQual)
library(ChIPexoQualExample)
library(parallel)

context("ExoData object is generated correctly with reads and files")

files = list.files(system.file("extdata",package = "ChIPexoQualExample"),
                   full.names = TRUE)
files = files[grep("bai",files,invert = TRUE)]
files = files[3]
reads = readGAlignments(files,param = NULL)
data("exoExample")

depths = nreads(exoExample)
reads_exo = ExoData(reads = reads,mc.cores = 2L)

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

test_that("Non-error with save_reads argument",
          {
              expect_s4_class(
                  ExoData(file = files,save_reads = TRUE),
                  "ExoData"
              )
          })

test_that("At least one position in every island",
          {
              expect_true(
                  all(reads_exo$u >= 1)
              )
          })

test_that("Number of unique positions in islands is <= depth",
          {
              expect_true(
                  all(reads_exo$u <= reads_exo$d)
              )
              
          })

test_that("Writing correct nr. of reads in metadata",
          {
              expect_identical(
                  nreads(reads_exo),depths)
          })

test_that("Reads returns the same .MA_DF output",
          {
              expect_identical(
                  .MA_DF(reads_exo),.MA_DF(exoExample))
              
          })

test_that("Reads returns the same .ARC_URC_DF output",
          {
              expect_identical(
                  .ARC_URC_DF(reads_exo,FALSE),
                  .ARC_URC_DF(exoExample,FALSE))
          })

test_that("Reads returns the same .FSR_dist_DF output",
          {
              quantiles = c(.25,.5,.75)
              depth_values = seq_len(25)
              expect_identical(
                  .FSR_dist_DF(reads_exo,quantiles = quantiles,
                               depth_values = depth_values,
                               both_strand = FALSE),
                  .FSR_dist_DF(exoExample,quantiles,
                               depth_values,FALSE))
          })

test_that("Reads returns the same .region_comp_DF output",
          {
              depth_values = seq_len(10)
              expect_identical(
                  .region_comp_DF(reads_exo,depth_values),
                  .region_comp_DF(exoExample,depth_values)
              )
          })

