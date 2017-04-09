library(ChIPexoQual)
library(parallel)

context("ExoDataBlacklist returns errors")

data("exoExample")
data("blacklists")

test_that("Error on ExoDataBlacklist non-arguments",
          {
            expect_error(ExoDataBlacklist())
          })

test_that("Error when blacklist is not GRanges",
          {
            expect_error(ExoDataBlacklist(exoExample,"mm9"))
          })

test_that("Error when which.param is not beta1 or beta2",
          {
            expect_error(ExoDataBlacklist(exoExample,blacklists[["mm9"]],"beta3"))
          })

test_that("Error when nregions is bigger than regions in ExoData object",
          {
            expect_error(ExoDataBlacklist(exoExample,blacklists[["mm9"]],nregions = lenth(exo) + 3))
          })