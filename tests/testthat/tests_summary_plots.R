
library(ChIPexoQual)
library(ChIPexoQualExample)

context("Plotting functions errors and utilities")

data("exoExample")

test_that("MA_plot error when the length of names is different from input",
          {
              n <- 2
              u <- as.character(runif(n))
              expect_error(
                  MAplot(exoExample,names.input = u)
              )
          })

test_that("ARC_URC_plot error when the length of names  
          is different from input",
          {
              n <- 2
              u <- as.character(runif(n))
              expect_error(
                  ARCvURCplot(exoExample,names.input = u)
              )
          })

test_that("FSR_dist_plot error when the length of names is different 
          from input",
          {
              n <- 2
              u <- as.character(runif(n))
              expect_error(
                  FSRDistplot(exoExample,names.input = u)
              )
          })

test_that("region_comp_plot error when the length of names is different 
          from input",
          {
              n <- 2
              u <- as.character(runif(n))
              expect_error(
                  regionCompplot(exoExample,names_input = u)
              )

          })

test_that("param_dist_boxplot error when the length of names is different 
          from input",
          {
              n <- 2
              u <- as.character(runif(n))
              expect_error(
                  paramDistBoxplot(exoExample,names_input = u)
              )
          })

test_that("param_dist_boxplot error when sort.as.numeric is not logical",
          {
              n <- 2
              u <- as.character(runif(n))
              expect_error(
                  paramDistBoxplot(exoExample,names_input = u , sort.as.numeric = "NO")
              )
          })