
library(ChIPexoQual)
library(ChIPexoQualExample)

context("Plotting functions errors and utilities")

data("exampleExoData")

test_that("Generate the names correctly",
          {
              n = 5
              expect_identical(
                  .generate_names(NULL,NULL,5),
                  paste("Sample:",1:n,sep = " ")
              )
              u = as.character(runif(n))
              expect_identical(
                  .generate_names(NULL,u,n),
                  u
              )
              expect_identical(
                  .generate_names(u,NULL,n),
                  u
              )
          })

test_that("MA_plot error when the length of names is different from input",
          {
              n = 2
              u = as.character(runif(n))
              expect_error(
                  MA_plot(exampleExoData,names_input = u)
              )
              expect_error(
                  MA_plot(exampleExoData[[1]],names_input = u)
              )
              expect_error(
                  MA_plot(exampleExoData[[1]],
                          exampleExoData[[2]],
                          exampleExoData[[3]],
                          names_input = u)
              )
          })

test_that("ARC_URC_plot error when the length of names  
          is different from input",
          {
              n = 2
              u = as.character(runif(n))
              expect_error(
                  ARC_URC_plot(exampleExoData,names_input = u)
              )
              expect_error(
                  ARC_URC_plot(exampleExoData[[1]],names_input = u)
              )
              expect_error(
                  ARC_URC_plot(exampleExoData[[1]],
                          exampleExoData[[2]],
                          exampleExoData[[3]],
                          names_input = u)
              )
          })

test_that("FSR_dist_plot error when the length of names is different 
          from input",
          {
              n = 2
              u = as.character(runif(n))
              expect_error(
                  FSR_dist_plot(exampleExoData,names_input = u)
              )
              expect_error(
                  FSR_dist_plot(exampleExoData[[1]],names_input = u)
              )
              expect_error(
                  FSR_dist_plot(exampleExoData[[1]],
                          exampleExoData[[2]],
                          exampleExoData[[3]],
                          names_input = u)
              )
              
          })

test_that("region_comp_plot error when the length of names is different 
          from input",
          {
              n = 2
              u = as.character(runif(n))
              expect_error(
                  region_comp_plot(exampleExoData,names_input = u)
              )
              expect_error(
                  region_comp_plot(exampleExoData[[1]],names_input = u)
              )
              expect_error(
                  region_comp_plot(exampleExoData[[1]],
                          exampleExoData[[2]],
                          exampleExoData[[3]],
                          names_input = u)
              )
              
          })

test_that("param_dist_boxplot error when the length of names is different 
          from input",
          {
              n = 2
              u = as.character(runif(n))
              expect_error(
                  param_dist_boxplot(exampleExoData,names_input = u)
              )
              expect_error(
                  param_dist_boxplot(exampleExoData[[1]],names_input = u)
              )
              expect_error(
                  param_dist_boxplot(exampleExoData[[1]],
                          exampleExoData[[2]],
                          exampleExoData[[3]],
                          names_input = u)
              )
              
          })