source("../io.R", chdir = TRUE)
library(testthat)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)

run_standard <- function(filename){
  run_standard <- read_csv(filename, col_types = "cddcd")
  run_standard <- run_standard %>%
    mutate(test_date = as.Date(test_date))
  return(run_standard)
}

test_that("run_data", {
  expect_equal(load_run("test_input_run1.csv"), 
               run_standard("test_output_run1.csv"))
})