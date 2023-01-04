library(knitr)
library(rmarkdown)
library(tidyverse)
library(readxl)
library(janitor)
setwd("/data/programming/Work/Covid_CAP_Lab/covidreporter")
source("io.R")
source_file <- "/data/work/Projects - Active/CAP/Covid/GMU_Virus/data_pipeline/test_shiny_2/example_indiv_report_input.xlsx"
indiv_report_source <- load_indiv_report_source(source_file)
for (row in 1:nrow(indiv_report_source)){
    barcode         <- indiv_report_source[[row, "barcode"]]
    first_name      <- indiv_report_source[[row, "first_name"]]
    last_name       <- indiv_report_source[[row, "last_name"]]
    test_result     <- indiv_report_source[[row, "test_result"]]
    test_date       <- indiv_report_source[[row, "test_date"]]
    collection_date <- indiv_report_source[[row, "collection_date"]]

    render("individual_report.Rmd", 
           output_file = paste0("tests/", first_name, "_", last_name, 
                                "_", test_date, ".pdf"),
           params=list(barcode = barcode,
                       first_name = first_name,
                       last_name = last_name,
                       test_result = test_result,
                       test_date = test_date,
                       collection_date = collection_date))
}
