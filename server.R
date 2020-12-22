#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(openxlsx)
library(rmarkdown)
source("io.R")
source("data_functions.R")

# cutpoints
cutpoint_rnasep   <<- 30
cutpoint_n1gene   <<- 36

# strings for output categories
positive_str      <<- "detected"     # N1 = Amp; P = Amp
negative_str      <<- "not detected" # N1 = no Amp; P = Amp
indeterminate_str <<- "indeterminate, pending repeat analysis" # P = no Amp
bad_str           <<- "specimen not adequate, need to recollect"
pending_str       <<- "pending"
not_received_str  <<- "sample not received"

# strings for status
final_str         <<- "final"
preliminary_str   <<- "preliminary"
scheduled_str     <<- "scheduled"

# string if report failed
fail_str          <<- "report generation fail"

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

    run_data_raw <- eventReactive(input$run_script, {
      file <- input$run
      run_check <- read_csv(file$datapath, skip = 24)
      t <- run_check
      validate(
        need(colnames(run_check)[1] == "Well", 
             "Corrupt data format in test run file!")
      )
      run_data_raw <- load_run(file$datapath)
    })
    
    run_data_qc <- reactive({
      qc <- check_run_data(run_data_raw())
    })
    
    run_data <- reactive({
      validate(need(run_data_qc()$qc != "FAIL", run_data_qc()$qc_str))
      run_data <- run_data_raw()
    })
    
    sample_manifest_raw <- eventReactive(input$run_script, {
      file <- input$sample_manifest
      sample_manifest_raw <- load_sample_manifest(file$datapath)
    })
    
    sample_manifest_qc <- reactive({
      qc <- check_sample_manifest(sample_manifest_raw(), run_data())
    })
    
    sample_manifest <- reactive({
      validate(need(sample_manifest_qc()$qc != "FAIL", 
                    sample_manifest_qc()$qc_str))
      sample_manifest <- sample_manifest_raw()
    })
    
    subject_info_raw <- eventReactive(input$run_script, {
      file <- input$subject_info
      subject_info <- load_subject_information(file$datapath)
    })
    
    subject_info_qc <- reactive({
      qc <- check_subject_info(subject_info_raw(), run_data())
    })
    
    subject_info <- reactive({
      validate(need(subject_info_qc()$qc != "FAIL", subject_info_qc()$qc_str))
      subject_info <- subject_info_raw()
    })
    
    covid_db_raw <- eventReactive(input$run_script, {  
      file <- input$covid_db
      covid_db_raw <- load_covid_db(file$datapath)
    })
    
    covid_db_qc <- reactive({
      qc <- check_covid_db_integrity(covid_db_raw())
    })
    
    covid_db <- reactive({
      validate(need(covid_db_qc()$qc != "FAIL", covid_db_qc()$qc_str))
      covid_db <- covid_db_raw()
    })
    
    sample_accession_raw <- eventReactive(input$run_script, {
      file <- input$sample_accession
        if (!is.null(file)){
        sample_accession_raw <- load_sample_manifest(file$datapath)
      }
    })
    
    sample_accession_qc <- reactive({
      qc <- check_sample_accession(sample_accession_raw())
    })
    
    sample_accession <- reactive({
      validate(need(sample_accession_qc()$qc != "FAIL", 
                    sample_accession_qc()$qc_str))
      sample_accession <- sample_accession_raw()
    })
    
    indiv_report_source_qc <- reactive({
      file <- input$indiv_report_source
      indiv_report_source <- NULL
      if (!is.null(file)){
        indiv_report_source <- load_indiv_report_source(file$datapath)
      }
      qc <- check_indiv_report_source(indiv_report_source)
    })
    
    run_results <- reactive({
      validate(
        need(!is.null(run_data()) & !is.null(sample_manifest()),
                      "Can't create results - problem loading run data or sample manifest!")
      )
      run_formatted <- format_run(run_data())
      run_results <- compute_run_results(run_formatted)
      run_results <- add_failed_manifest_qc_samples_to_results(run_results,
                                                               sample_manifest())
    })
    
    matched_results <- reactive({
      validate(
        need(!is.null(run_results()),
             "Can't match previous runs with current run!")
      )
      matched_results <- match_run_results_with_previous(run_results(), 
                                                         covid_db()$results)
      matched_results <- set_result_status(matched_results)
      if (!is.null(sample_accession())){
        matched_results <- add_failed_accession_qc_samples_to_results(matched_results,
                                                                      sample_accession())
      }
      matched_results
    })

    updated_covid_db_raw <- reactive({
      validate(
        need(!is.null(run_data()) &
               !is.null(covid_db()) &
               !is.null(run_results()),
             "Can't update database!")
      )
      covid_db_data <- add_run_to_covid_db(covid_db(), run_results(), 
                                           run_data())
    })
    
    updated_covid_db_qc <- reactive({
      qc <- updated_covid_db_raw()$qc
    })
    
    updated_covid_db <- reactive({
      validate(need(updated_covid_db_qc()$qc != "FAIL", 
                    updated_covid_db_qc()$qc_str))
      updated_covid_db <- updated_covid_db_raw()$covid_db
    })
    
    report_data <- reactive({
      validate(
        need(!is.null(matched_results()) & !is.null(subject_info()) &
               !is.null(updated_covid_db()),
             "Can't create report - result matching error or subject info not found!")
      )
      report_data <- select_results_for_report(matched_results(), 
                                               subject_info())
    })
    
    report_condensed <- reactive({
      validate(
        need(!is.null(report_data()) & !is.null(subject_info()),
             "Can't create condensed report!")
      )
      report_with_subject_info <- add_subject_info_to_report(report_data()$final, 
                                                             subject_info())
      report_condensed <- format_condensed_report(report_with_subject_info)
    })
    
    report_medicat <- reactive({
      validate(
        need(!is.null(report_data()) & !is.null(subject_info()),
             "Can't create medicat report!")
      )
      report_with_subject_info <- add_subject_info_to_report(report_data()$final, 
                                                             subject_info())
      report_medicat <- format_medicat_report(report_with_subject_info)
    })
    
    report_vdh <- reactive({
      validate(
        need(!is.null(report_data()) & !is.null(subject_info()),
             "Can't create vdh report!")
      )
      report_with_subject_info <- add_subject_info_to_report(report_data()$final, 
                                                             subject_info())
      vdh_template <- create_vdh_template()
      report_vdh <- format_vdh_report(vdh_template, report_with_subject_info)
    })
    
    output$run_file_qc <- renderText(
      paste0("Run file loaded: ", run_data_qc()$qc, " (", 
             run_data_qc()$qc_str, ")")
    )
    
    output$sample_manifest_qc <- renderText(
      paste0("Sample Manifest loaded: ", sample_manifest_qc()$qc, " (",
             sample_manifest_qc()$qc_str, ")")
    )
    
    output$covid_db_qc <- renderText(
      paste0("Covid DB integrity: ", covid_db_qc()$qc, " (",
             covid_db_qc()$qc_str, ")")
    )
    
    output$subject_info_qc <- renderText(
      paste0("Subject information loaded: ", subject_info_qc()$qc, " (",
             subject_info_qc()$qc_str, ")")
    )
    
    output$sample_accession_qc <- renderText(
      paste0("Sample accession information: ", sample_accession_qc()$qc, " (",
             sample_accession_qc()$qc_str, ")")
    )
    
    output$updated_covid_db_qc <- renderText(
      paste0("Updating Covid DB: ", updated_covid_db_qc()$qc, " (",
             updated_covid_db_qc()$qc_str, ")")
    )
    
    output$run_info <- renderText(
      paste0("Run ID: ", run_data()$run_id[[1]], " Test Date: ",
             run_data()$test_date[[1]]
             )
    )
    
    output$indiv_report_source_qc <- renderText(
      paste0("Individual Reports Source File: ", indiv_report_source_qc()$qc, 
             " (", indiv_report_source_qc()$qc_str, ")")
    )
    
    output$final_report <- renderTable({
        report_data()$final
    })
    
    output$prelim_report <- renderTable({
        report_data()$prelim
    })

    output$download_condensed_report <- downloadHandler(
        filename = function () {
            paste0("Report_",
                   format(Sys.time(), "%m%d%y_%k%M%S"),
                   ".xlsx")
        },
        content = function(file) {
            write.xlsx(report_condensed(), file)
        })
    
    output$download_qc_report <- downloadHandler(
        filename = function () {
            paste0("Report_QC_",
                   format(Sys.time(), "%m%d%y_%k%M%S"),
                   ".xlsx")
        },
        content = function(file) {
            datasets <- list("Status_Final" = report_data()$final, 
                             "Status_Preliminary" = report_data()$prelim)
            write.xlsx(datasets, file)
        })
    
    output$download_medicat_report <- downloadHandler(
      filename = function () {
        paste0("Report_Medicat_",
               format(Sys.time(), "%m%d%y_%k%M%S"),
               ".csv")
      },
      content = function(file) {
        write_delim(report_medicat(), file, delim = "|", col_names = FALSE)
      })
    
    output$download_vdh_report <- downloadHandler(
      filename = function () {
        paste0("Report_VDH_",
               format(Sys.time(), "%m%d%y_%k%M%S"),
               ".txt")
      },
      content = function(file) {
        write_delim(report_vdh(), file, delim = "|", col_names = FALSE)
      })
    
    output$download_covid_db <- downloadHandler(
      filename = function () {
        covid_db_runs <- updated_covid_db()$runs
        latest_run <- max(covid_db_runs$run_id)
        paste0("Covid_DB_up_to_run_",
               latest_run,
               "_",
               format(Sys.time(), "%m%d%y_%k%M%S"),
               ".xlsx")
      },
      content = function(file) {
        datasets <- list("results" = updated_covid_db()$results,
                         "runs" = updated_covid_db()$runs)
        write.xlsx(datasets, file)
      })
    
    output$download_indiv_reports <- downloadHandler(
      filename = function() {
        paste0("individual_reports", ".zip")
      },
      content = function(file) {
        validate(need(indiv_report_source_qc()$qc != "FAIL", 
                      indiv_report_source_qc()$qc_str))
        detach("package:openxlsx", unload=TRUE)
        library(zip)
        source_file = input$indiv_report_source
        indiv_report_source <- load_indiv_report_source(source_file$datapath)
        n <- nrow(indiv_report_source)
        reports <- c()
        withProgress(message = 'Creating Invidivual Report:', value = 0, {
          for (row in 1:nrow(indiv_report_source)){
            barcode         <- indiv_report_source[[row, "barcode"]]
            first_name      <- indiv_report_source[[row, "first_name"]]
            last_name       <- indiv_report_source[[row, "last_name"]]
            test_result     <- indiv_report_source[[row, "test_result"]]
            test_date       <- indiv_report_source[[row, "test_date"]]
            collection_date <- indiv_report_source[[row, "collection_date"]]
            
            incProgress(1/n, detail = paste(first_name, last_name))
            
            path <- paste0(first_name, "_", last_name, "_", test_date, ".pdf")
            render("individual_report.Rmd", 
                   output_file = path,
                   params=list(barcode = barcode,
                               first_name = first_name,
                               last_name = last_name,
                               test_result = test_result,
                               test_date = test_date,
                               collection_date = collection_date))
            
            reports <- c(reports, path)
          }})
          zip(file, reports)
      },
      contentType = "application/zip"
    )
})
