#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("GMU Clinical Proteomics Lab - Covid Report"),
    
    sidebarLayout(
        
        sidebarPanel(
            
            fileInput(
                "run",
                "New Test Run Data",
                multiple = FALSE,
                accept = ".csv",
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            textOutput("run_file_qc"),
            HTML("<br>"),
            
            fileInput(
                "sample_manifest",
                "New Sample Manifest",
                multiple = FALSE,
                accept = c(".xlsx", ".xls"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            textOutput("sample_manifest_qc"),
            HTML("<br>"),
            
            fileInput(
                "subject_info",
                "Subject Information",
                multiple = FALSE,
                accept = c(".xlsx", ".xls", ".xlm", ".xlsm"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            textOutput("subject_info_qc"),
            HTML("<br>"),
            
            fileInput(
                "covid_db",
                "Covid Database",
                multiple = FALSE,
                accept = c(".xlsx"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            textOutput("covid_db_qc"),
            HTML("<br>"),
            
            # fileInput(
            #     "sample_accession",
            #     "New Sample Accession Data (optional)",
            #     multiple = FALSE,
            #     accept = c(".xlsx", ".xls"),
            #     width = NULL,
            #     buttonLabel = "Browse...",
            #     placeholder = "No file selected"
            # ),
            # 
            # textOutput("sample_accession_qc"),
            # HTML("<br>"),
            
            actionButton(
                "run_script",
                "Update database and create report",
                icon = NULL,
                width = NULL
            ),
            
            HTML("<br><br>"),
            textOutput("updated_covid_db_qc"),
            textOutput("results_qc"),
            
            HTML("<br><br><br>"),
            
            downloadButton(
                "download_condensed_report", 
                "Download Condensed Report"),
            
            downloadButton(
                "download_qc_report", 
                "Download QC Report"),
            
            downloadButton(
                "download_medicat_report",
                "Download Medicat Report"
            ),
            
            downloadButton(
                "download_vdh_report",
                "Download VDH Report"
            ),
            
            downloadButton(
                "download_covid_db", 
                "Download Updated Covid DB"),

            HTML("<br>"),
            
            hr(style = "border-top: 1px solid #000000;"),
            
            HTML("<br><br><br>"),
            
            fileInput(
                "indiv_report_source",
                "Source File for Individual Reports",
                multiple = FALSE,
                accept = c(".xlsx", ".xls"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            textOutput("indiv_report_source_qc"),
            HTML("<br>"),
            
            downloadButton(
                "download_indiv_reports",
                "Create and Download Individual Reports"
            )

        ),

        mainPanel(
            
            textOutput("run_info"),
            h1("Test Results (Final)"),
            tableOutput("final_report"),
            h1("Test Results (Preliminary)"),
            tableOutput("prelim_report")
            )
        
    )
))
