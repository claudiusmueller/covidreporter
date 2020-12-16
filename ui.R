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
            
            fileInput(
                "sample_manifest",
                "New Sample Manifest",
                multiple = FALSE,
                accept = c(".xlsx", ".xls"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            fileInput(
                "subject_info",
                "Subject Information",
                multiple = FALSE,
                accept = c(".xlsx", ".xls", ".xlm", ".xlsm"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            fileInput(
                "covid_db",
                "Covid Database",
                multiple = FALSE,
                accept = c(".xlsx"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            fileInput(
                "sample_accession",
                "New Sample Accession Data (optional)",
                multiple = FALSE,
                accept = c(".xlsx", ".xls"),
                width = NULL,
                buttonLabel = "Browse...",
                placeholder = "No file selected"
            ),
            
            actionButton(
                "run_script",
                "Update database and create report",
                icon = NULL,
                width = NULL
            ),
            
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
            
            HTML("<br><br><br>"),
            
            textOutput("run_file_qc"),
            textOutput("sample_manifest_qc"),
            textOutput("covid_db_qc"),
            textOutput("subject_info_qc"),
            textOutput("sample_accession_qc"),
            textOutput("updated_covid_db_qc")

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