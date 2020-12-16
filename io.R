load_run <- function(filename){
  run_head <- read_csv(filename, col_names = c("line", "backup"),
                          col_types = cols(line = col_character(), 
                                           backup = col_character()))
  run <- read_csv(filename, skip = 24, col_types = "iccccccccdccddddcdcii")
  
  date_pattern <- "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2} \\w{2} \\w{3}"
  id_pattern <- "(?<=([Rr][Uu][Nn]))\\d+"
  
  run_head <- run_head %>%
    mutate(date = str_extract(line, date_pattern),
           id = str_extract(line, id_pattern),
           date_backup = str_extract(backup, date_pattern),
           id_backup = str_extract(backup, id_pattern))
  
  run_id <- run_head$id[1]
  run_date <- run_head$date[11]
  if (is.na(run_id)){
    run_id <- run_head$id_backup[1]
    run_date <- run_head$date_backup[11]
  }
  
  run <- run %>%
    clean_names() %>%
    rename(barcode = sample) %>%
    select(barcode, target, cq) %>%
    drop_na(barcode) %>%
    pivot_wider(names_from = target, values_from = cq) %>%
    clean_names()
  
  spread_repeated_barcodes <- function(barcode, rnasep, n1gene){
    df <- tibble(barcode = rep(barcode, times = length(rnasep)),
                 n1gene = n1gene,
                 rnasep = rnasep)
    return(df)
  }
  
  run <- pmap_dfr(run, spread_repeated_barcodes)
  
  run <- run %>%
    mutate(test_date = as_datetime(run_date),
           run_id = as.integer(run_id))
  return(run)
}


check_run_data <- function(run){
  qc = "PASS"
  qc_str = "No errors found"
  if(is.null(run)){
    qc = "FAIL"
    qc_str = "Could not load test run file!"
  } else if (typeof(run$run_id[[1]]) != "integer" | is.na(run$run_id[[1]])) {
    qc = "FAIL"
    qc_str = "No correct run ID identified!"
  }
  return(list(qc = qc, qc_str = qc_str))
}


load_sample_manifest <- function(filename){
  sample_manifest <- read_excel(filename, col_types = c("text", "date", "text", 
                                                        "text"))
  sample_manifest <- sample_manifest %>%
    clean_names()
  return(sample_manifest)
}


check_sample_manifest <- function(sample_manifest, run){
  qc = "PASS"
  qc_str = "No errors found"
  if (nrow(sample_manifest) == 0){
    qc = "FAIL"
    qc_str = "Sample Manifest is empty!"
  } else if(!any(run$barcode %in% sample_manifest$barcode)){
    qc = "FAIL"
    qc_str = "No Test Run barcodes found in Sample Manifest!"
  }
  return(list(qc = qc, qc_str = qc_str))
}


load_sample_accession <- function(filename){
  sample_accession <- read_excel(filename, col_types = c("text", "date", "text", 
                                                         "text"))
  sample_accession <- sample_accession %>%
    clean_names()
  return(sample_accession)
}


check_sample_accession <- function(sample_accession){
  qc = "PASS"
  qc_str = "No errors found"
  if (is.null(sample_accession)){
    qc = "PASS"
    qc_str = "No sample accession data loaded!"
  } else if (nrow(sample_accession) == 0){
    qc = "FAIL"
    qc_str = "Sample accession is empty!"
  }
  return(list(qc = qc, qc_str = qc_str))
}


load_subject_information <- function(filename){
  sample_information <- read_excel(filename,
                            col_types = c("text", "date", "date", "text",
                                          "text", "text", "date", "text",
                                          "text", "text", "text", "numeric",
                                          "text", "text", "text", "text"))
  sample_information <- sample_information %>%
    clean_names() %>%
    rename(collection_time = collection_time_hh_mm_am_pm,
           collection_date = date) %>%
    mutate(collection_time = hm(paste(hour(collection_time),
                                      minute(collection_time),
                                      sep = ":")),
           collection_date = date(collection_date),
           dob = date(dob))
  return(sample_information)
}


check_subject_info <- function(subject_info, run){
  qc = "PASS"
  qc_str = "No errors found"
  run <- run %>%
    mutate(barcode_check = ifelse(barcode %in% subject_info$barcode,
                                  TRUE, FALSE))
  if (!any(run$barcode_check)) {
    qc = "FAIL"
    qc_str = "No test run barcode matches subject information file!"
  } else if (!all(run$barcode_check)) {
    qc = "PASS"
    qc_str = paste("Samples with no subject information:",
                   paste(run$barcode[run$barcode_check == FALSE], 
                         collapse = ", "))
  }
  return(list(qc = qc, qc_str = qc_str))
}


load_covid_db <- function(filename){
  # TODO: replace with a database queries
  results <- read_excel(filename, sheet = "results",
                        col_types = c("text", "numeric", "numeric", "numeric",
                                      "text"))
  runs <- read_excel(filename, sheet = "runs",
                     col_types = c("numeric", "date"))
  return(list(results = results, runs = runs))
}


check_covid_db_integrity <- function(covid_db){
  runs <- covid_db$runs
  results <- covid_db$results
  qc = "PASS"
  qc_str = "No database errors found"
  if (nrow(runs) == 0){
    qc = "FAIL"
    qc_str = "Covid DB is empty!"
  } else if (unique(runs$run_id) != runs$run_id){
    qc = "FAIL"
    qc_str = "Duplicate run ID's found in 'runs' table!"
  } else if (unique(results$run_id) != runs$run_id){
    qc = "FAIL"
    qc_str = "Run ID mismatch between 'results' and 'runs' table!"
  }
  return(list(qc = qc, qc_str = qc_str))
}