load_run <- function(filename){
  
  # TODO: remove this line once readr has been fixed. Readr ver2.0.0 currently 
  # crashes R. This line forces readr to use the first edition of readr.
  readr::local_edition(1)
  # ========================================================
  
  run_head <- read_csv(filename, col_names = c("line", "backup"),
                          col_types = cols(line = col_character(), 
                                           backup = col_character()))
  
  date_pattern <- "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2} \\w{2} \\w{3}"
  id_pattern <- "(?<=([Rr][Uu][Nn]))\\d+"
  
  run_head <- run_head %>%
    mutate(date = str_extract(line, date_pattern),
           id = str_extract(line, id_pattern),
           date_backup = str_extract(backup, date_pattern),
           id_backup = str_extract(backup, id_pattern),
           correct_date_pos = str_detect(line, "Run Start Date"),
           correct_id_pos = str_detect(line, "File Name"),
           correct_well_pos = str_detect(line, "^Well"))
  
  id_pos <- which(run_head$correct_id_pos == TRUE)
  date_pos <- which(run_head$correct_date_pos == TRUE)
  well_pos <- which(run_head$correct_well_pos == TRUE)
  
  run_id <- run_head$id[id_pos]
  run_date <- run_head$date[date_pos]
  if (is.na(run_id)){
    run_id <- run_head$id_backup[id_pos]
    run_date <- run_head$date_backup[date_pos]
  }
  
  run_data <- read_csv(filename, skip = well_pos - 1)
                  # col_types = "iccccccccdccddddcdcii")
  
  run_data <- run_data %>%
    clean_names() %>%
    rename(barcode = sample) %>%
    select(barcode, target, cq) %>%
    drop_na(barcode) %>%
    mutate(barcode = as.character(barcode),
           target = as.character(target),
           cq = as.double(cq)) %>%
    pivot_wider(names_from = target, values_from = cq) %>%
    clean_names()
  
  spread_repeated_barcodes <- function(barcode, rnasep, n1gene){
    df <- tibble(barcode = rep(barcode, times = length(rnasep)),
                 n1gene = n1gene,
                 rnasep = rnasep)
    return(df)
  }
  
  run_data <- pmap_dfr(run_data, spread_repeated_barcodes)
  
  run_data <- run_data %>%
    mutate(test_date = as_datetime(run_date),
           run_id = as.integer(run_id))
  return(run_data)
}


check_run_data <- function(run_data){
  if (is.null(run_data)){
    qc = "FAIL"
    qc_str = "No test run data loaded!"
  } else if(is.na(run_data)){
    qc = "FAIL"
    qc_str = "Could not load test run file!"
  } else if (typeof(run_data$run_id[[1]]) != "integer" | 
             is.na(run_data$run_id[[1]])) {
    qc = "FAIL"
    qc_str = "No correct run ID identified!"
  } else {
    qc = "PASS"
    qc_str = "No errors found"
  }
  return(list(qc = qc, qc_str = qc_str))
}


load_sample_manifest <- function(filename){
  sample_manifest <- read_excel(filename, sheet = "Sample_Manifest")
  sample_manifest <- sample_manifest %>%
    clean_names() %>%
    mutate(date = as_date(date),
           barcode = as.character(barcode))
  return(sample_manifest)
}


check_sample_manifest <- function(sample_manifest, run){
  if (is.null(sample_manifest)){
    qc = "FAIL"
    qc_str = "No sample manifest loaded!"
  } else if (is.na(sample_manifest)){
    qc = "FAIL"
    qc_str = "Failed to load sample manifest file!"
  } else if (is.na(run)){
    qc = "FAIL"
    qc_str = "No run data loaded (needed for sample manifest qc)"
  } else if (nrow(sample_manifest) == 0){
    qc = "FAIL"
    qc_str = "Sample manifest is empty!"
  } else if(!any(run$barcode %in% sample_manifest$barcode)){
    qc = "FAIL"
    qc_str = "No test run barcodes found in sample manifest!"
  } else {
    qc = "PASS"
    qc_str = "No errors found"
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
  subject_info <- read_excel(filename)
  subject_info <- subject_info %>%
    clean_names() %>%
    rename(barcode = barcode_number,
           dob = birth_date,
           netid = net_id) %>%
    mutate(collection_date = date(collection_date),
           dob = date(dob),
           barcode = as.character(barcode)) %>%
    select(barcode, g_number, first_name, last_name, dob, netid, barcode, 
           collection_date) %>%
    distinct(barcode, netid, last_name, first_name, .keep_all = TRUE)
  subject_info
  return(subject_info)
}


check_subject_info <- function(subject_info, run, covid_db){
  if (is.null(subject_info)){
    qc = "FAIL"
    qc_str = "No subject information loaded!"
  } else if (is.na(subject_info)){
    qc = "FAIL"
    qc_str = "Failed to load subject information file!"
  } else if (is.na(run)){
    qc = "FAIL"
    qc_str = "No run data loaded (needed for sample manifest qc)"
  } else {
    known_barcodes <- c(unique(subject_info$barcode), 
                        unique(covid_db$subject_info$barcode))
    run <- run %>%
      mutate(barcode_check = ifelse(barcode %in% known_barcodes,
                                    TRUE, FALSE))
    
    subject_info_dupl <- subject_info %>%
      group_by(barcode) %>%
      filter(n() > 1)
    
    subject_info_missing <- subject_info %>%
      mutate(across(.cols = everything(), ~as.character(.))) %>%
      rowwise() %>%
      mutate(flag = any(is.na(c_across(cols = everything())))) %>%
      ungroup() %>%
      filter(flag == TRUE)
    
    if (nrow(subject_info_dupl) > 0){
      qc = "FAIL"
      qc_str = paste("Found duplicate barcodes with mismatching subject information:",
                     paste(unique(subject_info_dupl$barcode), collapse = ", "))
    } else if (nrow(subject_info_missing) > 0){
      qc = "FAIL"
      qc_str = paste("Found barcodes with missing subject information:",
                     paste(subject_info_missing$barcode), collapse = ", ")
    } else if (!all(run$barcode_check)) {
      qc = "PASS"
      qc_str = paste("Samples with no subject information:",
                     paste(run$barcode[run$barcode_check == FALSE], 
                           collapse = ", "))
    } else {
      qc = "PASS"
      qc_str = "No errors found"
    }
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
  subject_info <- read_excel(filename, sheet = "subject_info",
                             col_types = c("text", "text", "text", "text",
                                           "text", "date", "date"))
  full <- left_join(results, runs, by = "run_id")
  return(list(results = results, runs = runs, subject_info = subject_info, 
              full = full))
}


check_covid_db_integrity <- function(covid_db){
  if (is.null(covid_db)){
    qc = "FAIL"
    qc_str = "No covid database loaded!"
  } else if (is.na(covid_db)){
    qc = "FAIL"
    qc_str = "Could not load covid database file!"
  } else {
    runs <- covid_db$runs
    results <- covid_db$results
    if (nrow(runs) == 0){
      qc = "FAIL"
      qc_str = "Covid DB is empty!"
    } else if (!identical(sort(unique(runs$run_id)), sort(runs$run_id))){
      qc = "FAIL"
      qc_str = "Duplicate run ID's found in 'runs' table!"
    } else if (!identical(sort(unique(results$run_id)), sort(runs$run_id))){
      qc = "FAIL"
      qc_str = "Run ID mismatch between 'results' and 'runs' table!"
    } else {
      qc = "PASS"
      qc_str = "No database errors found"
    }
  }
  return(list(qc = qc, qc_str = qc_str))
}


load_indiv_report_source <- function(filename){
  indiv_report_source <- read_excel(filename) %>%
    clean_names() %>%
    select(barcode, first_name, last_name, dob, collection_date, test_date, 
           test_result) %>%
    mutate(dob = as_date(dob),
           collection_date = as_date(collection_date),
           test_date = as_date(test_date))
  return(indiv_report_source)
}


check_indiv_report_source <- function(indiv_report_source){
  if (is.null(indiv_report_source)){
    qc = "FAIL"
    qc_str = "No report source file loaded!"
  } else if (is.na(indiv_report_source)){
    qc = "FAIL"
    qc_str = "Failed to load file!"
  } else if (nrow(indiv_report_source) < 1) {
    qc = "FAIL"
    qc_str = "Source file empty!"
  } else {
    qc = "PASS"
    qc_str = "No errors found"
  }
  return(list(qc = qc, qc_str = qc_str))
}


check_results <- function(matched_results, subject_info){
  check_barcodes <- matched_results %>%
    filter(barcode %in% subject_info$barcode) %>%
    left_join(subject_info, by = "barcode") %>%
    group_by(barcode) %>%
    filter(n() > 1)
  if (nrow(check_barcodes) > 0){
    qc = "FAIL"
    qc_str = paste("Barcodes assigned to more than one subject:", 
                   paste(unique(check_barcodes$barcode), collapse = ", "))
  } else {
    qc = "PASS"
    qc_str = "No errors found"
  }
  return(list(qc = qc, qc_str = qc_str))
}


check_controls <- function(run_data){
  negative_str <- "Negative"
  twist_str    <- "Twist positive"
  idt_str      <- "IDT positive"
  
  results <- run_data %>%
    mutate(result = case_when((rnasep >= cutpoint_rnasep | 
                                 is.na(rnasep)) ~ indeterminate_str,
                              (rnasep < cutpoint_rnasep & 
                                 n1gene < cutpoint_n1gene) ~ positive_str,
                              (rnasep < cutpoint_rnasep & 
                                 (n1gene >= cutpoint_n1gene | 
                                    is.na(n1gene))) ~ negative_str,
                              TRUE ~ fail_str))
  counts <- run_data %>%
    count(barcode)
  results <- left_join(results, counts, by = "barcode")
  
  if (!(negative_str %in% results$barcode)){
    negative_qc <- "FAIL"
    negative_qc_str <- "Could not find negative control!"
  } else if(results$n[results$barcode == negative_str] > 1){
    negative_qc <- "FAIL"
    negative_qc_str <- "More than one negative control found!"
  } else if (results$result[results$barcode == negative_str] == positive_str){
    negative_qc <- "FAIL"
    negative_qc_str <- "Negative control is positive!"
  } else {
    negative_qc <- "PASS"
    negative_qc_str <- "Negative control is negative"
  }
  
  if (!(twist_str %in% results$barcode)){
    twist_qc <- "FAIL"
    twist_qc_str <- "Could not find Twist positive control!"
  } else if(results$n[results$barcode == twist_str] > 1){
    twist_qc <- "FAIL"
    twist_qc_str <- "More than one Twist positive control found!"
  } else if (results$result[results$barcode == twist_str] == negative_str |
             results$result[results$barcode == twist_str] == indeterminate_str){
    twist_qc <- "FAIL"
    twist_qc_str <- "Twist positive control is negative!"
  } else {
    twist_qc <- "PASS"
    twist_qc_str <- "Twist positive is positive"
  }
  
  if ((!idt_str %in% results$barcode)){
    idt_qc <- "FAIL"
    idt_qc_str <- "Could not find IDT positive control!"
  } else if(results$n[results$barcode == idt_str] > 1){
    idt_qc <- "FAIL"
    idt_qc_str <- "More than one idt positive control found!"
  } else if (results$result[results$barcode == idt_str] == negative_str |
             results$result[results$barcode == idt_str] == indeterminate_str){
    idt_qc <- "FAIL"
    idt_qc_str <- "IDT positive control is negative!"
  } else {
    idt_qc <- "PASS"
    idt_qc_str <- "IDT positive control is positive"
  }
  
  return(list(negative = list(qc = negative_qc, qc_str = negative_qc_str),
              twist = list(qc = twist_qc, qc_str = twist_qc_str),
              idt = list(qc = idt_qc, qc_str = idt_qc_str)))
}


create_vdh_template <- function(){
  vdh_template <- tibble(
    Sending_Facility_Name = c(NA),
    Sending_Facility_CLIA = c(NA),
    Message_Control_ID = c(NA),
    PatientID = c(NA),
    SSN = c(NA),
    Last_Name = c(NA),
    First_Name = c(NA),
    Middle_Initial = c(NA),
    Street_Address = c(NA),
    Street_Address_2 = c(NA),
    City = c(NA),
    County_FIPS_Code = c(NA),
    State = c(NA),
    Zip = c(NA),
    Patient_Phone = c(NA),
    Race = c(NA),
    Ethnic_Group = c(NA),
    DOB = c(NA),
    Sex = c(NA),
    Message_Date_Time = c(NA),
    Specimen_ID = c(NA),
    Specimen_Type_Description = c(NA),
    Specimen_Source_Site_Text = c(NA),
    Result_Unit_ID = c(NA),
    Provider_ID = c(NA),
    Provider_Last_Name = c(NA),
    Provider_First_Name = c(NA),
    Ordering_Provider_Addr_1 = c(NA),
    Ordering_Provider_Addr_2 = c(NA),
    Ordering_Provider_City = c(NA),
    Ordering_Provider_State = c(NA),
    Ordering_Provider_Zip = c(NA),
    Ordering_Provider_County_FIPS_code = c(NA),
    Ordering_Provider_Phone = c(NA),
    Ordering_Facility_Name = c(NA),
    Ordering_Facility_Address_1 = c(NA),
    Ordering_Facility_Address_2 = c(NA),
    Ordering_Facility_City = c(NA),
    Ordering_Facility_State = c(NA),
    Ordering_Facility_Zip = c(NA),
    Ordering_Facility_County_FIPS_Code = c(NA),
    Ordering_Facility_Phone = c(NA),
    Observation_Date_Time = c(NA),
    Result_Status = c(NA),
    Specimen_Received_Date = c(NA),
    Order_Code = c(NA),
    Order_Code_Text_Description = c(NA),
    Order_Code_Naming_System = c(NA),
    Result_Value_Type = c(NA),
    Result_Test_Code = c(NA),
    Result_Test_Text_Description = c(NA),
    Result_Test_Naming_System = c(NA),
    Observation_Value = c(NA),
    Observation_Value_Result_Text = c(NA),
    Observation_Value_Result_Naming_System = c(NA),
    Test_Result_Status = c(NA),
    Performing_Lab_ID_Producer_ID = c(NA),
    Performing_Lab_ID_Producer_Text = c(NA),
    Performing_Lab_ID_Producer_Naming_System = c(NA),
    Date_Reported = c(NA),
    Performing_Lab_Street_Address_Line_1 = c(NA),
    Performing_Lab_Street_Address_Line_2 = c(NA),
    Performing_Lab_City = c(NA),
    Performing_Lab_State = c(NA),
    Performing_Lab_Zip = c(NA),
    Performing_Lab_County_FIPS_Code = c(NA),
    Specimen_Type_Identifier = c(NA),
    Specimen_Type_Naming_System = c(NA),
    Date_Test_Ordered = c(NA),
    EUA_Based_Test_Kit_Identification = c(NA),
    Model_Name_Based_Test_Kit_Identification = c(NA),
    Device_Identifier_Based_Test_Kit_Identification = c(NA),
    Model_Name_Based_Instrument_Identification = c(NA),
    Device_Identifier_Based_Instrument_Identification = c(NA),
    Instance_Based_Test_Kit_Identification = c(NA),
    Instance_Based_Instrument_Identification = c(NA),
    Patient_Age_Value = c(NA),
    Patient_Age_Units = c(NA),
    First_Test = c(NA),
    Employed_In_Healthcare = c(NA),
    Symptomatic = c(NA),
    Date_Of_Symptom_Onset = c(NA),
    Hospitalized = c(NA),
    ICU = c(NA),
    Congregate_Care_Setting = c(NA),
    Pregnant = c(NA))
  return(vdh_template)
}