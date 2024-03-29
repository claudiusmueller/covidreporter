compute_run_results <- function(run){
  results <- run %>%
    mutate(result = case_when((rnasep >= cutpoint_rnasep | 
                                 is.na(rnasep)) ~ indeterminate_str,
                              (rnasep < cutpoint_rnasep & 
                                 n1gene < cutpoint_n1gene) ~ positive_str,
                              (rnasep < cutpoint_rnasep & 
                                 (n1gene >= cutpoint_n1gene | 
                                    is.na(n1gene))) ~ negative_str,
                              TRUE ~ fail_str))
  return(results)
}


match_run_results_with_previous <- function(run_results, previous_results){
  previous_results_sel <- previous_results %>%
    filter(barcode %in% run_results$barcode)
  matched_results <- run_results %>%
    bind_rows(previous_results_sel) %>%
    arrange(run_id) %>%
    group_by(barcode) %>%
    mutate(occurance = 1:n())
  matched_results <- matched_results %>%
    pivot_wider(names_from = occurance,
                values_from = c(result, run_id, test_date, rnasep, n1gene)) %>%
    mutate(result_2 = ifelse("result_2" %in% names(.), result_2, NA),
           run_id_2 = ifelse("run_id_2" %in% names(.), run_id_2, NA),
           test_date_2 = ifelse("test_date_2" %in% names(.), test_date_2, NA),
           test_date_2 = as_datetime(test_date_2),
           rnasep_2 = ifelse("rnasep_2" %in% names(.), rnasep_2, NA),
           n1gene_2 = ifelse("n1gene_2" %in% names(.), n1gene_2, NA),
           result_3 = ifelse("result_3" %in% names(.), result_3, NA),
           run_id_3 = ifelse("run_id_3" %in% names(.), run_id_3, NA),
           test_date_3 = ifelse("test_date_3" %in% names(.), test_date_3, NA),
           test_date_3 = as_datetime(test_date_3),
           rnasep_3 = ifelse("rnasep_3" %in% names(.), rnasep_3, NA),
           n1gene_3 = ifelse("n1gene_3" %in% names(.), n1gene_3, NA))
  return(matched_results)
}


set_result_status <- function(matched_results){
  matched_results <- matched_results %>%
    mutate(result_final = 
             case_when(result_1 == negative_str ~ negative_str,
                       result_1 == bad_str ~ bad_str,
                       (result_1 == positive_str & 
                          result_2 == positive_str) ~ positive_str,
                       (result_1 == positive_str &
                          result_2 == bad_str) ~ bad_str,
                       (result_1 == positive_str & 
                          result_2 == negative_str &
                          result_3 == negative_str) ~ negative_str,
                       (result_1 == positive_str &
                          result_2 == negative_str &
                          result_3 == positive_str) ~ positive_str,
                       (result_1 == positive_str &
                          result_2 == negative_str &
                          result_3 == indeterminate_str) ~ bad_str,
                       (result_1 == positive_str &
                          result_2 == negative_str &
                          result_3 == bad_str) ~ bad_str,
                       (result_1 == positive_str &
                          result_2 == indeterminate_str & 
                          result_3 == negative_str) ~ bad_str,
                       (result_1 == positive_str &
                          result_2 == indeterminate_str & 
                          result_3 == positive_str) ~ positive_str,
                       (result_1 == positive_str &
                          result_2 == indeterminate_str & 
                          result_3 == bad_str) ~ bad_str,
                       (result_1 == indeterminate_str & 
                          result_2 == positive_str & 
                          result_3 == negative_str) ~ bad_str,
                       (result_1 == indeterminate_str &
                          result_2 == positive_str &
                          result_3 == positive_str) ~ positive_str,
                       (result_1 == indeterminate_str &
                          result_2 == negative_str) ~ negative_str,
                       (result_1 == indeterminate_str &
                          result_2 == indeterminate_str) ~ bad_str,
                       (result_1 == indeterminate_str &
                          result_2 == bad_str) ~ bad_str,
                       TRUE ~ pending_str),
           test_status = ifelse(result_final == pending_str, preliminary_str,
                                final_str))
  return(matched_results)
}


add_failed_manifest_qc_samples_to_results <- function(run_results,
                                                      sample_manifest){
  failed_samples <- sample_manifest %>%
    filter(str_detect(qc, "([Ff][Aa][Ii][Ll])") == TRUE) %>%
    select(barcode)

  run_results <- run_results %>%
    bind_rows(failed_samples) %>%
    mutate(result = ifelse(barcode %in% failed_samples$barcode, bad_str, 
                           result),
           test_date = ifelse(barcode %in% failed_samples$barcode, test_date[1],
                              test_date),
           test_date = as_datetime(test_date))
  return(run_results)
}


add_failed_accession_qc_samples_to_results <- function(matched_results, 
                                                       sample_accession){
  failed_samples <- sample_accession %>%
    filter(str_detect(qc, "([Ff][Aa][Ii][Ll])") == TRUE) %>%
    select(barcode)
  
  matched_results <- matched_results %>%
    bind_rows(failed_samples) %>%
    mutate(result_1 = ifelse(barcode %in% failed_samples$barcode, bad_str, 
                             result_1),
           test_status = ifelse(barcode %in% failed_samples$barcode, final_str,
                                test_status))
  return(matched_results)
}


select_results_for_report <- function(matched_results, subject_info){
  matched_results <- matched_results %>%
    filter(barcode %in% subject_info$barcode) %>%
    select(barcode, 
           run_id_1, test_date_1, rnasep_1, n1gene_1, result_1,
           run_id_2, test_date_2, rnasep_2, n1gene_2, result_2,
           run_id_3, test_date_3, rnasep_3, n1gene_3, result_3,
           result_final, test_status)
  report_final <- matched_results %>%
    filter(test_status == final_str) %>%
    select(-test_status)
  report_prelim <- matched_results %>%
    filter(test_status == preliminary_str) %>%
    select(-test_status)
  return(list(final = report_final, prelim = report_prelim))
}


format_report_for_shiny_table <- function(report_data){
  report_final <- report_data$final
  report_final <- report_final %>%
    mutate(test_date_1 = format(as_date(test_date_1), "%Y-%m-%d"),
           test_date_2 = format(as_date(test_date_2), "%Y-%m-%d"),
           test_date_3 = format(as_date(test_date_3), "%Y-%m-%d"))
  report_prelim <- report_data$prelim
  report_prelim <- report_prelim %>%
    mutate(test_date_1 = format(as_date(test_date_1), "%Y-%m-%d"),
           test_date_2 = format(as_date(test_date_2), "%Y-%m-%d"),
           test_date_3 = format(as_date(test_date_3), "%Y-%m-%d"))
  return(list(final = report_final, prelim = report_prelim))
    
}


add_subject_info_to_report <- function(report_final, subject_info){
  report_with_subject_info <- report_final %>%
    left_join(subject_info, by = "barcode")
    # mutate(collection_date = ymd_hms(collection_date))
  return(report_with_subject_info)
}


format_condensed_report <- function(report_with_subject_info){
  mymax <- function(x, y, z){
    max(x, y, z, na.rm = TRUE)
  }
  report_condensed <- report_with_subject_info %>%
    mutate(test_date_final = pmap_dbl(list(test_date_1, test_date_2, 
                                           test_date_3), mymax),
           test_date_final = as_datetime(test_date_final),
           test_date_final = as_date(test_date_final),
           dob = as_date(dob)) %>%
    select(barcode, collection_date, g_number, netid,
           first_name, last_name, dob, test_date_final, result_final) %>%
    rename(Barcode = barcode, 
           `Collection Date` = collection_date,
           `G#` = g_number,
           NetID = netid,
           `Last Name` = last_name,
           `First Name` = first_name,
           DOB = dob,
           `Test Date` = test_date_final,
           `Test Result` = result_final)
  return(report_condensed)
}


format_medicat_report <- function(report_with_subject_info){
  report_medicat <- report_with_subject_info %>%
    select(barcode, result_final) %>%
    rename(result = result_final)
  return(report_medicat)
}


add_run_to_covid_db <- function(covid_db, run_results, run_data, subject_info){
  runs <- covid_db$runs
  results <- covid_db$results
  db_subject_info <- covid_db$subject_info
  qc <- "FAIL"
  qc_str <- "Database not updated!"

  barcode_check <- subject_info %>%
    filter(barcode %in% db_subject_info$barcode) %>%
    select(barcode, netid) %>%
    left_join(db_subject_info, by = "barcode") %>%
    select(barcode, netid.x, netid.y) %>%
    mutate(netid_check = ifelse(netid.x == netid.y, TRUE, FALSE))
  
  if (!any(barcode_check$netid_check == FALSE)){
    db_subject_info <- db_subject_info %>%
      bind_rows(subject_info) %>%
      distinct(barcode, netid, g_number, .keep_all = TRUE) %>%
      mutate(dob = as_date(dob))
    
    if (!(run_data$run_id[[1]] %in% runs$run_id)){
      runs <- runs %>%
        bind_rows(list(run_id = run_data$run_id[[1]], 
                       test_date = run_data$test_date[[1]]))
      run_results_db_format <- run_results %>%
        select(-test_date)
      results <- results %>%
        bind_rows(run_results_db_format)
      qc <- "PASS"
      qc_str <- "Test run data added!"
    } else {
      qc <- "QUESTIONABLE"
      qc_str <- "Test run data already in database - no new test data added!"
    }
    
  } else {
    qc <- "FAIL"
    qc_str <- paste("Duplicate barcodes with different netid's in subject information and database: ",
                    paste(barcode_check$barcode[barcode_check$netid_check == FALSE], collapse = ", "))
  }
 
  covid_db <- list(results = results, runs = runs, 
                   subject_info = db_subject_info)
  qc <- list(qc = qc, qc_str = qc_str)
  return(list(covid_db = covid_db, qc = qc))
}


format_vdh_report <- function(vdh_template, report_with_subject_info){
  report_renamed <- report_with_subject_info %>%
    select(barcode, result_final, g_number, last_name, first_name, dob,
           collection_date) %>%
    rename(Specimen_ID = barcode,
           Observation_Value = result_final,
           PatientID = g_number,
           Last_Name = last_name,
           First_Name = first_name,
           DOB = dob,
           Observation_Date_Time = collection_date)
  
  report_time <- format(Sys.time(), "%Y%m%d%k%M%S")
  
  test_str <- "SARS corona virus 2 RNA [Presence] in saliva specimens by extraction-free nucleic acid amplification probe detection"
  
  report_vdh <- full_join(vdh_template, report_renamed) %>%
    drop_na(PatientID) %>%
    mutate(across(everything(), ~ifelse(is.na(.), "", .))) %>%
    mutate(Sending_Facility_Name = "GMU Clinical Proteomics Laboratory",
           Sending_Facility_CLIA = "49D2002076",
           Message_Date_Time = report_time,
           DOB = as_date(DOB),
           DOB = format(DOB, "%Y%m%d"),
           Specimen_Type_Description = "Saliva specimen",
           Ordering_Facility_Name = "George Mason University",
           Ordering_Facility_Address_1 = "4400 University Drive",
           Ordering_Facility_City = "Fairfax",
           Ordering_Facility_State = "VA",
           Ordering_Facility_Zip = 22030,
           Ordering_Facility_Phone = 7039939444,
           Observation_Date_Time = as_datetime(Observation_Date_Time),
           Observation_Date_Time = paste0(format(Observation_Date_Time, 
                                                 "%Y%m%d%k%M%S"), "-0000"),
           Result_Status = "F",
           Specimen_Received_Date = Observation_Date_Time,
           Order_Code = "COTEST PCR",
           Order_Code_Text_Description = test_str,
           Order_Code_Naming_System = "Local",
           Result_Value_Type = "ST",
           Result_Test_Code = "COTEST PCR",
           Result_Test_Text_Description = test_str,
           Result_Test_Naming_System = "Local",
           Observation_Value_Result_Text = Observation_Value,
           Observation_Value_Result_Naming_System = "Local",
           Test_Result_Status = "F",
           Performing_Lab_ID_Producer_Text = "GMU Clinical Proteomics Laboratory",
           Performing_Lab_ID_Producer_Naming_System = "CLIA",
           Date_Reported = paste0(report_time, "-0000"),
           Performing_Lab_Street_Address_Line_1 = "10920 George Mason Circle",
           Performing_Lab_City = "Manassas",
           Performing_Lab_State = "VA",
           Performing_Lab_Zip = 20110,
           Specimen_Type_Identifier = 119342007,
           Specimen_Type_Naming_System = "SCT")
  return(report_vdh)
}