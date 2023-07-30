### Read in Cerner Data
getCernerData <- function(path = "/storage1/fs1/zaydmanm/Active/nspies/BMP/Data/20220713"){

  bmp_raw <-
    list.files(path, pattern = ".tsv", full.names = T) %>%
      map(read_delim, .progress = T) %>%
        bind_rows() %>%
        janitor::clean_names() %>%
        filter(patient_location_facility == "BJH") %>%
        mutate(index_raw = row_number())

  bmp_raw

}

### Deidentify Data ###
deidentifyData <- function(bmp_raw){

  bmp_deidentified <-
    bmp_raw %>%
      mutate(age_at_draw = drawn_dt_tm - patient_dob) %>%
      select(-x1, -patient_name, -patient_dob, -ordering_provider_name,
             -ordering_provider_id, -volume_units, -normal_low, -normal_high,
             -critical_low, -critical_high, -result_units)

  ### Comments may contain PHI and are not yet removed
  bmp_deidentified

}

getLUO <- function(deidentified_data){
  
  luo <- 
    deidentified_data %>% 
      filter(grepl("LUO", task_assay)) %>% 
      dplyr::select(matches("id"), matches("dt_tm"), task_assay, result_val) %>% 
      pivot_wider(everything(), names_from = task_assay, values_from = result_val, values_fn = last) %>%
      group_by(specimen_id) %>%
      fill(matches("LUO"), .direction = "downup") %>%
      distinct(patient_id, drawn_dt_tm, performed_dt_tm, verified_dt_tm, .keep_all = T)
  
  luo
  
}

keepNumeric <- function(bmp_deidentified = deidentified_data){

  bmp_deidentified <-
    bmp_deidentified %>%
      dplyr::filter(task_assay %in% c("Sodium", "Chloride", "Potassium Plas", "CO2 Totl", "BUN", "Creatinine", "Calcium", "Glucose"))

  non_numeric <- bmp_deidentified %>% dplyr::filter(is.na(as.numeric(result_val))) %>% select(task_assay, result_val) %>% group_by(task_assay) %>% count(result_val) %>% arrange(desc(n))
  write_delim(non_numeric, "../Results/QC/non_numeric_results.tsv", delim = "\t")

  bmp_deidentified_numeric <- bmp_deidentified %>% mutate(result_val = as.numeric(str_replace_all(result_val, "<|>|,", "")))

  output <- bmp_deidentified_numeric %>% dplyr::filter(!is.na(result_val))

  output

}

rectangularizeResults <- function(bmp_numeric) {

  duplicates <-
    bmp_numeric %>%
     dplyr::group_by(specimen_id, orig_order_dt_tm, drawn_dt_tm, received_dt_tm, performed_dt_tm, verified_dt_tm, task_assay, long_text) %>%
      dplyr::summarise(n = dplyr::n(), max = max(result_val), min = min(result_val), .groups = "drop") %>%
      dplyr::filter(n > 1L)

  duplicates <- duplicates %>% mutate(range = max - min)

  write_delim(duplicates, "../Results/QC/duplicated_bmps.tsv")

  bmp_wide <-
    bmp_numeric %>%
      arrange(drawn_dt_tm, performed_dt_tm, verified_dt_tm) %>%
      pivot_wider(
        id_cols = c(specimen_id),
        names_from = task_assay,
        values_from = result_val,
        values_fill = NA,
        values_fn = list) %>%
      janitor::clean_names()

  output <-
    bmp_wide %>%
      ungroup() %>%
      arrange(patient_id, drawn_dt_tm) %>%
      mutate(label = "Patient",
             index = row_number(),
             anion_gap = sodium - chloride - co2_totl)

  output

  }

fillBMP <- function(bmp_wide){

    print("Filling missing results from same specimen ID...")
    bmp_filled <-
      bmp_wide %>%
        select(index, matches("_id"), matches("dt_tm"), any_of(lab_strings_bmp)) %>%
          group_by(specimen_id) %>%
          fill(any_of(lab_strings_bmp), .direction = "downup")

  bmp_filled_distinct_verified <-
    bmp_filled %>%
      distinct(patient_id, drawn_dt_tm, performed_dt_tm, verified_dt_tm, .keep_all = T) %>%
        ungroup() %>%
          filter(if_all(.cols = all_of(lab_strings_bmp), function(x) !is.na(x)))
  qsave(bmp_filled_distinct_verified, "_targets/objects/bmp_filled_distinct_verified")

  bmp_filled_distinct_performed <-
    bmp_filled %>%
      distinct(patient_id, drawn_dt_tm, performed_dt_tm, .keep_all = T) %>%
        ungroup() %>%
          filter(if_all(.cols = all_of(lab_strings_bmp), function(x) !is.na(x)))

  qsave(bmp_filled_distinct_performed, "_targets/objects/bmp_filled_distinct_performed")

  bmp_filled_distinct_drawn <-
    bmp_filled %>%
      distinct(patient_id, drawn_dt_tm, .keep_all = T) %>%
        ungroup() %>%
          filter(if_all(.cols = all_of(lab_strings_bmp), function(x) !is.na(x)))

  qsave(bmp_filled_distinct_drawn, "_targets/objects/bmp_filled_distinct_drawn")

  bmp_filled_distinct_verified

}

### Calculates Delta Matrices ###
calculateDeltas <- function(bmp_no_NA, time_diff_hours = 28) {
  
  print("Adding Patient ID and Drawn_DT_TM from metadata...")
  input <-
    bmp_no_NA %>% 
      arrange(patient_id, drawn_dt_tm, verified_dt_tm) %>% 
      mutate(hours_since_last_draw = 
               ifelse(patient_id == dplyr::lag(patient_id), as.numeric(drawn_dt_tm - dplyr::lag(drawn_dt_tm)) / 60 / 60, NA))
  
  print("Adding `recent_prior` flag to calculate deltas if samples are <24 hours (time_diff_hours) apart...")
  input <-
    input %>%
      mutate(recent_prior = if_else(hours_since_last_draw <= time_diff_hours, T, F, missing = F))
  
  print("Calculating Delta Matrix...")
  bmp_deltas <-
    as_tibble(input[,lab_strings_bmp] - dplyr::lag(input[,lab_strings_bmp])) %>%
      set_names(paste0(lab_strings_bmp, "_delta")) %>%
      bind_cols(index = input$index, recent_prior = input$recent_prior, .)
  
  library(DescTools)
  bmp_deltas[which(!bmp_deltas$recent_prior), which(names(bmp_deltas) %like% "%delta%")] <- NA

  print("Saving Delta Summaries to Results/Deltas/...")
  delta_summaries <-
    bmp_deltas %>%
      filter(recent_prior) %>%
      select(matches(paste0(lab_strings_bmp, collapse = "|"))) %>%
      lapply(., summary)
  write.table(capture.output(delta_summaries), file = here::here("../Results/Deltas/delta_summaries.txt"), sep = "\t", row.names = F, quote = F)
  
  print("Saving Delta Distributions Figures to Figures/Deltas/...")
  ggplot(
    data = bmp_deltas %>% filter(recent_prior) %>% pivot_longer(cols = matches(paste0(lab_strings_bmp, collapse = "|")), names_to = "Analyte", values_to = "Delta"),
    aes(x = Delta)) +
    geom_density() +
    facet_wrap( ~ Analyte, scales = "free") +
    theme(axis.text = element_text())
  ggsave(here::here("../Results/Deltas/delta_distributions.png"), width = 11, height = 11, dpi = 600)

  bmp_deltas
  
}



