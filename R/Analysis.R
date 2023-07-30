### Simulate IV Fluid Contamination
simulateContamination <- function(input, n = 1000) {

  set.seed(12345)

  output <- input %>% mutate(label = "Patient")

  percent_fluid <- rep(seq(0.01, 0.99, by = 0.01), each = n)

  library(foreach)
  foreach(i = 1:length(fluid_names)) %do% {

    tmp <-
      output %>%
      select(index, all_of(lab_strings_bmp), label) %>%
      slice_sample(n = length(percent_fluid), replace = T) %>%
      mutate(percent_fluid = 0)

    tmp2 <-
      tmp[, lab_strings_bmp] * (1 - percent_fluid) +
      as.data.frame(lapply(as.data.frame(matrix(
        rep(fluids[[i]][lab_strings_bmp], each = length(percent_fluid)), nrow = length(percent_fluid))),
        as.numeric)) * percent_fluid

    tmp2 <-
      tmp2 %>%
      mutate(across(c("sodium", "chloride", "co2_totl", "bun", "glucose"), ~round(.))) %>%
      mutate(across(c("potassium_plas", "calcium"), ~round(., 1))) %>%
      mutate(creatinine = round(creatinine, 2))

    tmp2 <-
      tmp2 %>%
      mutate(
        anion_gap = sodium - chloride - co2_totl,
        label = fluid_names[i],
        percent_fluid = percent_fluid,
        index = tmp$index
      )

    assign(paste0(fluid_names[i], "_mix"), tmp2)
  }

  bmp_mixed_all <- bind_rows(NS_mix, D5NS_mix, LR_mix, D5LR_mix)

  bmp_mixed_all <-
    bmp_mixed_all %>%
    ungroup() %>%
    mutate(anion_gap = sodium - chloride - co2_totl)

  bmp_mixed_all

}

simulateContaminationRow <- function(input, mix_ratio, fluid){
  
  fluid = eval(parse(text = fluid))
  
  input[names(fluid)] * (1 - mix_ratio) + fluid * mix_ratio
  
}

### Define DKA
defineDKA <- function(abg_raw, bmp_no_NA_verified){

  tar_load(abg_raw);tar_load(bmp_no_NA_verified)

  low_pH <- abg_raw %>% dplyr::filter(grepl("pH", task_assay) & result_val < 7.3)
  high_ketones <- abg_raw %>% dplyr::filter((grepl("Butyrate", task_assay) & result_val > 0.4) | (grepl("Ketone", task_assay) & result_val > 1))

  dka <-
    bind_rows(low_pH, high_ketones) %>%
    mutate(task_assay =
             case_when(
               task_assay == "pH Art gPOC" ~ "ph_art",
               task_assay == "pH Art" ~ "ph_art",
               task_assay == "Ketones, Ur" ~ "ketones",
               task_assay == "Ketones Ur Ql" ~ "ketones",
               task_assay == "B OH Butyrate" ~ "ketones")) %>%
    pivot_wider(id_cols = matches("dt_tm|id"), names_from = "task_assay", values_from = "result_val") %>%
    unnest(cols = c(ph_art, ketones)) %>%
    distinct()

  dka_filled <-
    dka %>%
    arrange(patient_id, drawn_dt_tm) %>%
    group_by(patient_id) %>%
    fill(everything(), .direction = "down")

  combined <-
    bind_rows(dka_filled, bmp_no_NA_verified) %>%
    arrange(patient_id, drawn_dt_tm) %>%
    group_by(patient_id) %>%
    fill(everything(), .direction = "down")

  dka_bmp <-
    combined %>%
    dplyr::filter(
      glucose > 250 &
        (ph_art < 7.3 | ketones > 1) &
        co2_totl < 22 &
        anion_gap > 12)

  dka_bmp

}

### Make Input For DKA vs Dextrose UMAP Model
makeHyperglycemiaInput <- function(dka_results, contam_sim){
  
 dka_input <- 
   dka_results %>% ungroup() %>% 
    select(all_of(lab_strings_bmp)) %>% distinct() %>% 
    mutate(label = "DKA")
 
 getSimilarGlucose <- function(glucose_map, sim){
   
   contam_sim %>% 
     filter(label %in% c("D5NS", "D5LR") & near(glucose_map, glucose, tol = 10)) %>% 
     slice_sample(n = 1)
   
 }
  
 matched_d5 <- bind_rows(lapply(dka_input$glucose, function(x) getSimilarGlucose(x, contam_sim)))
 
 input <- bind_rows(dka_input, matched_d5 %>% mutate(label = "D5")) %>% select(-percent_fluid)
 
 ggplot(input, aes(glucose, fill = label)) + geom_boxplot(alpha = 0.5) + theme(legend.title = element_blank())
 ggsave(here::here("../Figures/DKA_D5_matching_glucose_boxplots.png"), width = 11, height = 11, dpi = 300)
 
 input
 
}

### Identify Retrospective Contamination
identifyContamination <- function(bmp_results_no_NA, bmp_deltas, metadata) {

  print("Joining required data...")
  output <-
    metadata %>% select(index, patient_id, drawn_dt_tm) %>%
    left_join(., bmp_results_no_NA %>% select(index, any_of(lab_strings_bmp))) %>%
    left_join(., bmp_deltas %>% select(-matches("delta_cv"))) %>%
    arrange(patient_id, drawn_dt_tm) %>%
    group_by(patient_id)

  ##### Define Moving Functions #####
  lastMoveTowards <- function(data, col, fluid) {
    data$compare <- fluid[[col]]

    output_list <-
      ifelse(
        (lag(data[, col]) > data[, "compare"] &
           data[, paste0(col, "_delta")] < 0 &
           abs(data[, paste0(col, "_delta")]) <= abs(lag(data[, col]) - data[, "compare"]) &
           data[, "recent_prior"])
        |
          (lag(data[, col]) < data[, "compare"] &
             data[, paste0(col, "_delta")] > 0 &
             abs(data[, paste0(col, "_delta")]) <= abs(lag(data[, col]) - data[, "compare"]) &
             data[, "recent_prior"]),
        T, F)

    output_list[, 1]
  }
  nextMoveAway <- function(data, col, fluid) {
    data$compare <- fluid[[col]]

    output_list <-
      ifelse(
        (data[, col] > data[, "compare"] &
           lead(data[, paste0(col, "_delta")] > 0) &
           lead(data[, "recent_prior"]))
        |
          (data[, col] < data[, "compare"] &
             lead(data[, paste0(col, "_delta")]) < 0 &
             lead(data[, "recent_prior"]))
        |
          (data[, col] == data[, "compare"]) &
          (lead(data[, col]) != data[, col]) &
          lead(data[, "recent_prior"]),
        T, F)

    output_list[, 1]
  }

  ##### Apply Dummy Data To Test Workflow #####
  library(future);library(doFuture);library(parallel);

  test_data <- bind_rows(NS * 0.5, NS, NS * 0.5, NS * 0.5 + 1, NS, NS * 0.5 + 1)
  test_data <- tibble(test_data - lag(test_data)) %>% set_names(paste0(lab_strings_bmp, "_delta")) %>%
    bind_cols(recent_prior = T, test_data, .)

  moved_to <-
    foreach(fluid = fluids) %:%
    foreach(col = lab_strings_bmp, .combine = "bind_cols", .export = "lastMoveTowards") %dopar% {
      lastMoveTowards(test_data, col, fluid)
    }
  names(moved_to) <- fluid_names

  moves_away <-
    foreach(fluid = fluids) %:%
    foreach(col = lab_strings_bmp, .combine = "bind_cols", .export = "nextMoveAway") %dopar% {
      nextMoveAway(test_data, col, fluid)
    }
  names(moves_away) <- fluid_names

  moved_to <- moved_to %>% map(., ~set_names(., paste0(lab_strings_bmp, "_towards")))
  moves_away <- moves_away %>% map(., ~set_names(., paste0(lab_strings_bmp, "_away")))

  print("Counting number of labs that moved_to or moves_away per row...")
  to_sums <- lapply(moved_to, rowSums)
  away_sums <- lapply(moves_away, rowSums)

  contam_list_by_threshold <-
    foreach(fluid = fluid_names) %:%
    foreach(threshold = 1:9) %dopar% {
      to_sums[[fluid]] >= threshold & away_sums[[fluid]] >= threshold
    }

  contam_tables_by_threshold <-
    foreach(fluid = 1:length(fluids)) %dopar% {
      lapply(contam_list_by_threshold[[fluid]], table)
    }

  contam_long <-
    tibble(
      fluid = rep(fluid_names, each = 9),
      threshold = rep(1:9, times = 4),
      contaminated = unlist(lapply(flatten(contam_tables_by_threshold), function(x) x[2]))
    )

  ##### Apply Functions #####
  print("Assigning _moved_towards_ columns...")
  moved_to <-
    foreach(fluid = fluids) %:%
    foreach(col = lab_strings_bmp, .combine = "bind_cols", .export = "lastMoveTowards") %dopar% {
      lastMoveTowards(output, col, fluid)
    }
  names(moved_to) <- fluid_names

  print("Assigning _moves_away_ columns...")
  moves_away <-
    foreach(fluid = fluids) %:%
    foreach(col = lab_strings_bmp, .combine = "bind_cols", .export = "nextMoveAway") %dopar% {
      nextMoveAway(output, col, fluid)
    }
  names(moves_away) <- fluid_names

  moved_to <- moved_to %>% map(., ~set_names(., paste0(lab_strings_bmp, "_towards")))
  moves_away <- moves_away %>% map(., ~set_names(., paste0(lab_strings_bmp, "_away")))

  ##### Print Function Sanity Check #####
  function_check <-
    output %>%
    select(patient_id, recent_prior, sodium, sodium_delta) %>%
    bind_cols(., moved_to_NS = moved_to$NS$sodium_towards, moves_away_NS = moves_away$NS$sodium_away)
  write_delim(function_check, here("../Results/QC/MoveTo_MovesAway_sanity_check.tsv"), delim = "\t")

  print("Counting number of labs that moved_to or moves_away per row...")
  to_sums <- lapply(moved_to, rowSums)
  away_sums <- lapply(moves_away, rowSums)

  combined <- map2(moved_to, moves_away, ~bind_cols(.x, .y))

  ##### Assign Contam Columns By Rules #####
  NS_contam <-
    combined$NS %>%
    bind_cols(index = bmp_results_no_NA$index, recent_prior = bmp_deltas$recent_prior,  .) %>%
    mutate(NS =
             ifelse(
               recent_prior &
                 sodium_towards & sodium_away &
                 chloride_towards & chloride_away &
                 calcium_towards & calcium_away &
                 potassium_plas_towards & potassium_plas_away,
               T, F))

  D5NS_contam <-
    combined$D5NS %>%
    bind_cols(index = bmp_results_no_NA$index, recent_prior = bmp_deltas$recent_prior, ., glucose = bmp_results_no_NA$glucose) %>%
    mutate(D5NS =
             ifelse(
               recent_prior &
                 sodium_towards & sodium_away &
                 chloride_towards & chloride_away &
                 calcium_towards & calcium_away &
                 potassium_plas_towards & potassium_plas_away &
                 glucose_towards & glucose_away & glucose > 300,
               T, F))

  LR_contam <-
    combined$LR %>%
    bind_cols(index = bmp_results_no_NA$index, recent_prior = bmp_deltas$recent_prior,  .) %>%
    mutate(LR =
             ifelse(
               recent_prior &
                 sodium_towards & sodium_away &
                 chloride_towards & chloride_away &
                 calcium_towards & calcium_away,
               T, F))

  D5LR_contam <-
    combined$D5LR %>%
    bind_cols(index = bmp_results_no_NA$index, recent_prior = bmp_deltas$recent_prior,  ., glucose = bmp_results_no_NA$glucose) %>%
    mutate(D5LR =
             ifelse(
               recent_prior &
                 sodium_towards & sodium_away &
                 chloride_towards & chloride_away &
                 calcium_towards & calcium_away &
                 glucose_towards & glucose_away & glucose > 300,
               T, F))

  ##### Aggregate All Moved To/Away #####
  contam_list_by_threshold <-
    foreach(fluid = fluid_names) %:%
    foreach(threshold = 1:9) %dopar% {
      to_sums[[fluid]] >= threshold & away_sums[[fluid]] >= threshold
    }

  contam_tables_by_threshold <-
    foreach(fluid = 1:length(fluids)) %dopar% {
      lapply(contam_list_by_threshold[[fluid]], table)
    }

  contam_long <-
    tibble(
      fluid = rep(fluid_names, each = 9),
      threshold = rep(1:9, times = 4),
      contaminated = unlist(lapply(flatten(contam_tables_by_threshold), function(x) x[2]))
    )

  contam_long$contaminated <-
    contam_long$contaminated %>% replace_na(0)

  contam_final_cols <-
    bind_cols(index = NS_contam$index,
              NS = NS_contam$NS,
              D5NS = D5NS_contam$D5NS,
              LR = LR_contam$LR,
              D5LR = D5LR_contam$D5LR)

  output <- full_join(output, contam_final_cols, by = "index")

  ##### Print Retrospectively Identified Contaminated Results With Surrounding BMPs to Files #####
  print("Writing retrospectively identified contamination candidates to /NS_contam_spot_check.txt for manual review...")
  NS_indices <- which(output$NS)
  NS_indices <- sort(unique(c(NS_indices, NS_indices - 1, NS_indices + 1)))
  NS_contam_results <- bind_cols(output[NS_indices, ],
                                 num_toward = to_sums$NS[NS_indices],
                                 num_away = away_sums$NS[NS_indices])
  write.table(NS_contam_results, here("../Results/NS_contamination_retrospective.txt"), quote = F, sep = "\t", row.names = F)

  D5NS_indices <- which(output$D5NS)
  D5NS_indices <- sort(unique(c(D5NS_indices, D5NS_indices - 1, D5NS_indices + 1)))
  D5NS_contam_results <- bind_cols(output[D5NS_indices, ],
                                   num_toward = to_sums$D5NS[D5NS_indices],
                                   num_away = away_sums$D5NS[D5NS_indices])
  write.table(D5NS_contam_results, here("../Results/D5NS_contamination_retrospective.txt"), quote = F, sep = "\t", row.names = F)

  LR_indices <- which(output$LR)
  LR_indices <- sort(unique(c(LR_indices, LR_indices - 1, LR_indices + 1)))
  LR_contam_results <- bind_cols(output[LR_indices, ],
                                 num_toward = to_sums$LR[LR_indices],
                                 num_away = away_sums$LR[LR_indices])
  write.table(LR_contam_results, here("../Results/LR_contamination_retrospective.txt"), quote = F, sep = "\t", row.names = F)

  D5LR_indices <- which(output$D5LR)
  D5LR_indices <- sort(unique(c(D5LR_indices, D5LR_indices - 1, D5LR_indices + 1)))
  D5LR_contam_results <- bind_cols(output[D5LR_indices, ],
                                   num_toward = to_sums$D5LR[D5LR_indices],
                                   num_away = away_sums$D5LR[D5LR_indices])
  write.table(D5LR_contam_results, here("../Results/D5LR_contamination_retrospective.txt"), quote = F, sep = "\t", row.names = F)


  ##### Reassign Label Column #####
  print("Assigning ID'ed contamination candidates new labels...")
  output <-
    output %>%
    ungroup() %>%
    mutate(label = if_else(
      (NS | LR | D5NS | D5LR),
      "Contamination", "Patient"))

  qsave(contam_list_by_threshold, "_targets/objects/contam_lists")
  qsave(contam_tables_by_threshold, "_targets/objects/contam_tables")
  qsave(contam_final_cols, "_targets/objects/contam_final_cols")

  output <- output %>%
    select(index, label, all_of(fluid_names))

  output

}

### Make Enrichment Scores 
calculateEnrichmentScores <- function(results = unsup_train_kde, sim = sim_results_kde){
  
  results_bin <- 
    results %>% mutate(bin_x = round(UMAP1, 1), bin_y = round(UMAP2, 1)) %>% 
      select(bin_x, bin_y, kde_patients = kde_results_umap) %>% group_by(bin_x, bin_y) %>% 
      summarize(mean_patient = mean(kde_patients, na.rm = T))
  
  sim_bin <- 
    sim %>% mutate(bin_x = round(UMAP1, 1), bin_y = round(UMAP2, 1)) %>% 
      select(bin_x, bin_y, kde_sim = kde_results_umap) %>% group_by(bin_x, bin_y) %>% 
    summarize(mean_sim = mean(kde_sim, na.rm = T))
  
  enrichment_scores <- 
    full_join(results_bin, sim_bin, by = c("bin_x", "bin_y")) %>%
      mutate(enrichment = mean_sim / mean_patient) 
  
  fill_min_e <- enrichment_scores[which(ntile(enrichment_scores$enrichment, n = 1000) == 10)[[1]], "enrichment"][[1]]
  fill_max_e <- enrichment_scores[which(ntile(enrichment_scores$enrichment, n = 1000) == 990)[[1]], "enrichment"][[1]]
  
  enrichment_scores <- 
    enrichment_scores %>% 
      mutate(enrichment = 
               case_when(
                 is.na(mean_sim) & !is.na(mean_patient) ~ fill_min_e,
                 !is.na(mean_sim) & is.na(mean_patient) ~ fill_max_e,
                 T ~ enrichment
               ))

  enrichment_scores
  
}

### Assign Enrichment to Test Set
assignEnrichmentScores <- function(enrichment_scores, unsup_results_kde){
  
  unsup_results_enrich <- 
    unsup_results_kde %>% 
      mutate(bin_x = round(UMAP1, digits = 1), bin_y = round(UMAP2, digits = 1)) %>%
      left_join(enrichment_scores) %>% 
      mutate(enrich_ntile = ntile(enrichment, 1000))
  
  unsup_results_enrich
  
}

printResultsToReview <- function(test_results_enrich = test_results_enrich, bmp_deltas = bmp_deltas, bmp_no_NA = bmp_no_NA, n = 100){
  
  sample <- test_results_enrich %>% arrange(desc(enrichment)) %>% distinct(sodium, chloride, potassium_plas, glucose, calcium, bun, creatinine, co2_totl, .keep_all = T) %>% slice_head(n = n)
  
  sample
  
  library(foreach)
  full_review <- 
    foreach::foreach(index_keep = sample$index, .combine = bind_rows) %dopar% {
    
    tmp <- bmp_no_NA %>% filter(index == index_keep) %>% select(patient_id, drawn_dt_tm)
    
    same_patient <- bmp_no_NA %>% filter(patient_id == tmp$patient_id) %>% distinct(drawn_dt_tm, sodium, chloride, potassium_plas, glucose, calcium, bun, creatinine, co2_totl)
    
    row <- which(same_patient$drawn_dt_tm == tmp$drawn_dt_tm)[[1]]
    if(row == 1){ rows = c(1,2,3) }
    if(row == 2){ rows = c(1,2,3,4) }
    if(row > 2){ rows = seq(row - 2, row + 2) }

    full = same_patient %>% slice(rows) %>% left_join(., bmp_no_NA) %>% distinct(drawn_dt_tm, sodium, chloride, potassium_plas, glucose, calcium, bun, creatinine, co2_totl, .keep_all = T)
    
    full
    
    }
  
  sample <- sample %>% mutate(review = T)

  full_review <- full_review %>% left_join(sample %>% select(index, review))
    
  write_delim(
    full_review %>% arrange(patient_id, drawn_dt_tm) %>% select(patient_id, drawn_dt_tm, all_of(lab_strings_bmp), review),
    file = here::here("../Results/enrichment_to_review.txt"), delim = '\t')
  
  
}

assessReviewedFlags <- function(){
  
  reviewed_enrich <- read_delim(here::here("../Results/enrich_sample_reviewed.txt"))
  
  reviewed_enrich <- reviewed_enrich %>% mutate(contam = ifelse(grepl("REAL|DKA", NS), F, T))
  reviewed_enrich %>% filter(!is.na(NS) & !contam) %>% summary(anion_gap)
  
  reviewed_proj <- read_delim(here::here(("../Results/projection_reviewed.txt")))
  
  reviewed_proj %>% count(ns)

  }

### Get Projection
getNSprojection <- function(unsup_results_embed){
  
  indexes <- unsup_results_embed %>% filter(index != 1057279 & UMAP1 > 4.085 & UMAP2 > -2.5 & UMAP2 < -1.5) %>% select(index) %>% pluck(1)

  library(foreach)
  review <- 
    foreach(index_keep = indexes, .combine = bind_rows) %dopar% {
    
    tmp <- bmp_no_NA %>% filter(index == index_keep) %>% select(patient_id, drawn_dt_tm)
    
    same_patient <- bmp_no_NA %>% filter(patient_id == tmp$patient_id) %>% distinct(drawn_dt_tm, sodium, chloride, potassium_plas, glucose, calcium, bun, creatinine, co2_totl)
    
    row <- which(same_patient$drawn_dt_tm == tmp$drawn_dt_tm)[[1]]
    if(row == 1){ rows = c(1,2,3) }
    if(row == 2){ rows = c(1,2,3,4) }
    if(row > 2){ rows = seq(row - 2, row + 2) }
    
    full = same_patient %>% slice(rows) %>% left_join(., bmp_no_NA) %>% distinct(drawn_dt_tm, sodium, chloride, potassium_plas, glucose, calcium, bun, creatinine, co2_totl, .keep_all = T)
    
    full
    
  }
  
  review <- review %>% mutate(review = ifelse(index %in% indexes, T, F))
  
  review <- 
    review %>% arrange(patient_id, drawn_dt_tm) %>% 
      mutate(same_after = ifelse(patient_id == dplyr::lead(patient_id), T , F), 
             same_prior = ifelse(patient_id == dplyr::lag(patient_id), T, F))
    
  write_delim(review, file = here::here("../Results/projection_to_review.txt"), delim = '\t')
  
}

### Assign Comments to Flags
assignComments <- function(bmp_numeric){
  
  comments <- 
    bmp_numeric %>% 
      dplyr::select(matches("id"), matches("dt_tm"), long_text) %>%
      mutate(
        code_comment = grepl("code 7|code blue|code ", long_text, ignore.case = T),
        contamination_comment = grepl("IV |i\\.v\\.|contam", long_text, ignore.case = T),
        mislabel_comment = grepl("mislabel", long_text, ignore.case = T),
        unlikely_comment = grepl("unlikely results", long_text, ignore.case = T),
        glycolysis_comment = grepl("glycolysis", long_text, ignore.case = T)
      )
  
  comments
  
}

### Color by Clinical Variables
assignClinicalVariables <- function(bmp_numeric){
  
  clinical <- 
    bmp_numeric %>% 
      dplyr::select(patient_id, specimen_id, drawn_dt_tm, performed_dt_tm, verified_dt_tm,
             age_at_draw, patient_sex, patient_location_facility, patient_location_nursing,
             collection_priority, collection_method, collector_role, long_text)
  
  clinical <- 
    clinical %>%
      mutate(
        icu = grepl("ICU|CTI|I1|I0", patient_location_nursing),
        emergency = grepl("EM|ED", patient_location_nursing))
  
  clinical

}