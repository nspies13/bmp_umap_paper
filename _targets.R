library(targets)
library(tarchetypes) 
library(tidyverse)
library(tidymodels)
library(pins)

tar_option_set(
  packages = c("tidyverse", "tidymodels", "qs", "future", "pins", "foreach", "parallelly", "embed"), format = "qs", deployment = "worker", 
  memory = "transient",  garbage_collection = T, storage = "worker", retrieval = "worker", 
  envir = globalenv(), tidy_eval = T, error = "continue", iteration = "list")

future::plan(future.callr::callr, workers = 8)
tar_source()
tar_source("../Helpers/SetGlobals.R")
tar_source("../Retrospective/R/retrospective_contamination.R")
tar_source("../Simulation/R/simulate_contamination.R")
tar_source("../Preprocessing/R/preprocessing.R")
tar_source("../Supervised/R/assess_models.R")
tar_source("../Supervised/R/assess_predictions.R")

##### Prep ML Inputs #####
load_inputs <- list(
  
  # Read preprocessed data from LIS extract board.
  tar_target(preprocessed_bmp_inputs, list(list(data = data_board %>% pin_read("BJH_bmps_with_error_flags") %>% 
                                                  drop_na(any_of(lab_strings_bmp)) %>% 
                                                  rename(patient_id = person_id, specimen_id = container_id), label = "inpatient"))),
  tar_target(train_cohort, train_cohorts),
  tar_target(panels_tar, panels),
  tar_target(contam_sim, map2(fluids, fluid_names, ~makeFullContaminationTibble(preprocessed_bmp_inputs[["data"]], mix_ratios = rep(seq(0.01, 1, by = 0.01), each = 1000), fluid = .x, fluid_name = .y)) %>% bind_rows(), pattern = map(preprocessed_bmp_inputs)),

  # Make train-test splits
  tar_target(split, group_initial_split(preprocessed_bmp_inputs[["data"]], group = patient_id, prop = 0.8), pattern = map(preprocessed_bmp_inputs)),
  tar_target(train_input, list(train_split = group_initial_split(training(split), group = patient_id, prop = 0.8), train_cohort = train_cohort), pattern = map(split, train_cohort)),
  tar_target(train, list(train = makeUnsupervisedInput(training(train_input[["train_split"]]), fluids, class_balance = 0.5), train_cohort = train_input[["train_cohort"]]), pattern = map(train_input)),
  tar_target(validation_set, list(validation = testing(train_input[["train_split"]]), train_cohort = train_input[["train_cohort"]]), pattern = map(train_input)),
  tar_target(test, testing(split), pattern = map(split))
  
)

##### Build Unsupervised Models ######
build_unsup_pipelines <- tar_map(
  
  values = 
    expand_grid(
      panel_map = panels, 
      train_cohort_map = train_cohorts, 
      mode = c("unsupervised", "semisupervised"),
      class_balance = c(0.5),
      num_comp = c(2),
      metric = c("cosine"), # Prior grid: cosine, manhattan
      neighbors = c(50, 100), # Prior grid: 25, 50, 100 
      bandwidth = c(100), # Prior grid: 2, 10, 100
      min_dist = c(0), # Prior grid: 0, 0.01
      init = c("spca"), # Prior grid: spectral, spca, agspectral 
      set_op_mix_ratio = c(1), # Prior grid: 0, 0.5, 1 
      dens_scale = c(1), # Prior grid: 0, 0.5, 1 
      fluid_list_map = c("TopTen")) %>% 
    mutate(label = paste(panel_map, train_cohort_map, mode, class_balance, num_comp, fluid_list_map, metric, neighbors, bandwidth, min_dist, init, set_op_mix_ratio, dens_scale, sep = "_")),
  names = "label",
  
  ### Build Models
  tar_target(unsup_wf_sets, makeUnsupRecipes(train[[1]][["train"]], num_comp = num_comp, mode = mode, metric = metric, neighbors = neighbors, bandwidth = bandwidth, min_dist = min_dist, init = init, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale))

)

apply_models <- list(
  
  ### Load Final Model
  tar_target(model, qs::qread("_targets/objects/unsup_wf_sets_BMP_BJH_unsupervised_0.5_2_TopTen_cosine_100_100_0_spca_1_1") %>% pluck("results_umap")),
  
  ### Apply Models to Real Train/Test Set #####
  tar_target(validation_umap, applyRecipes(model, validation_set[[1]][["validation"]])),
  tar_target(test_umap, applyRecipes(model, test[[1]])),
  tar_target(sim_umap, applyRecipes(model, contam_sim[[1]])),

  ### Likelihood Ratio Grids
  tar_target(likelihood_ratio_umap_grid, makeEnrichmentGrid(validation_umap, sim_umap)),
  tar_target(validation_ratio_umap, assignEnrichmentScores(validation_umap, grid = likelihood_ratio_umap_grid)),
  tar_target(umap_threshold, getRatioThreshold(validation_ratio_umap)),
  
  tar_target(sim_ratio_umap, assignEnrichmentScores(sim_umap, grid = likelihood_ratio_umap_grid) %>% mutate(prediction = factor(ifelse(likelihood_ratio >= umap_threshold, 1, 0)))),
  tar_target(test_ratio_umap, assignEnrichmentScores(test_umap, grid = likelihood_ratio_umap_grid) %>% mutate(prediction = factor(ifelse(likelihood_ratio >= umap_threshold, 1, 0)))),

  tar_target(assess_umap, assessUnsupPerformance(test = test_ratio_umap, sim = sim_ratio_umap, label = "final_clinchem")),
  tar_target(manual_review, compareToManualReview())
  
)

list(load_inputs, apply_models)
  