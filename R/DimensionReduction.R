# Build Tidymodels Unsupervised Recipe
buildUnsupModels <- function(input, neighbors = 25, bandwidth = 1, connectivity = 1, metric = "cosine", min_dist = 0.001, dens_scale = 1, mix_ratio = 0.1, init = init) {

  library(tidymodels)
  library(embed)
  set.seed(12345)

  input <- input %>% dplyr::select(any_of(lab_strings_bmp))

  unsup_recipe <-
    recipe(~ ., data = input) %>%
    embed::step_umap(
      all_numeric_predictors(),
      neighbors = neighbors,
      num_comp = 2,
      min_dist = min_dist,
      seed = c(12, 23),
      keep_original_cols = T,
      options = list(
        bandwidth = bandwidth,
        local_connectivity = connectivity,
        metric = metric,
        min_dist = min_dist,
        dens_scale = dens_scale,
        init = init,
        set_op_mix_ratio = mix_ratio,
        n_threads = parallel::detectCores() - 2,
        n_sgd_threads = parallel::detectCores() - 2,
        fast_sgd = F,
        ret_model = T,
        verbose = F)
    ) %>%
    prep()

  unsup_recipe

  library(vetiver);library(bundle);library(pins)

  bundled <- bundle::bundle(unsup_recipe)

  bundled

}

# Plot Embedding
plotUMAP <- function(model){

  model <- bundle::unbundle(model)

  model %>% pluck("template") %>% slice_sample(prop = 0.01) %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(alpha = 0.01, stroke = 0, shape = 20)

}

# Apply Recipe
applyUnsupRecipes <- function(recipe, raw_input){

  library(tidymodels)

  recipe <- bundle::unbundle(recipe)

  cols <- recipe %>% pluck("template") %>% names

  new_data <- raw_input %>% dplyr::select(index, any_of(cols)) %>% na.omit

  output <- recipe %>% prep() %>% bake(new_data = new_data)

  output <- bind_cols(new_data %>% dplyr::select(-any_of(names(output))), output)

}

# Apply Kernel Density Estimate
getKDE <- function(input, cols, name = "kde", eps = 1){

  kde <- bind_cols(input, name = dbscan::pointdensity(input %>% select(cols), eps = eps, type = "density"))

  kde

}

# Build Hyperglycemia Model
buildHyperglycemiaModel <- function(hyperglycemia_input){
  
  library(tidymodels)
  library(embed)
  set.seed(12345)
  
  input <- hyperglycemia_input %>% dplyr::select(any_of(lab_strings_bmp))
  
  unsup_recipe <-
    recipe(~ ., data = input) %>%
    embed::step_umap(
      all_numeric_predictors(),
      neighbors = 35,
      num_comp = 2,
      min_dist = 1,
      seed = c(12, 23),
      keep_original_cols = T,
      options = list(
        n_threads = parallel::detectCores() - 2,
        n_sgd_threads = parallel::detectCores() - 2,
        dens_scale = 1,
        set_op_mix_ratio = 1,
        bandwidth = 100,
        local_connectivity = 2,
        init = "spca",
        metric = "cosine",
        fast_sgd = F,
        ret_model = T,
        verbose = T)
    ) %>%
    prep()
  
  unsup_recipe
  
  library(vetiver);library(bundle);library(pins)
  
  bundled <- bundle::bundle(unsup_recipe)
  
  bundled
  
}