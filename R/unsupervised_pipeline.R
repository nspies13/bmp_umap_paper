makeUnsupervisedInput <- function(input = train_input[[1]][["data"]], fluids_to_include = fluids, class_balance = 0.05){
  
  mixes <- seq(0.01, 1, by = 0.01)
  mix_ratios <- rep(mixes, each = 10000)
  
  contam_sim_all <-
    c(
      map2(fluids_to_include[which(grepl("D5|SW|hyper", names(fluids_to_include)))], 
           names(fluids_to_include[which(grepl("D5|SW|hyper", names(fluids_to_include)))]), 
           ~makeFullContaminationTibble(input, 
                                        mix_ratios = mix_ratios/2, 
                                        fluid = .x, 
                                        fluid_name = .y)),
      map2(fluids_to_include[which(!grepl("D5", names(fluids_to_include)))], 
           names(fluids_to_include[which(!grepl("D5", names(fluids_to_include)))]), 
           ~makeFullContaminationTibble(input, 
                                        mix_ratios = mix_ratios, 
                                        fluid = .x, 
                                        fluid_name = .y)))
  
  rows_per_group <- (nrow(input)*class_balance)/(length(mixes)*length(fluids_to_include))
  
  contam <- contam_sim_all %>% bind_rows() %>% group_by(label, mix_ratio) %>% slice_sample(n = floor(rows_per_group)) %>% ungroup()
  
  output <- bind_rows(contam, input) %>% ungroup()
  
  output

}

makeUnsupRecipes <- function(train, mode = c("unsupervised"), metric = "cosine", neighbors = 25, min_dist = 0, bandwidth = 1, init = "spca", set_op_mix_ratio = 0, dens_scale = 0, num_comp = 2){
  
  if(mode == "unsupervised"){
    
    base_results_rec <-     
      recipe(train %>% drop_na(any_of(!!lab_strings))) %>%
      update_role(everything(), new_role = "metadata") %>%
      update_role_requirements("metadata", bake = F) %>%
      update_role(any_of(c(!!!lab_strings)), new_role = "predictor") %>%
      step_normalize(all_predictors())
    
    results_umap_rec <- 
      base_results_rec %>%
      embed::step_umap(any_of(c(!!!lab_strings_bmp_no_gap)), keep_original_cols = F, metric = metric, num_comp = num_comp, neighbors = neighbors, min_dist = min_dist, 
                       options = list(bandwidth = bandwidth, local_connectivity = 100, init = init, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale,
                                      fast_sgd = F, approx_pow = F, n_threads = 62, n_sgd_threads = 62, tmpdir = "/scratch1/fs1/zaydmanm/nspies/umap/", ret_model = T, verbose = F)) %>% 
      prep() %>%
      butcher::butcher(verbose = T) %>% 
      bundle::bundle()
    
  }
  
  if(mode == "semisupervised"){
    
    base_results_rec <-     
      recipe(train %>% drop_na(any_of(!!lab_strings))) %>%
      update_role(everything(), new_role = "metadata") %>%
      update_role_requirements("metadata", bake = F) %>%
      update_role(any_of(c(!!!lab_strings)), new_role = "predictor") %>%
      step_normalize(all_predictors())
    
    results_umap_rec <- 
      base_results_rec %>%
      embed::step_umap(any_of(c(!!!lab_strings_bmp_no_gap)), keep_original_cols = F, outcome = vars(label), metric = metric, num_comp = num_comp, neighbors = neighbors, min_dist = min_dist, target_metric = metric, target_weight = 0.25,
                       options = list(bandwidth = bandwidth, local_connectivity = 100, init = init, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale,
                                      fast_sgd = F, approx_pow = F, n_threads = 62, n_sgd_threads = 62, tmpdir = "/scratch1/fs1/zaydmanm/nspies/umap/", ret_model = T, verbose = F)) %>% 
      prep() %>%
      butcher::butcher(verbose = T) %>% 
      bundle::bundle()
    
  }
  
  if(mode == "supervised"){
      
      base_results_rec <-     
        recipe(train %>% drop_na(any_of(!!lab_strings)) %>% mutate(label = ifelse(is.na(label), "Real", label), mix_ratio = ifelse(is.na(mix_ratio), 0, mix_ratio))) %>%
        update_role(everything(), new_role = "metadata") %>%
        update_role_requirements("metadata", bake = F) %>%
        update_role(any_of(c(!!!lab_strings)), new_role = "predictor") %>%
        step_normalize(all_predictors())
      
      results_umap_rec <- 
        base_results_rec %>%
        embed::step_umap(any_of(c(!!!lab_strings_bmp_no_gap)), keep_original_cols = F, outcome = vars(label), metric = metric, num_comp = num_comp, neighbors = neighbors, min_dist = min_dist, 
                         options = list(bandwidth = bandwidth, local_connectivity = 100, init = init, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale,
                                        fast_sgd = F, approx_pow = F, n_threads = 62, n_sgd_threads = 62, tmpdir = "/scratch1/fs1/zaydmanm/nspies/umap/", ret_model = T, verbose = F)) %>% 
        prep() %>%
        butcher::butcher(verbose = T) %>% 
        bundle::bundle()
      
    }
  
  results_pca_rec <- 
    base_results_rec %>%
      step_pca(any_of(c(!!!lab_strings_bmp_no_gap)), num_comp = 8, keep_original_cols = F, options = list(center = T, scale. = T)) %>%
      prep() %>%
      bundle::bundle()
  
  list("results_umap" = results_umap_rec, "results_pca" = results_pca_rec)
  
}

makeUnsupDeltaRecipes <- function(train, metric = "cosine", neighbors = 25, bandwidth = 1, init = "spca", set_op_mix_ratio = 0, dens_scale = 0){
  
  base_deltas_rec <-
    recipe(train %>% drop_na(matches("delta"))) %>%
      update_role(everything(), new_role = "metadata") %>%
      update_role_requirements("metadata", bake = F) %>%
      update_role(matches("delta_prior"), new_role = "predictor") %>%
      step_normalize(matches("delta_prior")) 

  deltas_umap_rec <-
    base_deltas_rec %>%
      embed::step_umap(matches("delta_prior"), keep_original_cols = F, metric = metric, num_comp = 2, neighbors = neighbors, min_dist = 0, 
                       options = list(bandwidth = bandwidth, local_connectivity = 100, n_trees = 200, search_k = 400 * neighbors, init = init, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale,
                                      fast_sgd = T, n_threads = 62, n_sgd_threads = 62, tmpdir = "/scratch1/fs1/zaydmanm/nspies/umap/", ret_model = T, verbose = F)) %>%
      prep() %>%
      bundle::bundle()

  deltas_pca_rec <-
    base_deltas_rec %>%
      step_pca(matches("delta_prior"), num_comp = 8, keep_original_cols = F, options = list(center = T, scale. = T)) %>%
      prep() %>%
      bundle::bundle()

  list("deltas_umap" = deltas_umap_rec, "deltas_pca" = deltas_pca_rec)
  
}

makeUnsupWorkflowSets <- function(train, metric = "cosine", neighbors = 25, bandwidth = 1, init = "spca", set_op_mix_ratio = 0, dens_scale = 0){
  
  train <- train %>% mutate(target = factor(ifelse(label == "Patient", 0, 1)), mix_ratio = ifelse(target == 0, 0, mix_ratio)) %>% drop_na(any_of(!!lab_strings), matches("delta"))
  
  model = list(rf = rand_forest(mode = "classification", engine = "ranger", mtry = 2, trees = 100)) 
  recipes = makeUnsupRecipes(train)
  
  wf_set <- workflow_set(recipes, model, cross = T) 
  
  prepped <- prepRecipes(wf_set)
  
  list(wf_set = wf_set, prepped = prepped)
  
}

applyRecipes <- function(recipe, raw_input){
  
  library(tidymodels)
  
  recipe <- bundle::unbundle(recipe)

  output <- recipe %>% bake(new_data = raw_input)
  
  output <- bind_cols(raw_input, output %>% select(matches("PC|UMAP")))
  
  output
  
}

plotEmbeddings <- function(input, type = "UMAP"){
  
  if(type == "UMAP"){
  theme_set(theme_ns)
   p <- 
     ggplot() + 
      geom_point(data = input %>% filter(!is.na(label)), aes(UMAP1, UMAP2, color = label), alpha = 0.25, shape = 20, stroke = 0, size = 3) +
      geom_point(data = input %>% filter(is.na(label)), aes(UMAP1, UMAP2), color = "grey50", alpha = 0.02, shape = 20, stroke = 0, size = 3) + 
      scale_color_manual(values = color_map_global) + 
      guides(color = guide_legend(title = "Fluid", reverse = T, override.aes = list(size = 8, stroke = 1, alpha = 1))) + 
      theme(axis.text = element_blank(), legend.key = element_rect(fill = NA, color = NA))
    
  }
  
  if(type == "PCA"){
  theme_set(theme_ns)
    p <- 
      ggplot() + 
        geom_point(data = input %>% filter(!is.na(label)), aes(PC1, PC2, color = label), alpha = 0.25, shape = 20, stroke = 0, size = 3) +
        geom_point(data = input %>% filter(is.na(label)), aes(PC1, PC2), color = "grey50", alpha = 0.02, shape = 20, stroke = 0, size = 3) + 
        scale_color_manual(values = color_map_global) + 
        guides(color = guide_legend(title = "Fluid", reverse = T, override.aes = list(size = 8, stroke = 1, alpha = 1))) + 
        theme(axis.text = element_blank(), legend.key = element_rect(fill = NA, color = NA))
  
  }

  p
  
}

makeEmbeddingPlots <- function(input = prepped, label = ""){

  theme_set(theme_ns)
  plots <- map2(input %>% map("template") %>% map(~slice_sample(.x, prop = 0.25)), rep(c("UMAP", "PCA"), times = 1), ~plotEmbeddings(.x, type = .y))
  
  library(ggpubr)
  ggarrange(plotlist = plots, ncol = 2, nrow = 1, common.legend = T, legend = "bottom") +
    annotate("text", label = "Results", x = 0.05, y = 0.98, hjust = 0, size = 4, fontface = "bold") 
#    annotate("text", label = "Deltas", x = 0.05, y = 0.52, hjust = 0, size = 4, fontface = "bold")
  ggsave(paste0("../../Figures/Unsupervised/TuningPlots/", label, "_tuning_plot.pdf"), width = 8.5, height = 4)
  
  
}

### Make Enrichment Scores 
makeEnrichmentGrid <- function(real_embed = tar_read(validation_umap_BMP_BJH_0.5_NS_LR_D5NS_D5LR_cosine_50_10_spca_1_1), sim = tar_read(sim_umap_BMP_BJH_0.5_NS_LR_D5NS_D5LR_cosine_50_10_spca_1_1)){
  
  library(KernSmooth)
  
  kde_input <- 
    bind_rows(
      D5 = sim %>% filter(mix_ratio > 0.05 & mix_ratio < 0.25 & grepl("D5|hyper", label)), 
      nonD5 = sim %>% filter(mix_ratio > 0.25 & !grepl("D5|hyper", label)))

  uncontam_real <- real_embed %>% filter(!contam_comment & !unlikely_comment & !code_comment & !credited_comment) %>% slice_sample(n = nrow(kde_input), replace = T)
  
  min_x <- min(c(kde_input$UMAP1, uncontam_real$UMAP1))
  max_x <- max(c(kde_input$UMAP1, uncontam_real$UMAP1))
  
  min_y <- min(c(kde_input$UMAP2, uncontam_real$UMAP2))
  max_y <- max(c(kde_input$UMAP2, uncontam_real$UMAP2))
  
  kde_uncontam <- 
    bkde2D(
      uncontam_real %>% select(UMAP1, UMAP2) %>% as.matrix(), 
      bandwidth = 0.2, 
      gridsize = c(10000, 10000),
      range.x = list(c(min_x, max_x), c(min_y, max_y)))
  
  kde_contam <- 
    bkde2D(
      kde_input %>% select(UMAP1, UMAP2) %>% as.matrix(), 
      bandwidth = 0.2, 
      gridsize = c(10000, 10000),
      range.x = list(c(min_x, max_x), c(min_y, max_y)))
  
  grid <- kde_contam
  grid$fhat <- (kde_contam$fhat + 1e-20)/(kde_uncontam$fhat + 1e-20)
  
  grid
  
}

### Assign Enrichment to Test Set
assignEnrichmentScores <- function(input, grid = likelihood_ratio_grid){
  
  bin_x <- cut(input$UMAP1, breaks = grid$x1, right = F) %>% as.numeric
  bin_y <- cut(input$UMAP2, breaks = grid$x2, right = F) %>% as.numeric
  
  enrichment <- map2(bin_x, bin_y, ~grid[["fhat"]][.x, .y]) %>% unlist()
  
  input$likelihood_ratio <- enrichment
  
  input
  
}

### Plot Enrichment
plotLikelihoodRatio <- function(input){
  
  cap1 <- quantile(input %>% pluck("likelihood_ratio"), probs = c(0.9995), na.rm = T)[[1]]
  
  test_overlay <- input %>% arrange(likelihood_ratio) %>% slice_tail(prop = 0.1) %>% mutate(likelihood_ratio = ifelse(likelihood_ratio > cap1, cap1, likelihood_ratio))
  
  ggplot(test_overlay, aes(UMAP1, UMAP2, fill = likelihood_ratio)) + 
    geom_point(color = "grey50", size = 3, stroke = 0.25, shape = 21, alpha = 1) +
    scico::scale_fill_scico(direction = -1, palette = "grayC", trans = scales::log1p_trans(), limits = c(0, cap1), breaks = c(0, cap1), labels = c("Lowest", "Highest")) +
    guides(fill = guide_colorbar(title = "Likelihood Ratio", direction = "horizontal", title.position = "top", ticks = F, title.theme = element_text(size = 10, face = "bold", hjust = 0.5), label.theme = element_text(size = 8))) +
    ggtitle("Patient Embedding")
  
}

### Get Threshold That Would Match Current Alarm Rate from Validation Data
getRatioThreshold <- function(input = tar_read(validation_ratio_umap)){
  
  threshold <- 
    input %>% 
      pluck("likelihood_ratio") %>%
        quantile(., probs = c(0.997), na.rm = T) %>% 
        pluck(1)
  
  library(geomtextpath)
  ggplot(input %>% mutate(flagged = ifelse(unlikely_comment | contam_comment, "Flagged", "Real")), aes(likelihood_ratio, fill = flagged)) + 
    geom_density(color = NA, alpha = 0.5, show.legend = F) +
    geom_textdensity(aes(label = flagged), hjust = 0.4, vjust = -0.1, size = 4, fontface = "bold", show.legend = F) + 
    geom_vline(xintercept = threshold, linetype = "dashed") +
    annotate("text", x = threshold + 10, y = 0.7, label = "Threshold", angle = -90, color = "black", alpha = 1) + 
    scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c(0.1, 1, 10, 100, 1000)) +
    scale_fill_manual(values = c("grey20", "grey80")) + 
    xlab("Likelihood Ratio") + ylab("Relative Proportion") +
    theme(axis.text.y.left = element_blank())
  ggsave("../../Figures/Unsupervised/TestPerformance/final_clinchem_UMAP/threshold_density_current_flags.svg", width = 6, height = 3)
  
  
  threshold
  
}

### Aggregate Metrics 
aggregateTestCommentMetrics <- function(input = tar_read(test_ratio_umap_BMP_BJH_0.5_NS_LR_D5NS_D5LR_cosine_50_10_spca_0_0), threshold = 1,
                                        train_cohort = "BJH", panel = "BMP", fluids = "NS_LR_D5NS_D5LR", algorithm = "UMAP", class_balance = 0.5, 
                                        init = "spca", metric = "cosine", neighbors = 50, bandwidth = 10, set_op_mix_ratio = 1, dens_scale = 1){
  
  preds <- 
    input %>% 
      mutate(target = as_factor(as.numeric(contam_comment)),
             pred = as_factor(ifelse(likelihood_ratio >= threshold, 1, 0), levels = c(0, 1)))
  
  conf_mat <- 
    preds %>% 
      conf_mat(target, pred)
  
  metrics <- 
    preds %>%
      metrics_binary(target, estimate = pred, event_level = "second")
  
  list(preds = preds, conf_mat = conf_mat, metrics = metrics, 
       train_cohort = train_cohort, panel = panel, fluids = fluids, 
       algorithm = algorithm, class_balance = class_balance,
       init = init, metric = metric, neighbors = neighbors, 
       bandwidth = bandwidth, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale)
  
}

### Aggregate Assessment 
assessUnsupPerformance <- function(test = tar_read(test_ratio_umap_BMP_BJH_unsupervised_0.5_2_TopTen_cosine_100_100_0_spca_1_1), sim = tar_read(sim_ratio_umap_BMP_BJH_unsupervised_0.5_2_TopTen_cosine_100_100_0_spca_1_1), grid = tar_read(likelihood_ratio_umap_grid_BMP_BJH_TopTen_cosine_25_10_0_spca_1_1), label = "review_final", algorithm = "UMAP"){
                         
  out_dir = paste0("../../Figures/Unsupervised/TestPerformance/", label, "_", algorithm) 
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  
  ##### Plot Embeddings on Real and Sim Data #####
  theme_quad <- 
    theme(title = element_text(size = 14, face = "bold.italic", hjust = 0.5), axis.line = element_line(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.key = element_blank(), legend.background = element_blank(), axis.title.x.bottom = element_text(size = 12, face = "bold.italic"), axis.title.y.left = element_blank(), axis.text = element_blank())
  
  theme_set(theme_quad)
  
  x1 <- quantile(test$UMAP1, probs = c(0.001, 0.999))
  y1 <- quantile(test$UMAP2, probs = c(0.001, 0.999))
  cap1 <- quantile(test$likelihood_ratio, probs = c(0.99995), na.rm = T)[[1]]
  
  library(KernSmooth)
  kde <- 
    bkde2D(
      test %>% select(UMAP1, UMAP2) %>% as.matrix(), 
      bandwidth = 0.2, 
      gridsize = c(10000, 10000)) 
  
  bin_x <- cut(test$UMAP1, breaks = kde$x1, right = F) %>% as.numeric
  bin_y <- cut(test$UMAP2, breaks = kde$x2, right = F) %>% as.numeric
  
  test$kde <- map2(bin_x, bin_y, ~kde[["fhat"]][.x, .y]) %>% unlist()
  
  test_overlay <- test
  #%>% mutate(likelihood_ratio = ifelse(likelihood_ratio > cap1, cap1, likelihood_ratio)) %>% slice_sample(prop = 0.05, weight_by = 1/kde) %>% mutate(alpha = (105 - ntile(kde, n = 100))/100)
  
  x2 <- quantile(sim %>% filter(mix_ratio < 0.95) %>% pluck("UMAP1"), probs = c(0.001, 0.999))
  y2 <- quantile(sim %>% filter(mix_ratio < 0.95) %>% pluck("UMAP2"), probs = c(0.001, 0.999))
  cap2 <- quantile(sim %>% filter(mix_ratio < 0.75) %>% pluck("likelihood_ratio"), probs = c(0.95), na.rm = T)[[1]]
  
  gg_test_ratio <- 
    ggplot(test, aes(UMAP1, UMAP2, color = kde)) + 
    scattermore::geom_scattermore() +
    scico::scale_color_scico(palette = "oslo", breaks = c(0.05, 1), end = 0.95, labels = c("Lowest", "Highest")) +
    guides(alpha = F, color = guide_colorbar(title = "Patient Density", direction = "horizontal", title.position = "top", ticks = F, title.theme = element_text(size = 8, face = "bold", hjust = 0.5), label.theme = element_text(size = 8))) +
    theme(plot.background = element_blank(), legend.position = c(0.25, 0.25), axis.title.y.left = element_text(size = 12, face = "bold.italic"), axis.title.x.bottom = element_text(size = 12, face = "bold.italic"))
  ggsave(paste0(out_dir, "/test_set_likelihood_ratio_embedding.png"), gg_test_ratio, width = 8.5, height = 6, dpi = 600)
  ggsave(paste0(out_dir, "/test_set_likelihood_ratio_embedding.pdf"), gg_test_ratio, width = 8.5, height = 6)
  ggsave(paste0(out_dir, "/test_set_likelihood_ratio_embedding.svg"), gg_test_ratio, width = 8.5, height = 6)
  
  
  sim_plot <- sim %>% group_by(label, mix_ratio) %>% mutate(likelihood_ratio = ifelse(likelihood_ratio > cap2, cap2, likelihood_ratio))
  
  sim_fluid <- 
  ggplot() + 
    scattermore::geom_scattermore(data = sim_plot, aes(UMAP1, UMAP2, color = label)) +
    xlim(x2[[1]] - 0.25, x2[[2]] + 0.2) +
    ylim(y2[[1]] - 0.25, y2[[2]] + 0.1) +
    scale_color_manual(values = color_map_global) +
    guides(color = guide_legend(title = "Fluid", direction = "horizontal", title.position = "top", reverse = T,
                                    ticks = F, title.theme = element_text(size = 10, face = "bold", hjust = 0.5), nrow = 2,  
                                    label.theme = element_text(size = 8, face = "bold"), override.aes = list(size = 8, alpha = 1, stroke = 1, shape = 21))) +
        theme(axis.title.y.left = element_text(size = 12, face = "bold.italic"),legend.position = "none")
  ggsave(paste0(out_dir, "/sim_set_fluid_embedding_grid.pdf"), sim_fluid, width = 8.5, height = 6)
  ggsave(paste0(out_dir, "/sim_set_fluid_embedding_grid.svg"), sim_fluid, width = 8.5, height = 6)
  ggsave(paste0(out_dir, "/sim_set_fluid_embedding_grid.png"), sim_fluid, width = 8.5, height = 6, dpi = 600)
              
  sim_mix_ratio <-           
      ggplot() + 
        scattermore::geom_scattermore(data = sim_plot, aes(UMAP1, UMAP2, color = mix_ratio)) +
        xlim(x2[[1]] - 0.25, x2[[2]] + 0.2) +
        ylim(y2[[1]] - 0.25, y2[[2]] + 0.1) +
        scico::scale_color_scico(direction = -1, end = 0.95, palette = "davos", breaks = c(0, 1), limits = c(0, 1), labels = c(0, 1)) +
        guides(color = guide_colorbar(title = "Mixture Ratio", direction = "horizontal", title.position = "top", ticks = F, 
                                      title.theme = element_text(size = 10, face = "bold", hjust = 0.5),
                                      label.theme = element_text(size = 8, face = "bold", hjust = 0))) +
        theme(legend.position = c(0.25, 0.55))
              ggsave(paste0(out_dir, "/sim_set_mix_ratio_embedding_grid.pdf"), sim_mix_ratio, width = 8.5, height = 6)
              ggsave(paste0(out_dir, "/sim_set_mix_ratio_embedding_grid.svg"), sim_mix_ratio, width = 8.5, height = 6)
              ggsave(paste0(out_dir, "/sim_set_mix_ratio_embedding_grid.png"), sim_mix_ratio, width = 8.5, height = 6, dpi = 600)
   
  sim_likelihood <-            
      ggplot() + 
        scattermore::geom_scattermore(data = sim_plot, aes(UMAP1, UMAP2, color = likelihood_ratio)) +    
        xlim(x2[[1]] - 0.25, x2[[2]] + 0.2) +
        ylim(y2[[1]] - 0.25, y2[[2]] + 0.1) +
        scico::scale_color_scico(direction = 1, begin = 0.15, palette = "bilbao", trans = scales::log1p_trans(), oob = scales::squish, limits = c(0, cap1*2), breaks = c(0, cap1*2), labels = c("Lowest", "Highest")) + 
        guides(color = guide_colorbar(title = "Likelihood Ratio", direction = "horizontal", title.position = "top", ticks = F, 
                                                title.theme = element_text(size = 10, face = "bold", hjust = 0.5),
                                                label.theme = element_text(size = 8, face = "bold", hjust = 0))) + 
        theme(legend.position = c(0.25, 0.55))
    ggsave(paste0(out_dir, "/sim_set_likelihood_ratio_embedding_grid.pdf"), sim_likelihood, width = 8.5, height = 6)
    ggsave(paste0(out_dir, "/sim_set_likelihood_ratio_embedding_grid.svg"), sim_likelihood, width = 8.5, height = 6)
    ggsave(paste0(out_dir, "/sim_set_likelihood_ratio_embedding_grid.png"), sim_likelihood, width = 8.5, height = 6, dpi = 600)
    
  library(ggpubr)  
  gg_sim <- ggarrange(plotlist = list(sim_fluid, sim_mix_ratio, sim_likelihood), nrow = 1, ncol = 3, align = "hv")
  ggsave(paste0(out_dir, "/sim_set.pdf"), gg_sim, width = 13, height = 4)
  ggsave(paste0(out_dir, "/sim_set.svg"), gg_sim, width = 13, height = 4)
  ggsave(paste0(out_dir, "/sim_set.png"), gg_sim, width = 13, height = 4, dpi = 600)
  
      
  ##### Sensitivity Analysis #####
  theme_set(theme_ns)
  library(geomtextpath)
    
  gg_sens <- 
    sim %>% 
      group_by(label, mix_ratio) %>% 
      count(prediction) %>%
      mutate(prop = n/sum(n)) %>%
      filter(prediction == 1) %>%
      ggplot(aes(mix_ratio, prop, color = label)) +
        geom_textpath(aes(label = label), linewidth = 2, fontface = "bold", hjust = 0.25, vjust = -0.2) +
        scale_color_manual(values = color_map_global) +
        xlab("Mixture Ratio") + ylab("Sensitivity") +
        ggtitle("Sensitivity in Simulated Contamination") + 
        theme(legend.position = "none")
  ggsave(paste0(out_dir, "/sensitivity_analysis.pdf"), gg_sens, width = 8.5, height = 4)
  ggsave(paste0(out_dir, "/sensitivity_analysis.svg"), gg_sens, width = 8.5, height = 4)
  ggsave(paste0(out_dir, "/sensitivity_analysis.png"), gg_sens, width = 8.5, height = 4, dpi = 600)
                         
  sens_tables <- 
     sim %>% 
      group_by(label, mix_ratio) %>% 
      count(prediction) %>%
      mutate(prop = n/sum(n)) %>%
      filter(prediction == 1 & mix_ratio %in% c(0.1, 0.25, 0.5)) %>% 
      select(label, mix_ratio, prop)
  
  ##### Comparison to Reported Technologist Comments #####
  test <- 
    test %>%
      mutate(truth = factor(as.numeric(contam_comment|unlikely_comment), levels = c(0, 1))) 
  
  conf_mat <- 
    test %>%
      conf_mat(truth, prediction)
  
  metrics <- 
    test %>%
      metrics_binary(truth, estimate = prediction, event_level = "second")
  
  roc <- 
    test %>% roc_curve(truth, likelihood_ratio, event_level = "second") 
  
  gg_roc <- 
    ggplot(test %>% roc_curve(truth, likelihood_ratio, event_level = "second"), aes(1-specificity, sensitivity)) + 
      geom_path() +
      coord_equal(ratio = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
      annotate("text", x = 0.85, y = 0.15, label = paste0("AUC = ", sprintf("%.2f", round(roc_auc(test, truth, likelihood_ratio, event_level = "second") %>% pluck(".estimate"), digits = 2))), size = 4, fontface = "bold") +
      xlim(0,1) + ylim(0,1) + 
      xlab("1 - Specificity") + ylab("Sensitivity") +
      ggtitle("ROC of Embedding Likelihood Ratio", subtitle = "Gold Standard: Technologist Contamination Flag")
  
  gg_pr <- 
    ggplot(test %>% pr_curve(truth, likelihood_ratio, event_level = "second") , aes(recall, precision)) + 
      geom_path() +
      coord_equal(ratio = 1) +
      annotate("text", x = 0.15, y = 0.15, label = paste0("AUC = ", sprintf("%.2f", round(pr_auc(test, truth, likelihood_ratio, event_level = "second") %>% pluck(".estimate"), digits = 2))), size = 4, fontface = "bold") +
      xlim(0,1) + ylim(0,1) + 
      xlab("Sensitivity") + ylab("Positive Predictive Value") +
      ggtitle("PR Curve of Embedding Likelihood Ratio", subtitle = "Gold Standard: Technologist Contamination Flag")
  
  gg_curves <- ggpubr::ggarrange(gg_roc, gg_pr)
  ggsave(paste0(out_dir, "/test_set_roc_pr_curves.pdf"), gg_curves, width = 13, height = 6)
  
  set.seed(12345)
  discrepancies_to_review <- 
    test %>% 
      filter(truth != prediction & drawn_dt_tm > "2020-01-01") %>%
      drop_na(matches("prior|post")) %>%
      slice_sample(n = 200, by = prediction)
  write_csv(discrepancies_to_review, paste0(out_dir, "/discrepancy_sample.csv"))
  
  concordance_to_review <-
    test %>% 
      filter(truth == prediction & drawn_dt_tm > "2020-01-01") %>%
      drop_na(matches("prior|post")) %>%
      slice_sample(n = 50, by = prediction)
  write_csv(discrepancies_to_review, paste0(out_dir, "/concordant_sample.csv"))
  
  out <- bind_rows(discrepancies_to_review, concordance_to_review)
  
  one <- 
    out[sample(1:nrow(out)),] %>% 
      reformatReviews() 
  
  two <- 
    one %>%
      left_join(tmp %>% select(epic_id, drawn_dt_tm, lab_strings_bmp), by = c("epic_id", lab_strings_bmp)) %>%
      mutate(drawn_dt_tm = ifelse(time_point == "Current", drawn_dt_tm.x, drawn_dt_tm.y)) %>%
      select(epic_id, drawn_dt_tm, time_point, any_of(lab_strings_bmp), fluids_infusing_at_draw, relevant_diagnoses, final_prediction, confidence, minutes_of_assessment, comments, contaminating_fluid) %>% 
      mutate(drawn_dt_tm = as_datetime(drawn_dt_tm) + 5 * 60 * 60)
    write_csv(two, paste0(out_dir, "/final_subset_samples_to_review.csv"))
    
  consecutive_pos <- 
    test %>%
      arrange(drawn_dt_tm) %>% 
      filter(prediction == 1) %>%
      slice_tail(n = 100) %>%
      reformatReviews() %>%
      left_join(test %>% select(any_of(lab_strings_bmp), matches("comment"), likelihood_ratio, prediction) %>% mutate(flagged = unlikely_comment|contam_comment))
  write_csv(consecutive_pos, paste0(out_dir, "/final_positives.csv"))
  
  
  consecutive <- 
    test %>%
      arrange(drawn_dt_tm) %>%
      slice_tail(n = 1000) %>%
      reformatReviews() %>%
      left_join(test %>% select(any_of(lab_strings_bmp), matches("comment"), likelihood_ratio, prediction) %>% mutate(flagged = unlikely_comment|contam_comment))
  write_csv(consecutive, paste0(out_dir, "/final_consecutives.csv"))
  
  
  # 
  # list(test = test, sim = sim, label = label, grid = grid, threshold = threshold, algorithm = "UMAP", 
  #      init = init, metric = metric, neighbors = neighbors, bandwidth = bandwidth, set_op_mix_ratio = set_op_mix_ratio, dens_scale = dens_scale,
  #      patient_likelihood_plot = gg_test_ratio, sim_likelihood_plot = gg_sim_ratio, roc_pr_curves = gg_curves, sensitivity_analysis_plot = gg_sens,
  #      sensitivity_tables = sens_tables, conf_mat = conf_mat, metrics = metrics, discrepancies_to_review = discrepancies_to_review)
                  
}

makeUnsupVetiverObject <- function(recipe = tar_read(unsup_wf_sets_BMP_BJH_0.5_NS_LR_D5NS_D5LR_cosine_50_10_spca_1_1) %>% pluck("results_umap"), metadata = tar_read(assess_umap_BMP_BJH_0.5_NS_LR_D5NS_D5LR_cosine_50_10_spca_1_1)){
  
  recipe = bundle::unbundle(recipe) %>% butcher::butcher(verbose = T)
  
  library(vetiver)
  proto <- metadata[["test"]] %>% select(any_of(lab_strings)) %>% slice_sample(n = 100)
  
  v <- 
    vetiver::vetiver_model(
      model = recipe, 
      model_name = "UMAP_ClinChem_Final_lean",
      save_prototype = proto,
      metadata = metadata[which(!grepl("plot", names(metadata)))],
      versioned = T
    )
  
  model_board %>%
    vetiver_pin_write(
      vetiver_model = v,
      check_renv = T
    )
  
  v
  
}

plotTechnologistFlaggedResults <- function(input = tar_read(test_ratio_umap)){
  
  library(ggpointdensity)
  
  input %>% 
    filter(contam_comment) %>%
    ggplot(aes(UMAP1, UMAP2)) + 
    geom_pointdensity() + 
    scale_color_viridis_c()
    
  
}

reformatReviews <- function(review = discrepancies_to_review){
  
  library(foreach)
  foreach(i = 1:nrow(review)) %do% {
    
    list(
      bind_rows(review[i,] %>% select(epic_id, drawn_dt_tm, matches("prior") & !matches("delta")) %>% mutate(time_point = "Prior", .after = drawn_dt_tm) %>% setNames(gsub("_prior", "", names(.))) %>% mutate(fluids_infusing_at_draw = " ", relevant_diagnoses = " ", final_prediction = " ", confidence = " ", minutes_of_assessment = " ", comments = " "),
                review[i,] %>% select(epic_id, drawn_dt_tm, lab_strings_bmp) %>% mutate(time_point = "Current", .after = drawn_dt_tm) %>% mutate(fluids_infusing_at_draw = "-", relevant_diagnoses = "-", final_prediction = "Contaminated/Real", contaminating_fluid = "Pick a fluid", confidence = "[0-1]", minutes_of_assessment = "-", comments = "-"), 
                review[i,] %>% select(epic_id, drawn_dt_tm, matches("post") & !matches("delta")) %>% setNames(gsub("_post", "", names(.))) %>% mutate(time_point = "Post", .after = drawn_dt_tm) %>% mutate(fluids_infusing_at_draw = " ", relevant_diagnoses = " ", final_prediction = " ", confidence = " ", minutes_of_assessment = " ", comments = " "),
                as_data_frame(rep("-", times = 2))
                )
                
    )
    
  } %>% bind_rows()
    
}

assessOperationalImpact <- function(test = tar_read(test_ratio_umap_BMP_BJH_unsupervised_0.5_2_TopTen_cosine_100_100_0_spca_1_1)){
  
  times_raw <- qs::qread("../Preprocessing/_targets/objects/bmp_comment_flags")
  
  times <-  
    times_raw %>% 
      rename(specimen_id = container_id) %>% 
      filter(specimen_id %in% test$specimen_id) %>%
      group_by(specimen_id) %>% 
      reframe(drawn_dt_tm = min(drawn_dt_tm), first_performed = min(perform_dt_tm), last_resulted = max(perform_result_updt_dt_tm))
  
  times_flagged <- 
    times %>% 
      left_join(test %>% transmute(specimen_id, prediction, contam_comment, unlikely_comment, target = contam_comment | unlikely_comment, hemolyzed_comment, contaminated = (prediction == 1 | contam_comment | unlikely_comment))) %>% 
      mutate(minutes_to_verify = as.numeric((last_resulted - first_performed) / 60), date = as_date(drawn_dt_tm))
  
  times_flagged %>%
    group_by(target, prediction) %>%
    summarise(quantile(minutes_to_verify, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
  
  times_flagged %>% 
    group_by(target) %>%
    arrange(minutes_to_verify) %>%
    ggplot() + 
      xlim(0, 60) +
      stat_ecdf(data = . %>% filter(!target & prediction == 0), aes(minutes_to_verify, color = "Neither")) +
      stat_ecdf(data = . %>% filter(target & prediction == 0), aes(minutes_to_verify, color = "Technologist Only")) +
      stat_ecdf(data = . %>% filter(!target & prediction == 1), aes(minutes_to_verify, color = "Model Only")) +
      stat_ecdf(data = . %>% filter(target & prediction == 1), aes(minutes_to_verify, color = "Both")) + 
      theme_minimal()
  
  cumulative_delay <- 
    times_flagged %>% 
      arrange(date) %>%
      pivot_longer(matches("comment|target|contam"), names_to = "Comment", values_to = "Flag") %>%
      filter(Flag) %>%
      group_by(Comment, date) %>%
      summarise(sum = sum(minutes_to_verify)) %>%
      ungroup() %>% 
      summarise(median_30d = runmed(sum, k = 31), .by = c(Comment))
      
  ggplot(cumulative_delay, aes(date, median_30d, color = Comment)) + 
    geom_point()

}

compareToPOC <- function(preds = bind_cols(preds, test %>% select(UMAP1, UMAP2, likelihood_ratio) %>% mutate(umap_pred = factor(ifelse(likelihood_ratio > tar_read(umap_threshold), 1, 0), levels = c(0, 1)))), poc = read_delim("../../../bmp_umap/Data/abg_poc_ketones_input.txt")){
  
  preds <- preds %>% mutate(index = row_number())
  preds_to_join <- preds %>% transmute(index, patient_id, start = drawn_dt_tm, end = drawn_dt_tm + (4 * 60 * 60)) %>% as.data.frame()
  
  poc_glucoses <- poc %>% filter(task_assay == "Glucose rPOC" & drawn_dt_tm > "2016-12-31" & patient_id %in% preds_to_join$patient_id)
  
  poc_to_join <- poc_glucoses %>% select(patient_id, result_val, drawn_dt_tm) %>% convertAllTimesToCST() %>% as.data.frame()
  
  library(sqldf)
  joined <- sqldf("select a.*, b.*
                  from preds_to_join as a
                  left join poc_to_join as b
                    on a.patient_id = b.patient_id and
                      b.drawn_dt_tm > a.start and
                      b.drawn_dt_tm < a.end")
  
  preds_joined <- 
    joined %>% select(-patient_id) %>% as_tibble() %>% 
      mutate(poc_time_diff = as.numeric(drawn_dt_tm - start), poc_glucose = result_val) %>% select(index, poc_glucose, poc_time_diff) %>%
      group_by(index) %>% arrange(poc_time_diff) %>% slice_head(n = 1) %>% left_join(., preds) %>% ungroup() %>% drop_na(poc_glucose) %>% 
      mutate(glucose_diff = poc_glucose - glucose)
  
  poc_diff_threshold <- quantile(preds_joined$glucose_diff, probs = c(0.05), na.rm = T)[[1]]
  
  d5_poc_comparison <- 
    preds_joined %>% 
      mutate(label = case_when(contam_comment & umap_pred == 1 & glucose > 300 ~ "BOTH", 
                              umap_pred == 1 & glucose > 300 ~ "UMAP", 
                              contam_comment & glucose > 300 ~ "TECH", 
                              T ~ "NONE")) %>% mutate(label = factor(label, levels = c("BOTH", "UMAP", "TECH", "NONE"))) 
  
  d5_confirm_table <- 
    d5_poc_comparison %>%
      group_by(label) %>%
      mutate(confirmed = glucose_diff < poc_diff_threshold) %>%
      count(confirmed) %>% mutate(PPV = round(n/sum(n), digits = 2)) %>%

  
  library(geomtextpath)
  library(ggpp)
  ggplot(d5_poc_comparison, aes(glucose_diff, fill = fct_rev(label))) +
    geom_vline(xintercept = poc_diff_threshold, linetype = "dashed", color = "grey50") +
    geom_density(show.legend = F, adjust = 2, alpha = 0.5) +
    geom_textdensity(aes(label = label, hjust = label), adjust = 2, linewidth = 0, vjust = -0.2, text_smoothing = 30, fontface = "bold") + 
    scale_hjust_manual(values = c(0.1, 0.45, 0.7, 0.8)) +
    xlim(-1000, 100) +
    scale_fill_manual(values = c("NONE" = "grey50", "TECH" = scico::scico(1, begin = 0.2)[[1]], "UMAP" = scico::scico(1, begin = 0.5)[[1]], "BOTH" = scico::scico(1, begin = 0.8)[[1]])) +
    xlab("Glucose Difference (mg/dL)") + ylab("Relative Density") + ggtitle("Serum vs Capillary Glucose Differences", subtitle = "Subsequent POC glucose within 2 hours of BMP.") +
    theme(axis.text.y.left = element_blank())
  ggsave("../../Figures/Unsupervised/TestPerformance/final_UMAP/poc_confirmation_distributions.svg", width = 8, height = 6)
    
  
  d5_confirm_review <- 
    d5_poc_comparison %>%
      group_by(label) %>%
      filter(glucose > 300 & glucose_diff < poc_diff_threshold) %>%
      drop_na(any_of(matches("prior|post"))) %>%
      slice_sample(n = 25) %>% 
      reformatReviews()
  write_delim(d5_confirm_review, "../../Results/Predictions/Unsupervised/UMAP_D5_confirmation_review.tsv", delim = "\t")
  
}

reviewPowerCalculations <- function(){
  
  pwrss::pwrss.z.2props(p1 = 0.75, p2 = 0.5, margin = 0.05, alpha = 0.05, power = 0.80, alternative = "equivalent")
  
}

  
makePrePostUMAPBoxplots <- function(input = test, group_col = "prediction", analytes_to_show = c("glucose", "calcium", "chloride")){
    
    mapper = bind_cols(Analyte = analytes_to_show %>% unlist(), Group = group_col)
    
    makeBoxPlots <- function(current_col, group_col){
      
      cols_to_keep = c(paste0(current_col, "_delta_prior"), paste0(current_col, "_delta_post"))
      
      gg_input = 
        input %>% 
        select(all_of(cols_to_keep), all_of(group_col)) %>% setNames(c("Prior", "Post", "Group")) %>% 
        drop_na() %>%
        pivot_longer(!matches("Group"), names_to = "time", values_to = "value") %>% 
        mutate(time = factor(time, levels = c("Prior", "Post"), labels = c("Change from Prior", "Change to Post")),
               group = factor(Group, levels = c(0, 1), labels = c("Real", "Contaminated"))) %>% 
        slice_sample(n = 1000, by = c(Group, time))
      
      scale_y = c(quantile(gg_input$value, probs = 0.01, na.rm = T), quantile(gg_input$value, probs = 0.9995, na.rm = T))
      
      library(ggpubr)
      library(rstatix)
      means <- 
        gg_input %>% group_by(time) %>% do( tidy(t.test(data = ., value ~ fct_rev(group)))) %>% mutate(ci_label = paste0("95% CI: [", round(conf.low, digits = 1), ", ", round(conf.high, digits = 1), "]"))
      
      stats <- 
        gg_input %>% group_by(time) %>% 
        rstatix::t_test(value ~ group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj") %>%
        add_xy_position(x = "time", group = "group") %>% mutate(ci_label = means$ci_label)
      
      ggplot() + 
        geom_boxplot(data = gg_input, aes(x = time, y = value, fill = group), outlier.shape = 20, outlier.alpha = 0.1, outlier.stroke = NA, outlier.size = 4, alpha = 0.75) + 
        guides(fill = guide_legend(title = "Contamination Flag", title.theme = element_text(face = "bold"), title.position = "left", direction = "horizontal", title.hjust = 0.5)) + 
        ylab("Result Deltas") +
        ggtitle(analyte_labels[[current_col]]) +
        scale_fill_manual(values = c("grey80", "darkred")) + 
        coord_cartesian(ylim = c(scale_y[1], scale_y[2])) + 
        stat_pvalue_manual(stats, label = "{ci_label}", bracket.nudge.y = -0.1, y.position = scale_y[2]*0.95, bracket.size = 1.1, bracket.shorten = 1.2, tip.length = 0.02, color = "black", fontface = "bold") +
        theme(plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5), axis.ticks = element_blank(), legend.background = element_blank(), legend.title = element_text(size = 14, face = "bold.italic"), legend.position = c(0.5, 0.1), axis.title = element_blank(), axis.text.x = element_text(size = 14, face = "bold", color = "black"))
      
    }
    
    boxplots = furrr::future_map2(mapper$Analyte, mapper$Group, ~makeBoxPlots(.x, .y))
    for(i in c(1,4,7)) {boxplots[[i]] = boxplots[[i]] + theme(axis.title.y = element_text(size = 18, face = "bold", angle = 90, margin = margin(0, 8, 0, 0)))}
    library(ggpubr)
    ggarrange(plotlist = boxplots, ncol = 3, nrow = 3, common.legend = T, legend = "bottom")
  
}

plotAnalyteMap <- function(embed = tar_read(test_ratio_umap)){
  
  theme_set(theme_ns)
  
  umap_norm <-
    embed %>% mutate(across(any_of(lab_strings_bmp), ~scale(.), .names = "{col}_norm")) %>%
      dplyr::select(matches("UMAP|norm")) %>%
      pivot_longer(cols = matches("norm"), names_to = "Analyte", values_to = "Scaled Result") %>%
      mutate(Analyte = factor(str_replace_all(Analyte, "_norm", ""), levels = c("sodium", "chloride", "glucose", "calcium",
                                                                                "potassium_plas",  "co2_totl", "creatinine",
                                                                                "bun", "anion_gap")))
  gg_results_by_analyte <-
    ggplot(umap_norm %>% group_by(Analyte), aes(UMAP1, UMAP2, color = `Scaled Result`)) +
      geom_point(shape = ".") +
      scico::scale_color_scico(palette = "berlin", limits = c(-3, 3), breaks = c(-2.75, 2.65), labels = c("Lowest", "Highest"), oob = scales::squish) +
      guides(color = guide_colorbar(title = "Result Value", ticks = F, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 10, face = "bold"), title.position = "top", title.hjust = 0.5)) +
      facet_wrap(~Analyte, ncol = 3, labeller = labeller(Analyte = analyte_labels)) +
      theme(panel.spacing = unit(0, "lines"), strip.text = element_text(size = 12, face = "bold.italic"), legend.direction = "horizontal", legend.position = "bottom", legend.key.size = unit(0.25, "in"), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())
    ggsave(here::here("../../Figures/Unsupervised/Supplemental/UMAP_by_analyte_berlin.png"), gg_results_by_analyte, width = 8, height = 8, dpi = 600)
    ggsave(here::here("../../Figures/Unsupervised/Supplemental/UMAP_by_analyte.pdf"), gg_results_by_analyte, width = 8, height = 8)
  
  fluid_coordinates <- 
    list(
      "D5" = list(xmin = 0.2, xmax = 2, ymin = -0.6, ymax = 0.6),
      "NS" = list(xmin = -2, xmax = -0.5, ymin = -0.3, ymax = 0.2),
      "LR" = list(xmin = 0, xmax = 0.6, ymin = -1.1, ymax = -0.3)
    )
  
  ns_area <- 
    umap_norm %>% filter(Analyte %in% c("chloride", "calcium") & between(UMAP1, fluid_coordinates[["NS"]][["xmin"]], fluid_coordinates[["NS"]][["xmax"]]) & between(UMAP2, fluid_coordinates[["NS"]][["ymin"]], fluid_coordinates[["NS"]][["ymax"]]))
  
  ns_zoom <- 
    ns_area %>% 
    ggplot(aes(UMAP1, UMAP2, color = `Scaled Result`)) +       
    geom_point(shape = 20, stroke = 1, size = 4, alpha = 1) +
    scico::scale_color_scico(palette = "berlin", limits = c(-3, 3), breaks = c(-2.75, 2.65), labels = c("Lowest", "Highest"), oob = scales::squish) +
    guides(color = guide_colorbar(title = "Result Value", ticks.colour = NA, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 10, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
    facet_wrap(~Analyte, ncol = 1, labeller = labeller(Analyte = analyte_labels)) +
    theme(strip.text = element_text(size = 14, hjust = 0.1, face = "bold.italic"), legend.position = "none", legend.direction = "horizontal", legend.key.size = unit(0.5, "in"), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  ggsave(here::here("../../Figures/Unsupervised/Supplemental/ns_zoom_only.pdf"), ns_zoom, width = 6, height = 6, dpi = 600)
  
  lr_area <- 
    umap_norm %>% filter(Analyte %in% c("calcium") & between(UMAP1, fluid_coordinates[["LR"]][["xmin"]], fluid_coordinates[["LR"]][["xmax"]]) & between(UMAP2, fluid_coordinates[["LR"]][["ymin"]], fluid_coordinates[["LR"]][["ymax"]]))
  
  lr_zoom <- 
    lr_area %>% 
    ggplot(aes(UMAP1, UMAP2, color = `Scaled Result`)) +       
    geom_point(shape = 20, stroke = 1, size = 4, alpha = 1) +
    scico::scale_color_scico(palette = "berlin", limits = c(-3, 3), breaks = c(-2.75, 2.65), labels = c("Lowest", "Highest"), oob = scales::squish) +
    guides(color = guide_colorbar(title = "Result Value", ticks.colour = NA, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 10, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
    facet_wrap(~Analyte, ncol = 1, labeller = labeller(Analyte = analyte_labels)) +
    theme(strip.text = element_text(size = 14, hjust = 0.1, face = "bold.italic"), legend.position = "none", legend.direction = "horizontal", legend.key.size = unit(0.5, "in"), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  ggsave(here::here("../../Figures/Unsupervised/Supplemental/lr_zoom_only.pdf"), lr_zoom, width = 6, height = 3, dpi = 600)
  
  d5_area <- 
    umap_norm %>% filter(Analyte %in% c("glucose") & between(UMAP1, fluid_coordinates[["D5"]][["xmin"]], fluid_coordinates[["D5"]][["xmax"]]) & between(UMAP2, fluid_coordinates[["D5"]][["ymin"]], fluid_coordinates[["D5"]][["ymax"]]))
  
  d5_zoom <- 
    d5_area %>% 
    ggplot(aes(UMAP1, UMAP2, color = `Scaled Result`)) +       
    geom_point(shape = 20, stroke = 1, size = 4, alpha = 1) +
    scico::scale_color_scico(palette = "berlin", limits = c(-3, 3), breaks = c(-2.75, 2.65), labels = c("Lowest", "Highest"), oob = scales::squish) +
    guides(color = guide_colorbar(title = "Result Value", ticks.colour = NA, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 10, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
    facet_wrap(~Analyte, ncol = 1, labeller = labeller(Analyte = analyte_labels)) +
    theme(strip.text = element_text(size = 14, hjust = 0.1, face = "bold.italic"), legend.position = "none", legend.direction = "horizontal", legend.key.size = unit(0.5, "in"), axis.text = element_blank(), axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank(), axis.line = element_blank())
  ggsave(here::here("../../Figures/Unsupervised/Supplemental/d5_zoom_only.pdf"), d5_zoom, width = 6, height = 3, dpi = 600)
  

  
  ##### Plot Pre Post #####
  glucose_box_input <- embed %>% mutate(d5_area = ifelse(between(UMAP1, fluid_coordinates[["D5"]][["xmin"]], fluid_coordinates[["D5"]][["xmax"]]) & between(UMAP2, fluid_coordinates[["D5"]][["ymin"]], fluid_coordinates[["D5"]][["ymax"]]), "D5 Contamination", "Non-D5"))
  
  glucose_boxplots <- makePrePostUMAPBoxplots(glucose_box_input, d5_area, c("glucose"))
  
  
}

compareToManualReview <- function(predictions = tar_read(test_ratio_umap), reviewed = readxl::read_xlsx("../../Data/reviewed_concordant.xlsx"), positives = read_delim("../../Data/final_consecutive_positives_reviewed.csv")){
  
  reviewed_current <- reviewed %>% filter(time_point == "Current")
  
  reviewed_current_replaced <- 
    reviewed_current %>% 
      mutate(
          drawn_dt_tm = as_datetime(drawn_dt_tm),
          expert_review = 
               factor(case_when(
                 grepl("Cont|NS", final_prediction, ignore.case = T) ~ "Contaminated",
                 grepl("Real", final_prediction, ignore.case = T) ~ "Real",
                 grepl("K2EDTA", contaminating_fluid, ignore.case = T) ~ "K2EDTA",
                 grepl("TPN", contaminating_fluid, ignore.case = T) ~ "TPN",
                 grepl("Hemolyzed", contaminating_fluid, ignore.case = T) ~ "Hemolyzed")))
  
  joined <- 
    inner_join(predictions, reviewed_current_replaced) %>% 
      mutate(final_prediction = factor(final_prediction), 
             flagged = factor((contam_comment | unlikely_comment), levels = c(F, T), labels = c("Real", "Contaminated")),
             prediction_factor = factor(prediction, levels = c(0, 1), labels = c("Real", "Contaminated")))
  
  review_pr <- 
    pr_curve(joined, final_prediction, likelihood_ratio) %>% mutate(standard = "Expert Review")
  
  points <- review_pr %>% filter(near(.threshold, 39, tol = 0.1) | near(.threshold, 10, tol = 0.1) | near(.threshold, 100, tol = 1)) %>% mutate(color = scico::scico(3, direction = -1, palette = "bilbao", begin = 0.2, end = 0.9))
  
  ggplot(review_pr, aes(recall, precision)) + 
    geom_path(linewidth = 1.5, color = "black", alpha = 0.5) +
    geom_point(data = points, aes(recall, precision, color = color), size = 6) + 
    coord_equal(ratio = 1) +
    scale_color_identity() + 
    annotate("text", x = 0.95, y = 0.1, hjust = 1, label = paste0("AUC = ", sprintf("%.2f", round(pr_auc(joined, final_prediction, likelihood_ratio) %>% pluck(".estimate"), digits = 2))), size = 4, fontface = "bold") +
    xlim(0,1) + ylim(0,1) + 
    xlab("Sensitivity") + ylab("Positive Predictive Value")
  ggsave("../../Figures/Unsupervised/ManualReview/pr_curve_with_dots.svg", width = 8, height = 4)
  
  alluvial <- 
    joined %>%
      mutate(fill = factor(paste0(prediction_factor, flagged), labels = c("Both Flagged", "UMAP Flagged", "Current Flagged", "Both Real"))) %>%
      slice_sample(n = 50, by = fill) %>%
      count(prediction_factor, final_prediction, flagged, fill)

  library(ggalluvial)
  library(scico)
  fills = scico(13, palette = "davos", direction = -1)
  
  ggplot(alluvial, aes(axis2 = fill, axis1 = final_prediction, y = n)) + 
    geom_alluvium(aes(fill = fill, color = fill), alpha = 0.75, show.legend = F, width = 0.25) +
    geom_stratum(fill = "white", alpha = 0, show.legend = F, width = 0.25) +
    scale_x_discrete(expand = c(0, 1)) +
    coord_flip() +
    scale_fill_manual(values = c("Real" = fills[1], "Contaminated" = fills[12], "Both Flagged" = fills[13], "Current Flagged" = fills[7], "UMAP Flagged" = fills[9], "Both Real" = fills[2])) +
  scale_color_manual(values = c("Real" = fills[1], "Contaminated" = fills[12], "Both Flagged" = fills[13], "Current Flagged" = fills[7], "UMAP Flagged" = fills[9], "Both Real" = fills[2])) +
    theme(axis.line = element_blank(), axis.title.x.bottom = element_blank(), axis.text = element_blank())
  ggsave("../../Figures/Unsupervised/ManualReview/random_sample_sankey.svg", width = 8, height = 4)
  

  
  
  reviewed_positives <- positives %>% filter(time_point == "Current") %>% mutate(prediction = factor(prediction, levels = c(0, 1)))

  joined_pos <- 
    reviewed_positives %>%
    mutate(final_prediction = factor(final_prediction),
           flagged = factor(flagged, levels = c(F, T), labels = c("Real", "Contaminated")),
           prediction_factor = factor(prediction, levels = c(0, 1), labels = c("Real", "Contaminated")))
  
  pos_long <- joined_pos %>% filter(flagged == "Real" & final_prediction == "Contaminated") %>% pivot_longer(lab_strings_bmp, names_to = "analyte", values_to = "result") %>% filter()
  
  critical_ranges <- critical_ranges %>% setNames(janitor::make_clean_names(names(critical_ranges)))
  critical = map2(pos_long$analyte, pos_long$result, ~(.y < critical_ranges[[.x]][["min"]] | .y > critical_ranges[[.x]][["max"]])) %>% unlist() %>% table()
  
  alluvial_pos <- 
    joined_pos %>%
      mutate(fill = factor(paste0(final_prediction, flagged), labels = c("UMAP False Positives", "Current False Negatives", "Concordant Positives"))) %>%
      count(prediction_factor, final_prediction, flagged, fill)
  
  ggplot(alluvial_pos, aes(axis3 = "UMAP Flagged", axis2 = fill, y = n)) + 
    geom_alluvium(aes(fill = fill), show.legend = F, width = 0.25, alpha = 0.5) +
    geom_stratum(fill = "white", alpha = 0, show.legend = F, width = 0.25) +
    scale_x_discrete(expand = c(0, 1)) +
    coord_flip() +
    scale_fill_manual(values = c(fills[13], fills[8], fills[5])) +
    theme(axis.line = element_blank(), axis.title.x.bottom = element_blank(), axis.text = element_blank())
  ggsave("../../Figures/Unsupervised/ManualReview/consecutive_positives_sankey.svg", width = 8, height = 4)
  
  
predictions   
  
daily <- tar_read(preprocessed_bmp_inputs)[[1]][["data"]] %>% 
  mutate(day = floor_date(drawn_dt_tm, "year")) %>% 
  group_by(day) %>% 
  summarise(daily_flags_current = sum(unlikely_comment | contam_comment), n = n(), prop = daily_flags_current/n)

daily <- daily %>% mutate(group = factor(ifelse(day < "2019-10-01", 1, 2))) %>% summarise(n = mean(n), flag = mean(daily_flags_current), .by = group)

 monthly_flags <-  
  predictions %>% 
    mutate(month = floor_date(drawn_dt_tm, "month"), 
           year = floor_date(drawn_dt_tm, "year")) %>% 
    group_by(month) %>% 
    summarise(daily_flags_low = sum(likelihood_ratio > 100) / n(),
              daily_flags_mid = sum(likelihood_ratio > 39) / n(),
              daily_flags_high = sum(likelihood_ratio > 10) / n(),
              daily_flags_current = sum(unlikely_comment | contam_comment) / n()) %>% 
   pivot_longer(matches("flags"), names_to = "Threshold", values_to = "Flag Rate") %>% 
   mutate(Threshold = factor(Threshold, levels = c("daily_flags_high", "daily_flags_mid", "daily_flags_low")))
 
 
 
 ggplot(monthly_flags, aes(month, `Flag Rate`, color = Threshold, fill = Threshold)) + 
    geom_smooth(data = monthly_flags %>% filter(!grepl("current", Threshold)), span = 0.3, alpha = 0.25, linetype = "dashed", show.legend = F) + 
    geom_smooth(data = monthly_flags %>% filter(is.na(Threshold)), span = 0.3, se = F, color = "black", linewidth = 4, lineend = "round", show.legend = F) +  
    scale_y_continuous(limits = c(0, 0.006), breaks = seq(0, 0.006, by = 0.002), labels = c("0", "2", "4", "6")) + 
    scale_x_datetime(date_breaks = "year", date_minor_breaks = "month", labels = date_format("%Y")) + 
    scico::scale_color_scico_d(palette = "bilbao", begin = 0.2, end = 0.9) + 
    scico::scale_fill_scico_d(palette = "bilbao", begin = 0.2, end = 0.9) +
    xlab("Monthly Rolling Average") + ylab("Expected Alarms Per 1000 BMPs") +
    theme()
 ggsave("../../Figures/Unsupervised/TestPerformance/final_clinchem_UMAP/alarm_rate_over_time.svg", width = 6, height = 4)
 
 resampled <- predictions %>% slice_sample(prop = 0.1, weight_by = likelihood_ratio) %>% mutate(color = case_when(likelihood_ratio <= 10 ~ scico::scico(30, palette = "bilbao", begin = 0.05, end = 0.9)[1], between(likelihood_ratio, 10, 39) ~ scico::scico(3, palette = "bilbao", begin = 0.2, end = 0.9)[1], between(likelihood_ratio, 39, 100) ~ scico::scico(3, palette = "bilbao", begin = 0.2, end = 0.9)[2], likelihood_ratio >= 100 ~ scico::scico(3, palette = "bilbao", begin = 0.2, end = 0.9)[3]))
 ggplot() +
   geom_point(data = resampled, aes(UMAP1, UMAP2, color = color), alpha = 0.5, size = 2, shape = 20, stroke = 0.1) + 
   scale_color_identity() + 
   theme(legend.position = c(0.2, 0.2), axis.text = element_blank())
 ggsave("../../Figures/Unsupervised/TestPerformance/final_clinchem_UMAP/binned_predictions_for_pr_inset.svg", width = 4, height = 4)
  
}
