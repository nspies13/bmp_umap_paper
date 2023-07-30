makeDensityFigure <- function(embed){

ggplot() +
  ggpointdensity::geom_pointdensity(data = embed, aes(x = UMAP1, y = UMAP2), adjust = 5, size = 3, shape = 20, stroke = NA) +
  scale_color_viridis_c() +
  guides(color = guide_colorbar(title = "Number of Points", ticks.colour = NA, title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold"), title.position = "top", title.hjust = 0.5)) +
  theme(legend.direction = "horizontal", legend.position = c(0.2, 0.2), legend.key.size = unit(0.5, "in"), axis.text = element_blank())
ggsave(here::here("../Figures/umap_results_patient_density.png"), width = 11, height = 11)
ggsave(here::here("../Figures/umap_results_patient_density.svg"), width = 11, height = 11)

}

makeSimEmbeddingFigures <- function(embed){

  gg_label <-
    ggplot() +
      geom_point(data = embed, aes(x = UMAP1, y = UMAP2, color = label), alpha = 0.15, size = 1, shape = 20, stroke = NA) +
      scale_color_manual(values = color_map_global) +
      annotate("text", label = "NS", x = 5.25, y = -2.5, size = 12, fontface = "bold.italic", color = color_map_global["NS"]) +
      annotate("text", label = "LR", x = 3, y = -3, size = 12, fontface = "bold.italic", color = color_map_global["LR"]) +
      annotate("text", label = "D5NS", x = -8.15, y = 1.65, size = 12, fontface = "bold.italic", color = color_map_global["D5NS"]) +
      annotate("text", label = "D5LR", x = -9, y = 2, size = 12, fontface = "bold.italic", color = color_map_global["D5LR"]) +
      theme(legend.position = "none", axis.text = element_blank(), legend.background = element_blank())
  ggsave(here::here("../Figures/sim_embedding_by_label.png"), gg_label, width = 11, height = 11)
  ggsave(here::here("../Figures/sim_embedding_by_label.svg"), gg_label, width = 11, height = 11)


  gg_percent <-
    ggplot() +
      geom_point(data = embed, aes(x = UMAP1, y = UMAP2, color = percent_fluid), alpha = 0.15, size = 1, shape = 20, stroke = NA) +
      scale_color_viridis_c(option = "B") +
      annotate("text", label = "NS", x = 5.25, y = -2.5, size = 12, fontface = "bold.italic", color = color_map_global["NS"]) +
      annotate("text", label = "LR", x = 3, y = -3, size = 12, fontface = "bold.italic", color = color_map_global["LR"]) +
      annotate("text", label = "D5NS", x = -8.15, y = 1.65, size = 12, fontface = "bold.italic", color = color_map_global["D5NS"]) +
      annotate("text", label = "D5LR", x = -9, y = 2, size = 12, fontface = "bold.italic", color = color_map_global["D5LR"]) +
      guides(colour = guide_colourbar(title = "Mixture Ratio", title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold.italic"), title.position="top", title.hjust = 0.5)) +
      theme(legend.position = c(0.2, 0.1), legend.key.size = unit(0.5, "in"), legend.text = element_text(size = 12, face = "bold.italic"), legend.direction = "horizontal", axis.text = element_blank(), legend.background = element_blank())
  ggsave(here::here("../Figures/sim_embedding_by_percent.png"), gg_percent, width = 11, height = 11)
  ggsave(here::here("../Figures/sim_embedding_by_percent.svg"), gg_percent, width = 11, height = 11)

  gg_sim <- ggpubr::ggarrange(gg_label, gg_percent, nrow = 1, ncol = 2, labels = "AUTO", font.label = list(size = 24))
  ggsave(here::here("../Figures/sim_embedding_combined.png"), gg_sim, width = 22, height = 11)

}

makeLabwiseEmbeddings <- function(embed){

  umap_norm <-
    embed %>% mutate(across(any_of(lab_strings_bmp), ~scale(.), .names = "{col}_norm")) %>%
      dplyr::select(matches("UMAP|norm")) %>%
      pivot_longer(cols = matches("norm"), names_to = "Analyte", values_to = "Scaled Result") %>%
      mutate(Analyte = factor(str_replace_all(Analyte, "_norm", ""), levels = c("sodium", "chloride", "glucose", "calcium",
                                                                                "potassium_plas",  "co2_totl", "creatinine",
                                                                                "bun", "anion_gap")))
  gg_results_by_analyte <-
    ggplot(umap_norm %>% group_by(Analyte), aes(UMAP1, UMAP2, color = `Scaled Result`)) +
      geom_point(shape = 20, stroke = 1, size = 1, alpha = 0.75) +
      scale_color_viridis_c(limits = c(-3, 3), breaks = c(-2.5, 2.45), labels = c("Lowest", "Highest"), oob = scales::squish, option = "F", end = 0.95) +
      guides(color = guide_colorbar(title = "Result Value", title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold"), title.position = "top", title.hjust = 0.5)) +
      facet_wrap(~Analyte, ncol = 3, labeller = labeller(Analyte = analyte_labels)) +
      theme(strip.text = element_text(size = 18, face = "bold.italic"), legend.direction = "horizontal", legend.position = "bottom", legend.key.size = unit(0.5, "in"), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  ggsave(here::here("../Figures/UMAP_by_Analyte.png"), gg_results_by_analyte, width = 11, height = 11)
  ggsave(here::here("../Figures/UMAP_by_Analyte.svg"), gg_results_by_analyte, width = 11, height = 11)

  gg_zoom <- 
    umap_norm %>% filter(Analyte %in% c("chloride", "calcium", "potassium_plas", "co2_totl") & between(UMAP1, 3.5, 6.5) & between(UMAP2, -4, -1)) %>% 
      ggplot(aes(UMAP1, UMAP2, color = `Scaled Result`)) +       
      geom_point(shape = 20, stroke = NA, size = 4, alpha = 0.75) +
      scale_color_viridis_c(limits = c(-3, 3), breaks = c(-2.45, 2.38), labels = c("Lowest", "Highest"), oob = scales::squish, option = "E") +
      guides(color = guide_colorbar(title = "Result Value", ticks.colour = NA, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 10, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
      facet_wrap(~Analyte, ncol = 2, labeller = labeller(Analyte = analyte_labels)) +
      theme(strip.text = element_text(size = 14, hjust = 0.1, face = "bold.italic"), legend.position = "none", legend.direction = "horizontal", legend.key.size = unit(0.5, "in"), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  ggsave(here::here("../Figures/embed_zoom_only.pdf"), gg_zoom, width = 5.4, height = 6.2, dpi = 600)
    
  ggplot(embed, aes(UMAP1, UMAP2, color = ntile(kde_results_umap, 100))) +
    geom_point(size = 3, shape = 20, stroke = 0) +
    coord_cartesian(xlim = c(-2, 15)) + 
    scale_color_viridis_c(name = "Density", option = "D", direction = 1,
                          limit = c(0, 100), breaks = c(0, 100),
                          labels = c("Lowest", "Highest"), trans = scales::pseudo_log_trans(sigma = 1)) +
    guides(color = guide_colorbar(title = "Density", ticks.colour = NA, title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold"), title.position = "top", title.hjust = 0.5)) +
    theme(legend.direction = "horizontal", legend.position = "none", legend.key.size = unit(0.5, "in"), axis.text = element_blank(), legend.background = element_blank(), legend.box.background = element_blank()) +
    patchwork::inset_element(gg_zoom, 0.5, 0, 1, 1)
  ggsave(here::here("../Figures/embed_contam_zoom.pdf"), width = 11, height = 8, dpi = 600)

}

plotKDE <- function(embed_kde){
  
  gg_kde <-
    ggplot(embed_kde, aes(UMAP1, UMAP2, color = ntile(kde_results_umap, 10000))) +
    geom_point(size = 3, shape = 20, stroke = NA) +
    scale_color_viridis_c(name = "Density", option = "D", direction = 1,
                          limit = c(0, 10000), breaks = c(175, 7250),
                          labels = c("Lowest", "Highest"), trans = scales::sqrt_trans()) +
    guides(color = guide_colorbar(title = "Density", ticks.colour = NA, title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
    theme(legend.direction = "horizontal", legend.position = c(0.2, 0.2), legend.key.size = unit(0.5, "in"), axis.text = element_blank(), legend.background = element_blank(), legend.box.background = element_blank())
  
  gg_kde
  
  ggsave("../Figures/figure1.pdf", gg_kde, width = 8, height = 8)
  
}

plotEnrichmentScores <- function(enrichment_scores, unsup_train_kde, sim_results_kde){
  
  gg_enrich <- 
    ggplot(enrichment_scores, aes(bin_x, bin_y, fill = ntile(enrichment, n = 100))) + 
      geom_tile() + scale_fill_viridis_c(labels = c("Lowest", "Highest"), limits = c(90, 100), breaks = c(92, 98), oob = scales::squish) +
      guides(fill = guide_colorbar(title = "Enrichment Score", ticks.colour = NA, title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
      xlab("UMAP1 Binned") + ylab("UMAP2 Binned") + 
      theme(legend.direction = "horizontal", legend.position = c(0.2, 0.2), legend.key.size = unit(0.5, "in"), axis.text = element_blank(), legend.background = element_blank(), legend.box.background = element_blank())
  ggsave(here::here("../Figures/enrichment_scores.png"), gg_enrich, width = 11, height = 11, dpi = 600)
  
  gg_train <- 
    ggplot(unsup_train_kde, aes(UMAP1, UMAP2, color = ntile(kde_results_umap, 10000))) +
      geom_point(size = 2, shape = 20, stroke = NA) +
      scale_color_viridis_c(name = "Density", option = "D", direction = 1,
                          limit = c(0, 10000), breaks = c(175, 7250),
                          labels = c("Lowest", "Highest"), trans = scales::sqrt_trans()) +
      guides(color = guide_colorbar(title = "Patient\nDensity", ticks.colour = NA, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 8, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
      theme(legend.position = c(0.2, 0.2), plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  
  gg_sim <- 
    ggplot(sim_results_kde, aes(UMAP1, UMAP2, color = ntile(kde_results_umap, 10000))) +
      geom_point(size = 2, shape = 20, stroke = NA) +
      scale_color_viridis_c(name = "Density", option = "D", direction = 1,
                            limit = c(0, 10000), breaks = c(175, 7250),
                            labels = c("Lowest", "Highest"), trans = scales::sqrt_trans()) +
      guides(color = guide_colorbar(title = "Contamination\nDensity", ticks.colour = NA, title.theme = element_text(size = 12, face = "bold.italic"), label.theme = element_text(size = 8, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
      theme(legend.position = c(0.2, 0.2), plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  base <- ggpubr::ggarrange(gg_sim, gg_train, nrow = 2)
  
  gg_enrich <- 
    ggplot(enrichment_scores, aes(bin_x, bin_y, fill = ntile(enrichment, n = 100))) +
      geom_tile() + scale_fill_viridis_c(labels = c("Lowest", "Highest"), limits = c(90, 100), breaks = c(91.5, 98.5), oob = scales::squish) +
      guides(fill = guide_colorbar(title = "Enrichment Score", ticks.colour = NA, title.theme = element_text(size = 16, face = "bold.italic"), label.theme = element_text(size = 12, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
      theme(legend.direction = "horizontal", legend.position = c(0.3, 0.2), legend.key.size = unit(0.5, "in"), legend.background = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank())
    
  gg_combo <- ggpubr::ggarrange(base, gg_enrich, nrow = 1, widths = c(0.4, 0.6))
  ggsave(here::here("../Figures/enrichment_combined.png"), gg_combo, bg = "white", width = 11, height = 8, dpi = 600)

}

makeDKAembeddingFigure <- function(embed, input){
  
  gg_input <- embed %>% pluck("template") %>% bind_cols(label = input$label)
  
  ggplot() + 
    geom_point(data = gg_input %>% filter(label == "D5"), aes(UMAP1, UMAP2), color = viridis::viridis(9)[[2]], shape = 20, size = 4, stroke = NA, alpha = 0.5) + 
    geom_point(data = gg_input %>% filter(label == "DKA"), aes(UMAP1, UMAP2), color = viridis::viridis(9)[[6]],  shape = 20, size = 4, stroke = NA, alpha = 0.5)

  }

makeClinicalFigures <- function(unsup_results_test, unsup_results_kde, clinicals){
  
  gg_input <- 
    left_join(unsup_results_test, unsup_results_kde, by = "index") %>% 
    left_join(clinicals %>% group_by(patient_id, drawn_dt_tm, performed_dt_tm, verified_dt_tm) %>% slice_tail(n = 1)) %>% 
    distinct()
  
  gg_icu <- 
    ggplot() + 
      ggpointdensity::geom_pointdensity(data = gg_input %>% filter(icu), aes(UMAP1, UMAP2), size = 2, adjust = 0.1, stroke = NA, shape = 20) + 
      scale_color_viridis_c() + 
      guides(color = guide_colorbar(title = "Critical Care\nDensity", ticks.colour = NA, title.theme = element_text(size = 8, face = "bold.italic"), label.theme = element_text(size = 6, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
      theme(legend.position = "none", plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.background = element_blank())
  
  
  gg_emergency <- 
    ggplot() + 
      ggpointdensity::geom_pointdensity(data = gg_input %>% filter(emergency), aes(UMAP1, UMAP2), size = 2, adjust = 0.1, stroke = NA, shape = 20) +
      scale_color_viridis_c() + 
      guides(color = guide_colorbar(title = "Emergency\nDensity", ticks.colour = NA, title.theme = element_text(size = 8, face = "bold.italic"), label.theme = element_text(size = 6, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
      theme(legend.position = "none", plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.background = element_blank())
    
  
  gg_not_icu <-   
    ggplot() + 
      ggpointdensity::geom_pointdensity(data = gg_input %>% filter(!icu & !emergency), aes(UMAP1, UMAP2), size = 2, adjust = 0.1, stroke = NA, shape = 20) + 
      scale_color_viridis_c(labels = c("Lowest", "Highest"), breaks = c(300, 1600), oob = scales::squish) + 
      guides(color = guide_colorbar(title = "Patient\nDensity", ticks.colour = NA, title.theme = element_text(size = 8, face = "bold.italic"), label.theme = element_text(size = 6, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
      theme(legend.position = c(0.1, 0.1), plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.2, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.background = element_blank())
    
  gg_loc <- ggpubr::ggarrange(gg_not_icu, gg_icu, gg_emergency, nrow = 1)
  ggsave(here::here("../Figures/patient_location_densities.png"), gg_loc, width = 8.5, height = 4, dpi = 600, bg = "white")
  
}

makeCommentFigures <- function(unsup_results_test, unsup_results_kde, comments){
  
  gg_input <- 
    left_join(unsup_results_test, unsup_results_kde, by = "index") %>% 
    left_join(comments %>% group_by(patient_id, drawn_dt_tm, performed_dt_tm, verified_dt_tm) %>% slice_tail(n = 1)) %>% 
    distinct()
  
  gg_input %>% filter(code_comment) %>% ggplot(aes(UMAP1, UMAP2)) + ggpointdensity::geom_pointdensity(size = 2, adjust = 0.1, stroke = NA, shape = 20) + 
    scale_color_viridis_c()

  gg_input %>% filter(unlikely_comment) %>% ggplot(aes(UMAP1, UMAP2)) + ggpointdensity::geom_pointdensity(size = 2, adjust = 0.1, stroke = NA, shape = 20) + 
    scale_color_viridis_c()
  
  gg_input %>% filter(contamination_comment) %>% ggplot(aes(UMAP1, UMAP2)) + ggpointdensity::geom_pointdensity(size = 2, adjust = 0.1, stroke = NA, shape = 20) + 
    scale_color_viridis_c(labels = c("Lowest", "Highest"), breaks = c(12, 78), oob = scales::squish) +      
    guides(color = guide_colorbar(title = "Density", ticks.colour = NA, title.theme = element_text(size = 8, face = "bold.italic"), label.theme = element_text(size = 6, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
    theme(legend.position = c(0.2, 0.2), plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.background = element_blank())
  ggsave("../Figures/contamination_comment.png", width = 5, height = 5)
  
  gg_input %>% filter(glycolysis_comment) %>% ggplot(aes(UMAP1, UMAP2)) + ggpointdensity::geom_pointdensity(size = 2, adjust = 0.1, stroke = NA, shape = 20) + 
    scale_color_viridis_c()
  
    
}

makeLUOfigures <- function(unsup_results_test, unsup_results_kde, luo){
  
  gg_input <- 
    left_join(unsup_results_test, unsup_results_kde, by = "index") %>% 
    left_join(luo %>% group_by(patient_id, drawn_dt_tm, performed_dt_tm, verified_dt_tm) %>% slice_tail(n = 1)) %>% 
    distinct()
  
  gg_long <- 
    gg_input %>% 
      dplyr::select(UMAP1, UMAP2, matches("LUO-H|LUO-I|LUO-L")) %>% 
      pivot_longer(matches("LUO"), names_to = "lab", values_to = "value") %>%
      mutate(value = as.numeric(value))
  
  gg_input %>% mutate(value = as.numeric(`LUO-H`)) %>% select(UMAP1, UMAP2, value) %>%
    mutate(value = ifelse(value > 50, 50, value)) %>% 
    ggplot(aes(UMAP1, UMAP2, color = value)) + 
    geom_point(size = 2, shape = 20, stroke = NA, alpha = 0.5) + 
    scale_color_viridis_c(oob = scales::squish) + 
    guides(color = guide_colorbar(title = "Density", ticks.colour = NA, title.theme = element_text(size = 8, face = "bold.italic"), label.theme = element_text(size = 6, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
    theme(legend.position = "none", plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.background = element_blank())
  ggsave("../Figures/hemolysis_figure.png", width = 11, height = 11, dpi = 600)
  
  gg_long %>% ggplot(aes(UMAP1, UMAP2, color = value)) + 
    geom_point(size = 2, stroke = NA, shape = 20) +
    facet_wrap(~lab, scales = "free") + 
    scale_color_viridis_c(oob = scales::squish) + 
    guides(color = guide_colorbar(title = "Density", ticks.colour = NA, title.theme = element_text(size = 8, face = "bold.italic"), label.theme = element_text(size = 6, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) + 
    theme(legend.position = "none", plot.background = element_blank(), legend.direction = "horizontal", legend.key.height = unit(0.25, "in"), legend.key.width = unit(0.25, "in"), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.background = element_blank())
  ggsave("../Figures/luo_figures.png", width = 8.5, height = 5, dpi = 600)
  
}

makeNSprojectionFigure <- function(unsup_results_embed, unsup_results_test){
  
  indexes <- unsup_results_embed %>% filter(index != 1057279 & UMAP1 > 4.085 & UMAP2 > -2.5 & UMAP2 < -1.5) %>% select(index) %>% pluck(1)
  
  prior <- 
    unsup_results_test %>% 
    group_by(patient_id) %>% arrange(patient_id, drawn_dt_tm) %>% distinct(drawn_dt_tm, .keep_all = T) %>%
    mutate(across(all_of(lab_strings_bmp), ~dplyr::lag(.x) - .x, .names = "{.col}-Prior")) %>% filter(`potassium_plas-Prior` > 0)
  
  post <- 
    unsup_results_test %>% 
    group_by(patient_id) %>% arrange(patient_id, drawn_dt_tm) %>% distinct(drawn_dt_tm, .keep_all = T) %>%
    mutate(across(all_of(lab_strings_bmp), ~dplyr::lead(.x) - .x, .names = "{.col}-Post")) %>% filter(`potassium_plas-Post` > 0)
  
  input <- 
    left_join(unsup_results_test %>% filter(index %in% indexes), prior) %>% 
    left_join(post) %>% select(matches("potassium|calcium|chloride|co2")) %>% 
    mutate(across(!matches("Prior|Post"), ~.x * 0), patient = row_number())
  
  input_long <- 
    input %>% 
      na.omit()  %>%
      mutate(patient = fct_reorder(as_factor(patient), `chloride-Prior`, .na_rm = F)) %>%
      pivot_longer(cols = !matches("patient"), names_to = c("lab", "time"), values_to = "result", names_sep = "-") %>% 
      mutate(time = as_factor(ifelse(is.na(time), " ", time))) %>%
      mutate(time = factor(time, levels = c("Prior", " ", "Post"))) %>%
      mutate(lab = factor(lab, levels = c("chloride", "calcium", "potassium_plas", "co2_totl")))
  
  input_long %>% 
    ggplot(aes(time, result)) + 
    geom_line(aes(group = patient, color = patient), linewidth = 1, lineend = "round") +
    scale_color_viridis_d(direction = -1, begin = 0.75, end = 0.95, option = "E") + 
    geom_vline(xintercept = " ", linetype = "dashed") + 
    facet_wrap(~lab, scales = "free", labeller = labeller(lab = analyte_labels), nrow = 1) + 
    theme(strip.text = element_text(size = 14, face = "bold.italic"), legend.position = "none", axis.title = element_blank())
  ggsave("../Figures/deltas_fig2b.pdf", width = 12.5, height = 4.5, dpi = 600)
  
}

makeDextroseFigure <- function(embed, dka){
  
  embed <- unsup_results_embed
  
  gg_embed <- 
    embed %>% 
      ggplot(aes(UMAP1, UMAP2, color = glucose)) +       
      geom_point(shape = 20, stroke = 1, size = 4, alpha = 1) +
      scale_color_viridis_c(limits = c(100, 500), breaks = c(130, 460), labels = c("<100", ">500"), oob = scales::squish, option = "B", end = 0.8) +
      coord_cartesian(xlim = c(-12, -2)) +
      guides(color = guide_colorbar(title = "Glucose (mg/dL)", title.theme = element_text(size = 10, face = "bold.italic"), ticks.colour = NA, label.theme = element_text(size = 8, face = "bold.italic"), title.position = "top", title.hjust = 0.5)) +
      theme(strip.text = element_text(size = 18, face = "bold.italic"), legend.direction = "horizontal", legend.position = c(0.55, 0.95), legend.background = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())
  
  
  indexes <- unsup_results_embed %>% arrange(UMAP1) %>% slice_head(n = 100) %>% select(index) %>% pluck(1)
  
  prior <- 
    unsup_results_test %>% 
    group_by(patient_id) %>% arrange(patient_id, drawn_dt_tm) %>% distinct(drawn_dt_tm, .keep_all = T) %>%
    mutate(across(all_of(lab_strings_bmp), ~dplyr::lag(.x) - .x, .names = "{.col}-Prior"))
  
  post <- 
    unsup_results_test %>% 
    group_by(patient_id) %>% arrange(patient_id, drawn_dt_tm) %>% distinct(drawn_dt_tm, .keep_all = T) %>%
    mutate(across(all_of(lab_strings_bmp), ~dplyr::lead(.x) - .x, .names = "{.col}-Post"))
  
  input <- 
    left_join(unsup_results_test %>% filter(index %in% indexes), prior) %>% 
    left_join(post) %>% select(matches("glucose")) %>% 
    mutate(across(!matches("Prior|Post"), ~.x * 0), patient = row_number(), color = ifelse(`glucose-Prior` > 0, "DKA", "D5"))
  
  input_long <- 
    input %>% 
    mutate(patient = fct_reorder(as_factor(patient), `glucose-Prior`, .na_rm = F)) %>%
    pivot_longer(cols = !matches("patient|color"), names_to = c("lab", "time"), values_to = "result", names_sep = "-") %>% 
    mutate(time = as_factor(ifelse(is.na(time), " ", time))) %>%
    mutate(time = factor(time, levels = c("Prior", " ", "Post"))) %>%
    mutate(lab = factor(lab, levels = c("glucose")))
  
  gg_in <- 
    input_long %>% 
      na.omit() %>%
      ggplot(aes(time, result)) + 
      geom_vline(xintercept = " ", linetype = "dashed", color = "grey50") + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_line(aes(group = patient, color = patient), linewidth = 1, lineend = "round") +
      scale_color_viridis_d(begin = 0.7, end = 0.9, direction = -1, option = "B") +
      facet_wrap(~lab, scales = "free", labeller = labeller(lab = analyte_labels), nrow = 1) + 
      theme(strip.text = element_text(size = 14, face = "bold.italic"), legend.position = "none", axis.title = element_blank())
  ggsave("../Figures/deltas_fig5b.pdf", gg_in, width = 8, height = 4, dpi = 600)
  
  library(ggpubr)
  library(patchwork)
  
  gg_embed + 
    inset_element(gg_in, left = 0, right = 0.4, bottom = 0, top = 0.6)
  ggsave("../Figures/fig5.pdf", width = 8, height = 6, dpi = 600)
  
  input
 
}

makeValidationFigure <- function(contam_validation){
  
  library(ggpmisc)
  
  contam_validation <- read_delim(here::here("../Data/contamination_simulation_validation.txt"))
  
  control_mean <- 
    contam_validation %>% 
      filter(Fluid == "Control") %>% 
      summarise(across(!matches("ID|Fluid|Mixture", ignore.case = F), ~mean(.x))) %>%
      mutate(anion_gap = sodium - chloride - co2_totl)
  
  
  validation_long <- 
    contam_validation %>% 
      filter(Fluid != "Control") %>% 
      select(Fluid, Mixture, any_of(lab_strings_bmp)) %>%
      pivot_longer(!matches("Fluid|Mixture|ID", ignore.case = F), names_to = "Analyte", values_to = "Result")
  
  expected <- 
    map2(contam_validation %>% 
           filter(Fluid != "Control") %>% 
           select(Mixture) %>% pluck(1), 
         contam_validation %>% 
           filter(Fluid != "Control") %>% 
           select(Fluid) %>% pluck(1), 
         ~simulateContaminationRow(control_mean, .x, .y)) %>% 
    bind_rows() %>% 
    bind_cols(contam_validation %>% 
                filter(Fluid != "Control") %>% 
                select(Mixture, Fluid), .)
  
  expected_long <- 
    expected %>% select(-anion_gap) %>% pivot_longer(!matches("Fluid|Mixture|ID", ignore.case = F), names_to = "Analyte", values_to = "Result")
  
  formula = y ~ poly(x, degree = 1)
  
  ggplot() +
    geom_point(data = validation_long, aes(Mixture, Result, color = Fluid)) + 
    geom_line(data = expected_long, aes(Mixture, Result, color = Fluid)) +
    facet_grid(cols = vars(Fluid), rows = vars(Analyte), scales = "free", labeller = labeller(Analyte = analyte_labels)) +
    scale_color_manual(values = color_map_global) + 
    theme(legend.position = "none", strip.text = element_text(size = 12, face = "bold.italic"))
  ggsave("../Figures/validation_plots.pdf", width = 10, height = 10)
  
}
