# Global settings and variables

theme_ns <- theme(text = element_text(family = "Helvetica"),
                  title = element_text(face = "bold.italic", size = 14),
                  plot.background = element_rect(fill = "#ffffff"),
                  plot.subtitle = element_text(size = 12, face = "plain", hjust = 0),
                  plot.title = element_text(hjust = 0),
                  axis.title = element_text(size = 12, face = "bold", margin = margin(4,4,4,4)),
                  axis.title.x.bottom = element_text(margin = margin(8,0,0,0)),
                  axis.title.y.left = element_text(margin = margin(0,8,0,0)),
                  legend.title = element_text(face = "bold.italic", size = 12),
                  axis.line = element_line(),
                  axis.ticks = element_blank(),
                  panel.grid = element_blank(), 
                  panel.background = element_blank(),
                  strip.text = element_text(size = 10, face = "bold.italic"),
                  strip.background = element_blank())
theme_set(theme_ns)

assign("libs", .libPaths(), envir = globalenv())

assign("color_map_global",
       c("NS" = scico::scico(palette = "roma", direction = -1, n = 17)[15],
         "LR" = scico::scico(palette = "roma", direction = -1, n = 17)[12],
         "halfNS" = scico::scico(palette = "roma", direction = -1, n = 17)[14],
         "SW" = scico::scico(palette = "roma", direction = -1, n = 17)[10],
         "D5halfNS" = scico::scico(palette = "roma", direction = -1, n = 17)[7],
         "D5halfNSwK" = scico::scico(palette = "roma", direction = -1, n = 17)[6],
         "hyperNS" = scico::scico(palette = "roma", direction = -1, n = 17)[16],
         "D5NS" = scico::scico(palette = "roma", direction = -1, n = 17)[3],
         "D5LR" = scico::scico(palette = "roma", direction = -1, n = 17)[4],
         "D5W" = scico::scico(palette = "roma", direction = -1, n = 17)[1]),
       envir = globalenv())

assign("analyte_labels",
       c("glucose" = "Glucose",
         "sodium" = "Sodium",
         "chloride" = "Chloride",
         "potassium_plas" = "Potassium",
         "calcium" = "Calcium",
         "co2_totl" = "CO2",
         "bun" = "BUN",
         "creatinine" = "Creatinine",
         "anion_gap" = "Anion Gap"),
       envir = globalenv())

assign("metric_labels",
       c("classification_cost_2to1" = "Cost",
         "binary_cost" = "Binary Cost",
         "mcc" = "MCC",
         "pr_auc" = "auPRC",
         "roc_auc" = "auROC",
         "sensitivity" = "Sens",
         "specificity" = "Spec",
         "ppv" = "PPV",
         "npv" = "NPV",
         "accuracy" = "Accuracy"),
       envir = globalenv())

assign("lab_strings_bmp", c("sodium", "chloride", "bun", "calcium", "potassium_plas", "co2_totl", "creatinine", "glucose", "anion_gap"), envir = globalenv())
assign("lab_strings_bmp_no_gap", c("sodium", "chloride", "bun", "calcium", "potassium_plas", "co2_totl", "creatinine", "glucose"), envir = globalenv())
assign("lab_strings_cmp", c("sodium", "chloride", "potassium_plas", "co2_totl", "bun", "creatinine", "calcium", "glucose", "anion_gap", "albumin", "ast", "alt", "bili_totl", "protein_plas", "alk_phos"), envir = globalenv())
assign("lab_strings_cbc", c("wbc", "hgb", "hct", "plt", "rbc"), envir = globalenv())
assign("lab_strings", c(lab_strings_cmp, lab_strings_cbc), envir = globalenv())

fluids <- setNames(rep(0, length(lab_strings)), lab_strings)

SW <- fluids

NS <- fluids
NS[c("sodium", "chloride")] <- c(154, 154)

halfNS <- NS/2
hyperNS <- NS*3/(0.9)

D5NS <- NS
D5NS["glucose"] <- 5000

D5halfNS <- halfNS
D5halfNS['glucose'] <- 5000

D5quarterNS <- halfNS/2
D5quarterNS['glucose'] <- 5000

D5halfNSwK <- D5halfNS
D5halfNSwK[c("potassium_plas", "chloride")] <- c(20, 97)

LR <- fluids
LR[c("sodium", "chloride", "potassium_plas", "calcium", "anion_gap")] <- c(130, 109, 4, 5.4, 21)

D5LR <- LR
D5LR["glucose"] <- 5000

D5W <- fluids
D5W["glucose"] <- 5000


#assign("fluids", list(NS, LR, D5W), envir = globalenv())
#assign("fluid_names", c("NS", "LR", "D5W"), envir = globalenv())

#assign("fluids", list(NS, D5NS, LR, D5LR, D5W), envir = globalenv())
#assign("fluid_names", c("NS", "D5NS", "LR", "D5LR", "D5W"), envir = globalenv())

assign("fluids", list(NS, D5NS, LR, D5LR, D5W, D5halfNSwK, D5halfNS, halfNS, hyperNS, SW), envir = globalenv())
assign("fluid_names", c("NS", "D5NS", "LR", "D5LR", "D5W", "D5halfNSwK", "D5halfNS", "halfNS", "hyperNS", "SW"), envir = globalenv())
names(fluids) <- fluid_names

relabeller <- setNames(names(analyte_labels), analyte_labels)

panels = c("BMP")
train_cohorts = c("BJH")

task_assays_bmp <- c("Sodium", "Chloride", "Potassium Plas", "CO2 Totl", "BUN", "Creatinine", "Calcium", "Glucose", "Anion Gap")

# Means of the high and low level CV from recent QC
CV <- list("albumin" = (3.45+4.44)/2, "alk_phos" = (4.55+3.24)/2, "alt" = (10.71+4.71)/2, "ast" = (8.33+2.90)/2,
                 "bili_totl" = (12.50+4.35)/2, "bun" = (7.14+4.65)/2, "calcium" = (1.98+1.79)/2, "chloride" = (1.98+2.35)/2,
                 "co2_totl" = (7.98+6.25)/2, "creatinine" = (5.10+3.54)/2, "glucose" = (3.53+2.50)/2, "potassium_plas" = (2.38+1.52)/2,
                 "protein_plas" = (2.08+2.86)/2, "sodium" = (1.22+1.28)/2, "anion_gap" = 1)

# Calculated using https://biologicalvariation.eu/ at the mean CVs of the levels in QC. 
RCV_increase <- list(sodium = 0.032, chloride = 0.058, potassium_plas = 0.11, co2_totl = 0.208, creatinine = 0.155, bun = 0.418, calcium = 0.062, glucose = 0.145)
RCV_decrease <- list(sodium = 0.031, chloride = 0.055, potassium_plas = 0.099, co2_totl = 0.172, creatinine = 0.134, bun = 0.295, calcium = 0.059, glucose = 0.126)

# Aggregated from https://datainnovations.com/allowable-total-error-table 
TAE_absolute <- list(sodium = 4, chloride = 3, potassium_plas = 0.5, co2_totl = 2, creatinine = 0.3, bun = 2, calcium = 1, glucose = 6)
TAE_percent <- list(sodium = 0.05, chloride = 0.05, potassium_plas = 0.058, co2_totl = 0.1, creatinine = 0.15, bun = 0.09, calcium = 0.024, glucose = 0.1)

TAE_CLIA <- list(sodium = list(threshold = 4, type = "absolute", source = "CLIA"), 
                 chloride = list(threshold = 3, type = "absolute", source = "RPCA"), 
                 potassium_plas = list(threshold = 0.5, type = "absolute", source = "CLIA"), 
                 co2_totl = list(threshold = 2, type = "absolute", source = "RPCA"),
                 creatinine = list(threshold = 0.3, type = "absolute", source = "CLIA"),
                 bun = list(threshold = 2, type = "absolute", source = "CLIA"),
                 calcium = list(threshold = 1, type = "absolute", source = "CLIA"),
                 glucose = list(threshold = 6, type = "absolute", source = "CLIA"),
                 link = "https://datainnovations.com/allowable-total-error-table")

analyte_ranges <- 
  list(
    sodium = list(min = 120, max = 180),
    chloride = list(min = 80, max = 140),
    potassium_plas = list(min = 0, max = 8), 
    co2_totl = list(min = 0, max = 35),
    creatinine = list(min = 0, max = 3),
    bun = list(min = 0, max = 50),
    calcium = list(min = 0, max = 10),
    glucose = list(min = 0, max = 1000)
  )

critical_ranges <- 
  list(
    Sodium = list(min = 120, max = 160),
    Chloride = list(min = 75, max = 130),
    `Potassium Plas` = list(min = 2.5, max = 6.2),
    `CO2 Totl` = list(min = 10, max = 45),
    Calcium = list(min = 6.5, max = 14),
    Glucose = list(min = 54, max = 450)
  )
