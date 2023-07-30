# Global settings and variables
theme_ns <- theme(text = element_text(family = "Helvetica"),
                  title = element_text(face = "bold", size = 18),
                  plot.title = element_text(hjust = 0.5),
                  axis.title = element_text(size = 24, face = "bold.italic"),
                  legend.title = element_text(face = "bold.italic", size = 12),
                  axis.line = element_line(),
                  axis.ticks = element_blank(),
                  panel.grid = element_blank(),
                  panel.background = element_blank(),
                  strip.text = element_text(size = 18, face = "bold.italic"),
                  strip.background = element_blank())
theme_set(theme_ns)

assign("libs", .libPaths(), envir = globalenv())

assign("color_map_global",
       c("NS" = viridis::viridis(15)[14],
         "LR" = viridis::viridis(15)[10],
         "Non-D5" = viridis::viridis(15)[12],
         "D5NS" = viridis::viridis(15)[6],
         "D5LR" = viridis::viridis(15)[2],
         "D5W" = viridis::viridis(15)[4],
         "D5" = viridis::viridis(15)[5]),
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

model_board <- pins::board_folder("../Results/Models/", versioned = T)

assign("lab_strings_bmp", c("sodium", "chloride", "potassium_plas", "co2_totl", "bun", "creatinine", "calcium", "glucose", "anion_gap"), envir = globalenv())

NS <- as_tibble_row(c(154, 154, 0, 0, 0, 0, 0, 0, 0), .name_repair = "minimal") %>% set_names(lab_strings_bmp)
D5NS <- as_tibble_row(c(154, 154, 0, 0, 0, 0, 0, 5000, 0), .name_repair = "minimal") %>% set_names(lab_strings_bmp)
LR <- as_tibble_row(c(130, 109, 4, 0, 0, 0, 5.4, 0, 21), .name_repair = "minimal") %>% set_names(lab_strings_bmp)
D5LR <- as_tibble_row(c(130, 109, 4, 0, 0, 0, 5.4, 5000, 21), .name_repair = "minimal") %>% set_names(lab_strings_bmp)

assign("fluids", list(NS, D5NS, LR, D5LR), envir = globalenv())
assign("fluid_vars", list("NS" = c("sodium", "chloride", "calcium", "potassium_plas", "co2_totl"),
                          "D5NS" = c("sodium", "chloride", "calcium", "potassium_plas", "co2_totl", "glucose"),
                          "LR" = c("sodium", "chloride", "calcium", "co2_totl"),
                          "D5LR" = c("sodium", "chloride", "calcium", "co2_totl", "glucose")))
assign("fluid_names", c("NS", "D5NS", "LR", "D5LR"), envir = globalenv())

CV <- tibble_row("albumin" = (3.45+4.44)/2, "alk_phos" = (4.55+3.24)/2, "alt" = (10.71+4.71)/2, "ast" = (8.33+2.90)/2,
                 "bili_totl" = (12.50+4.35)/2, "bun" = (7.14+4.65)/2, "calcium" = (1.98+1.79)/2, "chloride" = (1.98+2.35)/2,
                 "co2_totl" = (7.98+6.25)/2, "creatinine" = (5.10+3.54)/2, "glucose" = (3.53+2.50)/2, "potassium_plas" = (2.38+1.52)/2,
                 "protein_plas" = (2.08+2.86)/2, "sodium" = (1.22+1.28)/2, "anion_gap" = 1)
