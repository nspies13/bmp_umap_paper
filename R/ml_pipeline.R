makeMLinputs <- function(bmp_no_NA_verified, contam_sim){
  
  library(foreach)
  
  foreach(ratio = seq(0.01, 0.5, by = 0.01)) %dopar% {
    
  }
  simulateContaminationRow(bmp_no_NA_verified[10,], mix_ratio, fluid)
  
}