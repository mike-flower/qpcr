#===============================================================================
# Comparative ct calculations using multiple reference targets
#===============================================================================
# (2^dctGOI / geomean[2^dctREFS]) # https://toptipbio.com/qpcr-multiple-reference-genes/

multiple_hk <- function(data) {
  
  # # Manual
  # data = ct
  
  
  # Calculate control mean ct for each target
  control_mean_ct <- lapply(data, function(x) {
    
    return(
      x[["summary"]] %>%
        filter(!!sym(unique_id) %in% calibrators) %>%
        summarise(ct = mean(ct, na.rm = TRUE)) %>%
        pull(ct)
    )
    
  })
  control_mean_ct
  
  
  # Calculate 2^dct (relative quantity, RQ) for each target
  comparative_ct <- lapply(data, function(x) {
    
    # Extract target
    t <- as.character(
      x[["summary"]] %>%
        slice(1) %>%
        pull(target)
    )
    
    # Extract summary and calculate 2^dct (relative quantity, RQ)
    y <- x[["summary"]] %>%
      select(-c(starts_with("n_"), sd)) %>%
      dplyr::mutate(cmct = control_mean_ct[[t]]) %>%
      dplyr::mutate(dct = cmct - ct) %>%
      dplyr::mutate(twodct = 2 ^ dct) # assumes 100% amplification efficiency, https://toptipbio.com/qpcr-multiple-reference-genes/
    
    # Output
    return(y)
    
  })
  
  # rbind targets, pivot wider and calculate geomean of HK genes
  comparative_ct <- rbindlist(comparative_ct) %>%
    pivot_wider(names_from = target,
                values_from = c(ct, cmct, dct, twodct)) %>%
    dplyr::mutate(analysis_method = analysis_method,
                  calibration_method = calibration_method,
                  calibrator = !!sym(unique_id) %in% calibrators) %>%
    relocate(all_of(c("analysis_method", "calibration_method", "calibrator")), .after = unique_id) %>%
    dplyr::mutate(calibrator = factor(calibrator, levels = c(TRUE, FALSE))) %>%
    rowwise() %>%
    dplyr::mutate(geomean_hk = case_when(exclude_failed_hk == "No" ~
                                           exp(mean(log(c_across(c(paste0("twodct_", hk)))), na.rm = TRUE)),
                                         TRUE ~
                                           exp(mean(log(c_across(c(paste0("twodct_", good_hk)))))))) %>%
    ungroup()
    
    
  # For each GOI calculate relative expression
  for (g in goi) {
    
    comparative_ct <- comparative_ct %>%
      dplyr::mutate(!!sym(g) := !!sym(paste0("twodct_", g)) / geomean_hk)
    
  }
  
  # Output
  return(comparative_ct)
  
}