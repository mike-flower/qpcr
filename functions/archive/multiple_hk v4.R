#===============================================================================
# Comparative ct calculations using multiple reference targets
#===============================================================================
# (2^dctGOI / geomean[2^dctREFS]) # https://toptipbio.com/qpcr-multiple-reference-genes/

multiple_hk <- function(ct_clean, summary) {
  
  # Manual
  ct_clean = ct_clean
  summary = summary
  
  
  # Calculate control mean ct for each target
  # (use summary rather than ct_clean, otherwise samples with more reps will contribute to mean ct more)
  control_mean_ct <- lapply(summary, function(x) {
    y <- x %>%
      filter(!!sym(unique_id) %in% calibrators) %>%
      summarise(ct = mean(ct, na.rm = TRUE)) %>%
      pull(ct)
    return(y)
  })
  control_mean_ct
  
  
  # Calculate 2^dct (relative quantity, RQ) for each target
  comparative_ct <- lapply(ct_clean, function(x) {
    
    # Extract target
    t <- as.character(
      x %>%
        slice(1) %>%
        pull(target)
    )
    
    # Extract summary and calculate 2^dct (relative quantity, RQ)
    y <- x %>%
      select(filename, well_short, target, !!sym(unique_id), ct) %>%
      dplyr::mutate(cmct = control_mean_ct[[t]]) %>%
      dplyr::mutate(dct = cmct - ct) %>%
      dplyr::mutate(twodct = 2 ^ dct) # assumes 100% amplification efficiency, https://toptipbio.com/qpcr-multiple-reference-genes/
    
    # Output
    return(y)
    
  })
  
  
  # rbind targets, pivot wider and calculate geomean of HK genes
  temp <- rbindlist(comparative_ct) %>%
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
  
  temp2 <- temp %>%
    filter(well_short == "A1")
  
  
  # For each GOI calculate relative expression
  for (g in goi) {
    
    comparative_ct <- comparative_ct %>%
      dplyr::mutate(!!sym(g) := !!sym(paste0("twodct_", g)) / geomean_hk)
    
  }
  
  # Output
  return(comparative_ct)
  
}



