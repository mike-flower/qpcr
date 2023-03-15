#===============================================================================
# Comparative ct calculations using classic 2^-ddct method
#===============================================================================
# https://toptipbio.com/delta-delta-ct-pcr/

twoddct <- function(data) {
  
  # # Manual
  # data = ct
  
  # Calculate control GEOmean ct for each SAMPLE
  # THE CLASSIC METHOD REQUIRES ONLY 1 HK GENE IS USED, BUT THIS IS A BODGE TO ALLOW >1 HK GENE
  if(exclude_failed_hk == "Yes") {
    
    geomean_hk <- rbindlist(lapply(data, extract_sublist, sub = "summary")) %>%
      select(-c(starts_with("n_"), sd)) %>%
      filter(target %in% good_hk,
             !!sym(unique_id) %notin% nc_samples) %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(geomean_hk = exp(mean(log(ct)))) %>%
      ungroup() 
    
  } else {
    
    geomean_hk <- rbindlist(lapply(data, extract_sublist, sub = "summary")) %>%
      select(-c(starts_with("n_"), sd)) %>%
      filter(target %in% hk,
             !!sym(unique_id) %notin% nc_samples) %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(geomean_hk = exp(mean(log(ct), na.rm = T))) %>%
      ungroup()
    
  }
  
  # 2^-ddct calculations
  comparative_ct <- lapply(data[names(data) %in% goi], function(x) {
    
    # # Manual
    # x = data[[1]]
    
    # Extract target
    t <- as.character(
      x[["summary"]] %>%
        slice(1) %>%
        pull(target)
    )
    
    # Extract summary and calculate dct
    y <- x[["summary"]] %>%
      select(-c(starts_with("n_"), sd)) %>%
      left_join(geomean_hk, by = unique_id) %>%
      dplyr::mutate(dct = ct - geomean_hk) %>%
      dplyr::mutate(analysis_method = analysis_method,
                    calibration_method = calibration_method,
                    calibrator = !!sym(unique_id) %in% calibrators) %>%
      relocate(all_of(c("analysis_method", "calibration_method", "calibrator")), .after = unique_id) %>%
      dplyr::mutate(calibrator = factor(calibrator, levels = c(TRUE, FALSE)))
    
    
    # Calculate average dct of calibrator samples
    dct_control_mean <- y %>%
      filter(!!sym(unique_id) %in% calibrators) %>%
      dplyr::summarise(dct = mean(dct)) %>%
      pull(dct)
    
    # Calculate 2^-ddct
    z <- y %>%
      dplyr:: mutate(dct_control_mean = dct_control_mean) %>%
      dplyr::mutate(ddct = dct - dct_control_mean) %>%
      dplyr::mutate(!!sym(t) := 2 ^ -ddct)
    
    # Output
    return(z)
    
  })
  
  # rbind comparative ct calculations
  comparative_ct <- rbindlist(comparative_ct)
  
  # Output
  return(comparative_ct)
  
}


