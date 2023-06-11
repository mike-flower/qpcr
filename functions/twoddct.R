#===============================================================================
# Comparative ct calculations using classic 2^-ddct method
#===============================================================================
# https://toptipbio.com/delta-delta-ct-pcr/

twoddct <- function(summary) {
  
  # Calculate control GEOmean ct for each SAMPLE
  # THE CLASSIC METHOD REQUIRES ONLY 1 HK GENE IS USED, BUT THIS IS A BODGE TO ALLOW >1 HK GENE
  if(exclude_failed_hk) {
    
    geomean_hk <- rbindlist(summary) %>%
      select(-c(total, tech_outliers, exclusions, included, sd)) %>%
      filter(target %in% good_hk,
             !!sym(unique_id) %notin% empty_samples) %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(geomean_hk = exp(mean(log(ct)))) %>%
      ungroup() 
    
  } else {
    
    geomean_hk <- rbindlist(summary) %>%
      select(-c(total, tech_outliers, exclusions, included, sd)) %>%
      filter(target %in% hk,
             !!sym(unique_id) %notin% empty_samples) %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(geomean_hk = exp(mean(log(ct), na.rm = T))) %>%
      ungroup()
    
  }
  
  
  # 2^-ddct calculations
  comparative_ct <- lapply(summary[names(summary) %in% goi], function(a) {
    
    # Extract target
    t <- as.character(
      a %>%
        slice(1) %>%
        pull(target)
    )
    
    # Split into list to calculate by calibrator group
    b <- a %>%
      split(f = as.factor(.[[calibrator_group_var]]))
    
    d <- lapply(b, function(c) {
      
      # Unlist calibrators
      cals <- unlist(calibrators, use.names = FALSE)
      
      # Calculate dct
      d <- c %>%
        select(-c(total, tech_outliers, exclusions, included, sd)) %>%
        left_join(geomean_hk, by = unique_id) %>%
        dplyr::mutate(dct = ct - geomean_hk) %>%
        dplyr::mutate(analysis_method = analysis_method,
                      calibration_method = calibration_method,
                      calibrator = !!sym(unique_id) %in% cals) %>%
        relocate(all_of(c("analysis_method", "calibration_method", "calibrator")), .after = unique_id) %>%
        dplyr::mutate(calibrator = factor(calibrator, levels = c(TRUE, FALSE)))
      
      # Calculate average dct of calibrator samples
      dct_control_mean <- d %>%
        filter(!!sym(unique_id) %in% cals) %>%
        dplyr::summarise(dct = mean(dct)) %>%
        pull(dct)
      dct_control_mean <- ifelse(is.nan(dct_control_mean), NA, dct_control_mean)
      
      # Calculate 2^-ddct
      d <- d %>%
        dplyr:: mutate(dct_control_mean = dct_control_mean) %>%
        dplyr::mutate(ddct = dct - dct_control_mean) %>%
        dplyr::mutate(!!sym(t) := 2 ^ -ddct)
      
      # Output
      return(d)
        
    })
    
    # Recombine the list of calibrator groups
    d <- rbindlist(d)
    
    # Put d samples back in their original order
    d <- d[match(unique(a[[unique_id]]), d[[unique_id]]),]
    
    # Output
    return(d)
    
    
    # # Extract summary and calculate dct
    # y <- x %>%
    #   select(-c(total, tech_outliers, exclusions, included, sd)) %>%
    #   left_join(geomean_hk, by = unique_id) %>%
    #   dplyr::mutate(dct = ct - geomean_hk) %>%
    #   dplyr::mutate(analysis_method = analysis_method,
    #                 calibration_method = calibration_method,
    #                 calibrator = !!sym(unique_id) %in% calibrators) %>%
    #   relocate(all_of(c("analysis_method", "calibration_method", "calibrator")), .after = unique_id) %>%
    #   dplyr::mutate(calibrator = factor(calibrator, levels = c(TRUE, FALSE)))
    # 
    # 
    # # Calculate average dct of calibrator samples
    # dct_control_mean <- y %>%
    #   filter(!!sym(unique_id) %in% calibrators) %>%
    #   dplyr::summarise(dct = mean(dct)) %>%
    #   pull(dct)
    # 
    # # Calculate 2^-ddct
    # z <- y %>%
    #   dplyr:: mutate(dct_control_mean = dct_control_mean) %>%
    #   dplyr::mutate(ddct = dct - dct_control_mean) %>%
    #   dplyr::mutate(!!sym(t) := 2 ^ -ddct)
    # 
    # # Output
    # return(z)
    
  })
  
  
  # rbind comparative ct calculations
  comparative_ct <- rbindlist(comparative_ct, fill = TRUE)
  
  # Output
  return(comparative_ct)
  
}
