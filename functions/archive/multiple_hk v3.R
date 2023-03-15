#===============================================================================
# Comparative ct calculations using multiple reference targets
#===============================================================================
# (2^dctGOI / geomean[2^dctREFS]) # https://toptipbio.com/qpcr-multiple-reference-genes/
# NEED RELATIVE EXPRESSION OF GOI FOR EVERY SAMPLE, SO CAN'T START FROM THE SUMMARY, NEED TO USE CT DATA

multiple_hk <- function(summary) {
  
  # Manual
  ct = ct
  
  # Remove outliers and exclusions
  ct_clean <- lapply(ct, function(x) {
    
    y <- x %>%
      {if (remove_outliers == "Yes")
      {dplyr::filter(., outlier != TRUE, is.na(!!sym(exclusion_var)))}
        else {dplyr::filter(., is.na(!!sym(exclusion_var)))}} # https://stackoverflow.com/questions/47624161/use-filter-in-dplyr-conditional-on-an-if-statement-in-r
    
    # Output
    return(y)
    
  })
  
  
  # Calculate control mean ct for each target
  control_mean_ct <- lapply(ct_clean, function(x) {
    
    y <- x %>%
      filter(!!sym(unique_id) %in% calibrators) %>%
      summarise(ct = mean(ct, na.rm = TRUE)) %>%
      pull(ct)
    return(y)
  })
  control_mean_ct
  
  
  # Calculate 2^dct (relative quantity, RQ) for each target
  comparative_ct <- lapply(ct_clean, function(x) {
    
    # Manual
    x = ct_clean[[1]]
    
    # Extract target
    t <- as.character(
      x %>%
        slice(1) %>%
        pull(target)
    )
    
    # Calculate 2^dct (relative quantity, RQ)
    y <- x %>%
      select(-c(total, outliers, exclusions, included, sd)) %>%
      dplyr::mutate(cmct = control_mean_ct[[t]]) %>%
      dplyr::mutate(dct = cmct - ct) %>%
      dplyr::mutate(twodct = 2 ^ dct) # assumes 100% amplification efficiency, https://toptipbio.com/qpcr-multiple-reference-genes/
    
    # Output
    return(y)
    
  })
  
  
  
}
  