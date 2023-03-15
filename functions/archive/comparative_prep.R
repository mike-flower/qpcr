#===============================================================================
# Comparative ct target prep
#===============================================================================
comparative_prep <- function(ct, settings) {
  
  # # Manual
  # ct = ct
  # settings = settings
  
  
  #===============================================================================
  # Average tech reps for each sample (removing outliers and exclusions)
  #===============================================================================
  summary <- lapply(ct, function(x) {
    
    # Target
    t = as.character(
      x %>%
        slice(1) %>%
        pull(target)
    )
    
    # Count outliers and exclusions
    oe <- x %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(total = n(),
                       outliers = sum(outlier == TRUE),
                       exclusions = sum(!is.na(!!sym(exclusion_var))))
    
    # Remove outliers and exclusions
    y <- x %>%
      {if (remove_outliers == "Yes")
      {dplyr::filter(., outlier != TRUE, is.na(!!sym(exclusion_var)))}
        else {dplyr::filter(., is.na(!!sym(exclusion_var)))}} # https://stackoverflow.com/questions/47624161/use-filter-in-dplyr-conditional-on-an-if-statement-in-r
    
    # Summarise
    z <- y %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise_at("ct",
                          list(ct_mean = mean,
                               sd = sd),
                          na.rm = T) %>%
      dplyr::rename(ct = ct_mean) %>%
      mutate_at(c("ct", "sd"),
                function(x) ifelse(is.nan(x), NA, x)) %>%
      mutate(target = t) %>%
      relocate(target)
    
    # Summarise number of tech reps included (has to be done separately)
    n_included <- y %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(included = n())
    
    # Merge mean, sd and n
    z <- z %>%
      left_join(oe, by = unique_id) %>%
      left_join(n_included, by = unique_id) %>%
      relocate(c(ct, sd), .after = "included")
    
    # Output
    return(z)
    
  })
  
  # Export summaries
  summary_rbind <- rbindlist(summary)
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "summary")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(summary_rbind, file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "summary", append = T, row.names = F)
  
  
  
  #===============================================================================
  # Pick calibrator samples
  #===============================================================================
  # Calibrator samples
  if(calibration_method == "Default" && 
     !is.na(calibrator_var) &&
     sum(!is.na(settings[[calibrator_var]])) > 0) {
    
    calibrators <- as.character(
      settings %>%
        filter(!is.na(!!sym(calibrator_var))) %>%
        distinct(!!sym(unique_id)) %>%
        pull(!!sym(unique_id)))
    
  } else if (calibration_method == "Default") {
    
    calibrators <- as.character(
      rbindlist(summary) %>%
        slice(which.min(ct)) %>%
        slice(1) %>%
        pull(!!sym(unique_id))
    )
    
  } else if (calibration_method == "Callibrator sample/s") {
    
    calibrators <- as.character(
      settings %>%
        filter(!is.na(!!sym(calibrator_var))) %>%
        distinct(!!sym(unique_id)) %>%
        pull(!!sym(unique_id)))
    
  } else if (calibration_method == "Lowest GOI ct") {
    
    calibrators <- as.character(
      rbindlist(summary) %>%
        slice(which.min(ct)) %>%
        slice(1) %>%
        pull(!!sym(unique_id)))
    
  } else if (calibration_method == "Highest GOI ct") {
    
    calibrators <- as.character(
      rbindlist(summary) %>%
        slice(which.max(ct)) %>%
        slice(1) %>%
        pull(!!sym(unique_id)))
    
  } else if (calibration_method == "Select sample/s") {
    
    calibrators <- calibrator_samples
    
  }
  
  calibrators
  
  
  
  #===============================================================================
  # Find HK genes where not every sample worked
  #===============================================================================
  # Find bad hk genes
  bad_hk <- as.character(
    rbindlist(summary) %>%
      select(-c(total, outliers, exclusions, included, sd)) %>%
      filter(target %in% hk,
             !!sym(unique_id) %notin% nc_samples,
             is.na(ct)) %>%
      distinct(target) %>%
      pull(target)
  )
  
  # Good hk genes
  good_hk <- setdiff(hk, bad_hk)
  
  # Target summary
  targets <- rbindlist(summary) %>%
    select(-c(total, outliers, exclusions, included, sd)) %>%
    dplyr::mutate(sample_type = ifelse(!!sym(unique_id) %in% nc_samples, "nc", "sample")) %>%
    group_by(target, sample_type) %>%
    dplyr::summarise(worked = sum(!is.na(ct)),
                     failed = sum(is.na(ct))) %>%
    pivot_wider(names_from = sample_type, values_from = c(worked, failed)) %>%
    dplyr::mutate(goi = ifelse(target %in% goi, TRUE, FALSE),
                  hk = ifelse(target %in% hk, TRUE, FALSE),
                  hk_included = ifelse(target %in% good_hk, TRUE, FALSE)) %>%
    relocate(c(goi, hk), .after = "target")
    
  # Export target summary
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "targets")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(as.data.frame(targets), file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "targets", append = T, row.names = F)
  
  # Output
  return(list(summary = summary,
              calibrators = calibrators,
              bad_hk = bad_hk,
              good_hk = good_hk,
              targets = targets))
  
}
    
  









