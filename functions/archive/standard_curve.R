#===============================================================================
# Standard curve calculations
#===============================================================================
standard_curve <- function(data) {
  
  # Manual
  data = ct
  
  
  #===============================================================================
  # Add sample_type (for plotting later)
  #===============================================================================
  data <- lapply(data, function(x) {
    
    # Extract data and add sample_type
    y <- x[["data"]] %>%
      mutate(sample_type = case_when(!is.na(.[[standard_var]]) ~ "Standard",
                                     TRUE ~ "Sample"))
    
    # Output
    x[["data"]] <- y
    return(x)
    
  })
  
  
  #===============================================================================
  # Standard curve
  #===============================================================================
  # Regression formulae
  formulae <- list("linear" = y ~ x,
                   "log" = y ~ log(x))
                   # "poly2" = y ~ poly(x, 2),
                   # "poly2raw" = y ~ poly(x, 2, raw = TRUE),
                   # "poly3" = y ~ poly(x, 3),
                   # "poly3raw" = y ~ poly(x, 3, raw = TRUE))
  
  # For each non-IPC target model using all formulae
  models <- lapply(data, function(x) {
    
    # Extract ct data
    y <- x[["data"]]
    
    # Extract standards and remove failed and outliers samples
    standards <- y %>%
      dplyr::filter(!is.na(!!sym(standard_var)),
                    !!sym(standard_var) > 0,
                    !is.na(ct),
                    outlier != TRUE)
    
    models <- lapply(setNames(names(formulae), names(formulae)), function(n) {
      
      # Extract formula
      f <- formulae[[n]]
      
      # Plot
      p <-
        ggplot(standards, aes(x = !!sym(standard_var), y = ct)) +
        geom_point() +
        stat_smooth(method = lm, formula = f, fullrange = T) +
        stat_poly_eq(formula = f,
                     aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
                     parse = TRUE, coef.digits = 3, f.digits = 3, p.digits = 3, 
                     rr.digits = 3) +
        labs(title = n,
             subtitle = c(f),
             x = standard_label,
             y = "Threshold cycle (Ct)")
      
      # Convert xy notation to real variables for lm function
      left <- gsub("y", "ct", as.character(f[2]))
      right <- gsub("x", standard_var, as.character(f[3]))
      f_long <- sprintf("%s ~ %s", left, right) # https://stackoverflow.com/questions/26381410/edit-and-reuse-the-formula-part-of-the-call-for-a-model-in-r
      f_long <- as.formula(f_long)
      f_long
      
      # Linear regression to get r2
      model <- lm(f_long, data = standards)
      r2 <- summary(model)$r.squared
      r2
      
      # Reverse the linear regression for predicting conc from ct
      f_reverse <- sprintf("%s ~ %s", right, left)
      model_reverse <- lm(f_reverse, data = standards)
      summary(model_reverse)$r.squared
      
      # Output
      return(list(plot = p,
                  model = model,
                  r2 = r2,
                  model_reverse = model_reverse))
      
    })
    
    # Output
    return(models)
    
  })
  
  
  #===============================================================================
  # Print plots of each model
  #===============================================================================
  for (t in names(models)) {
    
    # Extract plots
    plots <- lapply(models[[t]], extract_sublist, sub = "plot")
    
    # Print plots
    print(
      wrap_plots(plots) +
        plot_annotation(title = paste0(t, " standard curve"),
                        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    )
    
  }
  
  
  #===============================================================================
  # Export r2 for each model and find best model for each target
  #===============================================================================
  r2 <- lapply(setNames(names(models), names(models)), function(t) {
    
    # Extract r2
    r2 <- unlist(lapply(models[[t]], extract_sublist, sub = "r2"))
    
    # Find best model
    best_model <- r2[which.max(r2)]
    
    # r2 as a row
    r2 <- as.data.frame(t(r2)) %>%
      mutate(target = t,
             best_model = names(best_model)) %>%
      relocate(target)
    
    # Output
    return(r2)
    
  })
  
  # rbind r2
  r2 <- rbindlist(r2)
  r2
  
  # Export r2 table
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "r2")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(r2,
             file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "r2", append = T, row.names = F)
  
  
  
  
  #===============================================================================
  # Predict conc from ct
  #===============================================================================
  data <- lapply(data, function(x) {
    
    # Target name
    t <- as.character(
      x[["data"]] %>%
        slice(1) %>%
        pull(target)
    )
    
    # Extract best model
    best_model <- r2 %>%
      filter(target == t) %>%
      pull(best_model)
    
    # Extract reverse model
    model_reverse <- models[[t]][[best_model]][["model_reverse"]]
    
    # Add predictions of conc from ct
    y <- x[["data"]] %>%
      add_predictions(model_reverse, var = "predict_conc") %>%
      mutate(predict_conc = exp(predict_conc)) # https://stackoverflow.com/questions/75252211/predicted-values-from-a-log-model-using-add-predictions-are-a-long-way-off/75252264#75252264
    
    # Plot observed vs expected conc for standards
    p <- 
      y %>%
      filter(!is.na(!!sym(standard_var)),
             !!sym(standard_var) > 0,
             !is.na(ct)) %>%
      filter(outlier != TRUE) %>%
      ggplot(aes(x = !!sym(standard_var), y = predict_conc)) +
      geom_point(aes(fill = !!sym(standard_var)),
                 shape = 21, size = 5, alpha = 0.5, colour = "black") +
      geom_text_repel(aes(label = !!sym(standard_var))) +
      stat_smooth(method = lm, formula = y ~ x, fullrange = T) +
      stat_poly_eq(formula = y ~ x,
                   aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
                   parse = TRUE, coef.digits = 3, f.digits = 3, p.digits = 3, 
                   rr.digits = 3) +
      scale_x_continuous(trans='log2') +
      scale_y_continuous(trans='log2') +
      labs(title = t,
           subtitle = paste0("Predictive model: ", best_model, " [", c(formulae[[best_model]]), "]"),
           x = standard_label,
           y = "Observed concentration",
           fill = standard_label)
    
    # Output
    x[["data"]] <- y
    out <- append(x, list(plot_oe = p))
    return(out)
    
  })
  
  
  # Print observed vs expected plots
  plots <- lapply(data, extract_sublist, sub = "plot_oe")
  print(
    wrap_plots(plots) +
      plot_annotation(title = "Observed vs expected concentration (standards)",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
  
  
  
  #===============================================================================
  # Calculate total DNA in source volume
  #===============================================================================
  data <- lapply(data, function(x) {
    
    # Extract data and calculate total DNA yield
    y <- x[["data"]] %>%
      dplyr::mutate(total_DNA = predict_conc * !!sym(vol_var))
    
    # Output
    x[["data"]] <- y
    return(x)
    
  })
  
  
  
  #===============================================================================
  # Output data with predicted concentration and DNA yield
  #===============================================================================
  # Extract ct data
  predict <- rbindlist(lapply(data, extract_sublist, sub = "data"))
  
  # Export
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "predict")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(predict,
             file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "predict", append = T, row.names = F)
  
  
  #===============================================================================
  # DNA quality index (small / large)
  #===============================================================================
  # <1. DNA is not degraded
  # 1-10. Mild degredation, or PCR inhibition.
  # >10 / NA. Significant degredation
  
  # Calculate quality index
  quality <- predict %>%
    select(filename, well_short, !!sym(unique_id), sample_type, outlier, target, predict_conc) %>%
    pivot_wider(names_from = target,
                values_from = c(outlier, predict_conc)) %>%
    dplyr::mutate(quality = !!sym(paste0("predict_conc_", quality_numerator)) / !!sym(paste0("predict_conc_", quality_denominator)))
  
  # Count tech reps, excluding outliers
  quality_n <- quality %>%
    dplyr::filter(!!sym(paste0("outlier_", quality_numerator)) != TRUE,
                  !!sym(paste0("outlier_", quality_denominator)) != TRUE) %>%
    group_by(!!sym(unique_id)) %>%
    dplyr::summarise(n = n())
  
  # Summarise quality index (excluding outliers)
  quality_summary <- quality %>%
    dplyr::filter(!!sym(paste0("outlier_", quality_numerator)) != TRUE,
                  !!sym(paste0("outlier_", quality_denominator)) != TRUE) %>%
    group_by(!!sym(unique_id), sample_type) %>%
    dplyr::summarise(quality_mean = mean(quality, na.rm = TRUE),
                     quality_sd = sd(quality, na.rm = TRUE)) %>%
    mutate_at(vars(contains("quality")),
              function(x) ifelse(is.nan(x), NA, x)) %>%
    dplyr::rename(quality = quality_mean) %>%
    left_join(quality_n, by = unique_id) %>%
    ungroup() %>%
    mutate(quality_sem = quality_sd / sqrt(n)) %>%
    relocate(n, .after = !!sym(unique_id))
    
  # Sort by order in the settings file
  quality_summary <- quality_summary[match(unique(settings[[unique_id]]), quality_summary[[unique_id]]),]  
  
  # Export quality summary
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "quality")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(as.data.frame(quality_summary),
             file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "quality", append = T, row.names = F)
  
  # Plot quality index variability
  print(
    ggplot(quality_summary, aes(x = sample_type, y = quality)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width=0.1,height=0, size = 3) +
      geom_hline(yintercept = c(1, 10), linetype = "dashed", colour = "red") +
      labs(title = "Quality index",
           x = "Sample type",
           y = "Quality index") +
      theme(text = element_text(size = 20))
  )
  
  ### IN THE NEXT VERSION OF THE PROGRAM, AVERAGE THE CT TECH REPS, THEN PREDICT CONC AND YIELD OFF THAT
  
  
  
}