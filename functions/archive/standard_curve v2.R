#===============================================================================
# Standard curve calculations
#===============================================================================
standard_curve <- function(ct) {
  
  # Manual
  ct = ct
  
  
  #===============================================================================
  # Add sample_type (for plotting later)
  #===============================================================================
  y <- lapply(ct, function(x) {
    
    # Extract data and add sample_type
    x <- x %>%
      mutate(sample_type = case_when(!is.na(.[[standard_var]]) ~ "Standard",
                                     TRUE ~ "Sample"))
    
    # Output
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
  
  # Model each target using all formulae
  models <- lapply(y, function(z) {
    
    # Extract standards and remove failed and outliers samples
    standards <- z %>%
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
             x = paste0("Standard concentration (", conc_unit, ")"),
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
  # Plots of each model
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
  y <- lapply(y, function(z) {
    
    # Target name
    t <- as.character(
      z %>%
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
    z <- z %>%
      dplyr::mutate(best_model = best_model) %>%
      add_predictions(model_reverse, var = "predict_conc") %>%
      mutate(predict_conc = exp(predict_conc)) # https://stackoverflow.com/questions/75252211/predicted-values-from-a-log-model-using-add-predictions-are-a-long-way-off/75252264#75252264
    
    # Output
    return(z)
    
  })
  
  
  
  # Plot observed vs expected conc for standards
  plot <- lapply(y, function(z) {
    
    # Target name
    t <- as.character(
      z %>%
        slice(1) %>%
        pull(target)
    )
    
    # Extract best model
    best_model <- as.character(
      z %>%
        slice(1) %>%
        pull(best_model)
    )
    
    # Plot
    p <- 
      z %>%
      filter(!is.na(!!sym(standard_var)),
             !!sym(standard_var) > 0,
             !is.na(ct),
             outlier != TRUE) %>%
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
           x = paste0("Standard concentration (", conc_unit, ")"),
           y = "Observed concentration",
           fill = "Standard concentration")
        
    # Output
    return(p)
    
  })
  
  
  # Print observed vs expected plots
  print(
    wrap_plots(plot[2:3]) +
      plot_annotation(title = "Observed vs expected concentration (standards)",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
  
  
  
  #===============================================================================
  # Calculate total DNA in source volume
  #===============================================================================
  y <- lapply(y, function(z) {
    
    # Extract data and calculate total DNA yield
    z <- z %>%
      dplyr::mutate(total_DNA = predict_conc * !!sym(vol_var))
    
    # Output
    return(z)
    
  })
  
  
  
  #===============================================================================
  # Output data with predicted concentration and DNA yield
  #===============================================================================
  # Extract ct data
  predict <- rbindlist(y)
  
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
  
  if (quantifiler == "Yes") {
  
    # Calculate quality index
    quality <- rbindlist(y) %>%
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
        geom_label_repel(data = quality_summary %>% top_frac(0.10, quality),
                         aes(label = !!sym(unique_id)),
                         alpha = 0.5) +
        geom_hline(yintercept = c(1, 10), linetype = "dashed", colour = "red") +
        labs(title = "Quality index",
             x = "Sample type",
             y = "Quality index") +
        theme(text = element_text(size = 20))
    )
    
  }
  
  
  
  #===============================================================================
  # Summarise results for each unique sample
  #===============================================================================
  summary <- lapply(y, function(z){
    
    # Target name
    t <- as.character(
      z %>%
        slice(1) %>%
        pull(target)
    )
    
    # Summarise number of reps per unique sample (excluding failed samples and outliers)
    tech_n <- z %>%
      dplyr::filter(!is.na(ct),
                    outlier != TRUE) %>%
      group_by(!!sym(unique_id)) %>%
      dplyr::summarise(n = n())
    
    # Summarise data columns for each unique sample
    z <- z %>%
      dplyr::filter(!is.na(ct),
                    outlier != TRUE) %>%
      select(filename, !!sym(unique_id), target, ct, sample_type, predict_conc, total_DNA) %>%
      dplyr::mutate(across(c("ct", "predict_conc", "total_DNA"), ~as.numeric(.))) %>%
      group_by(!!sym(unique_id), sample_type) %>%
      dplyr::summarise_if(is.numeric, 
                          list(mean = mean, sd = sd),
                          na.rm = T) %>%
      mutate_if(is.numeric, function(x) {ifelse(is.nan(x), NA, x)}) %>%
      left_join(tech_n, by = unique_id) %>%
      mutate(ct_sem = ct_sd / sqrt(n),
             predict_conc_sem = predict_conc_sd / sqrt(n),
             total_DNA_sem = total_DNA_sd / sqrt(n),
             target = t) %>%
      relocate(c(target, n), .after = sample_type) %>%
      rename_at(.vars = vars(ends_with("_mean")),
                .funs = funs(sub("_mean", "", .)))
    
    # Output
    return(as.data.frame(z))
    
  })
  
  # Rbind the list into one table and add quality index (done separately as single value for sample, not different for each target)
  summary_rbind <- rbindlist(summary) %>%
    pivot_wider(names_from = target,
                values_from = -all_of(c(unique_id, "sample_type", "target")))
    # left_join(settings %>%
    #             select(all_of(c(unique_id, vol_var, standard_var))) %>%
    #             distinct(),
    #           by = unique_id) %>%
    # left_join(quality_summary %>% select(-c(n, sample_type)), by = unique_id) %>%
    # relocate(all_of(c(vol_var, standard_var)),
    #          .after = "sample_type")
  
  # Sort samples by the order in the settings file
  summary_rbind <- summary_rbind[match(unique(settings[[unique_id]]), summary_rbind[[unique_id]]),]
  
  # Export
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "summary")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(as.data.frame(summary_rbind),
             file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "summary", append = T, row.names = F)
  
  
  
  
  #===============================================================================
  # Plot concentration summary
  #===============================================================================
  print(
    rbindlist(summary) %>%
      ggplot(aes(x = target, y = predict_conc, fill = sample_type)) +
      geom_violin(trim = TRUE, position = position_dodge(width = 0.7)) +
      geom_boxplot(position = position_dodge(width = 0.7), width = 0.1) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                                 dodge.width = 0.7),
                 shape = 21, size = 3) +
      scale_y_continuous(trans="log2",
                         # expand = expansion(mult = c(NA, 0.2))
                         labels = function(z) {formatC(z, format = "e", digits = 1)} ) +
      labs(title = "Concentration",
           x = "Target",
           y = paste0("Concentration (", conc_unit, ")"),
           fill = "Sample type") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            # text = element_text(size = 15)
      )
  )
    
  
  
  #===============================================================================
  # Plot DNA yield summary
  #===============================================================================
  rbindlist(summary) %>%
    ggplot(aes(x = target, y = total_DNA)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1) +
    geom_jitter(width=0.1, height=0, size = 3) +
    scale_y_continuous(trans="log2",
                       # expand = expansion(mult = c(NA, 0.2))
                       labels = function(z) {formatC(z, format = "e", digits = 1)} ) +
    labs(title = "DNA yield",
         x = "Target",
         y = paste0("DNA yield (",
                    gsub("/.*$", "", conc_unit),
                    ")")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          # text = element_text(size = 15)
    )
  
  
  #===============================================================================
  # Plot concentration of each sample
  #===============================================================================
  # NEED TO PLOT EACH TARGET SEPARATELY AND THEN WRAP PLOT
  # WILL NEED STAT SUMMARY AS THERE ARE THREE TECH REPS FOR EACH SAMPLE
  
  plot <- lapply(y, function(z) {
    
    # Manual
    z = y[["Small"]]
    
    # Target name
    t <- as.character(
      z %>%
        slice(1) %>%
        pull(target)
    )
    
    # Plot
    p <- 
      ggplot(z, aes(x = !!sym(unique_id), y = predict_conc)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1) +
      geom_jitter(width=0.1, height=0, size = 3) +
      stat_summary(fun = mean, colour = "darkred", geom = "point", 
                   shape = 18, size = 4, show.legend = FALSE) + 
      stat_summary(fun = mean, colour = "darkred", geom = "text", show.legend = FALSE, 
                   vjust = -0.7,
                   aes(label = after_stat(ifelse(y > 0.01,
                                                 round(y, 2),
                                                 formatC(y, format = "e", digits = 1))))) +
      scale_y_continuous(trans="log2",
                         # expand = expansion(mult = c(NA, 0.2))
                         labels = function(z) {formatC(z, format = "e", digits = 1)} ) +
      labs(title = t,
           x = unique_id_label,
           y = paste0("Concentration (", conc_unit, ")")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            # text = element_text(size = 15)
      )
    
  })
  
  
  
  rbindlist(y) %>%
    ggplot(aes(x = !!sym(unique_id), y = predict_conc)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1) +
    geom_jitter(width=0.1, height=0, size = 3) +
    scale_y_continuous(trans="log2",
                       # expand = expansion(mult = c(NA, 0.2))
                       labels = function(z) {formatC(z, format = "e", digits = 1)} )
  
}