#===============================================================================
# qPCR stats
#===============================================================================
qpcr_stats <- function(expression, groups) {
  
  # # Manual
  # expression = expression
  # groups = grouping_vector
  # dependent_vars = dependent_vars
  
  
  #===============================================================================
  # Select settings data to annotate
  #===============================================================================
  s <- settings %>%
    select(any_of(c(unique_id, unname(groups)))) %>%
    distinct()
  
  
  #===============================================================================
  # Remove NA, annotate with settings data and set factor levels
  #===============================================================================
  expression_clean <- lapply(expression, function(x) {
    
    # Manual
    x = expression[["Small"]]
    
    # Remove NA
    y <- x %>%
      filter(!is.na(!!sym(unique_id)),
             !is.na(level))
    
    # Add settings
    y <- y %>%
      left_join(s, by = unique_id) %>%
      relocate(names(s))
    
    # Iterate through grouping variables removing NA and adding factor levels
    for (g in groups) {
      
      y <- y %>%
        filter(!is.na(!!sym(g))) %>%
        dplyr::mutate(!!sym(g) := factor(!!sym(g), levels = unique(settings[[g]]))) %>%
        droplevels()
      
    }
    
    # Output
    return(y)
    
  })
  
  
  #===============================================================================
  # t-test
  #===============================================================================
  # t-test each target
  plot <- lapply(expression_clean, function(x) {
    
    # Extract target
    t <- as.character(
      y %>%
        slice(1) %>%
        pull(target)
    )
    t
    
    # t-test
    ttest <- y %>%
      t_test(as.formula(paste0("level", " ~ ", groups[[1]]))) %>%
      add_y_position()
    ttest
    
    # Export t-test
    wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
    removeSheet(wb, sheetName = t)
    saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
    rm(wb)
    write.xlsx(data.frame(ttest) %>%
                 mutate(groups = as.character(groups)),
               file = file.path(out.dir, "qpcr.xlsx"),
               sheetName = t, append = T, row.names = F)
    
    # Plot
    p <-
      ggplot(y, aes(x=!!sym(groups[[1]]), y=level)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(outlier.shape = NA, width = 0.1) +
      geom_jitter(width = 0.1, height = 0, size = 3) +
      stat_pvalue_manual(ttest, 
                         label = ifelse("p.adj.signif" %in% names(ttest), "p.adj.signif", "p"),
                         tip.length = 0.01, hide.ns = T) +
      stat_summary(fun = mean, colour = "darkred", geom = "point", 
                   shape = 18, size = 4, show.legend = FALSE) + 
      stat_summary(fun = mean, colour = "darkred", geom = "text", show.legend = FALSE, 
                   vjust = -0.7,
                   aes(label = after_stat(round(y, 2)))) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = t,
           x = names(groups)[1],
           y = "Relative expression") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            # text = element_text(size = 15)
      )
    
    # Output
    return(p)
    
  })
  
  # Wrap plots
  print(
    wrap_plots(plot) +
      plot_annotation(title = "Relative expression",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
 
  
  
  #===============================================================================
  # ggstatsplot
  #===============================================================================
  # t-test each target
  plot <- lapply(expression_clean, function(x) {
    
    # # Manual
    # x = expression_clean[[1]]
    
    # Extract target
    t <- as.character(
      y %>%
        slice(1) %>%
        pull(target)
    )
    t
    
    # ggstatsplot
    # https://stackoverflow.com/questions/75147041/ggstatsplot-seems-to-think-my-data-has-a-differing-number-of-rows-but-it-doesn
    p <-
      tryCatch(
        ggbetweenstats(data = x,
                       x = !!sym(groups[[1]]), xlab = names(groups)[1],
                       y = level, ylab = "Relative expression",
                       title = t,
                       package = "rcartocolor", palette = "Vivid") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
        error = function(e)
          ggbetweenstats(data = x,
                         x = !!sym(groups[[1]]), xlab = names(groups)[1],
                         y = level, ylab = "Relative expression",
                         title = t,
                         package = "rcartocolor", palette = "Vivid",
                         pairwise.comparisons = FALSE) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      )
    
    # Output
    return(p)
    
  })
  
  # Wrap plots
  print(
    wrap_plots(plot) +
      plot_annotation(title = "Relative expression",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
  
   
  
}
  