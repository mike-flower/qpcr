#===============================================================================
# qPCR stats
#===============================================================================
qpcr_stats <- function(quantification, grouping_vars, dependent_vars, y_log, hide_ns) {
  
  # # Manual
  # quantification = quantification
  # grouping_vars = grouping_vector
  # dependent_vars = dependent_vector
  # y_log = y_log
  # hide_ns = hide_ns

  
  #===============================================================================
  # Select settings data to annotate
  #===============================================================================
  s <- settings %>%
    select(any_of(c(unique_id, unname(grouping_vars)))) %>%
    distinct()
  
  
  #===============================================================================
  # Remove NA, annotate with settings data and set factor levels
  #===============================================================================
  quantification_clean <- lapply(quantification, function(x) {
    
    # Remove NA
    y <- x %>%
      drop_na(any_of(c(unique_id, unname(dependent_vars))))
    
    # Add settings
    y <- y %>%
      left_join(s, by = unique_id) %>%
      relocate(names(s))
    
    # Iterate through grouping variables removing NA and adding factor levels
    for (g in grouping_vars) {
      
      y <- y %>%
        filter(!is.na(!!sym(g))) %>%
        dplyr::mutate(!!sym(g) := factor(!!sym(g), levels = unique(settings[[g]]))) %>%
        droplevels()
      
    }
    
    # Output
    return(y)
    
  })
  
  
  
  #===============================================================================
  # Stats on each dependent variable
  #===============================================================================
  # Perform ttest for each dependent variable
  ttest_plot <- lapply(setNames(names(dependent_vars), names(dependent_vars)), function(d) {
    
    # # Manual
    # d = names(dependent_vars)[1]
    
    #===============================================================================
    # t-test
    #===============================================================================
    # t-test each target
    ttest_plot <- lapply(quantification_clean, function(x) {
      
      # # Manual
      # x = quantification_clean[[1]]
      
      # Extract target
      t <- as.character(
        x %>%
          slice(1) %>%
          pull(target)
      )
      t
      
      # Identify groups with >1 sample (to keep for t-test)
      group_count <- x %>%
        group_by_at(unname(grouping_vars)) %>%
        dplyr::summarise(n_samples = n())
      
      # Do a ttest if there are enough groups with >1 sample
      # t-test
      if(nrow(group_count %>% filter(n_samples > 1)) > 1) {
        
        ttest <- x %>%
          left_join(group_count, by = unname(grouping_vars)) %>%
          filter(n_samples > 1) %>%
          droplevels() %>% # https://stackoverflow.com/questions/75227674/r-t-test-errors-with-not-enough-y-observations-but-there-are-plenty-of-x-an
          group_by_at(unname(grouping_vars)[-1]) %>% # First grouping variable is used as independent variable in t-test
          t_test(as.formula(paste0(dependent_vars[[d]], " ~ ", unname(grouping_vars)[1]))) %>%
          add_y_position()
        
      } else {
        
        ttest <- NA
        
      }
      
      ttest
      
      
      # Export t-test
      if (is.data.frame(ttest)) {
        rm_excel_sheet(file = file.path(out_dir, "qpcr.xlsx"),
                       sheetname = paste0("ttest_", dependent_vars[[d]], "_", t))
        write.xlsx(data.frame(ttest) %>%
                     mutate(groups = as.character(groups)),
                   file = file.path(out_dir, "qpcr.xlsx"),
                   sheetName = paste0("ttest_", dependent_vars[[d]], "_", t), append = T, row.names = F)
      }
      
      
      # Facet formula
      if (length(grouping_vars) > 1) {
        facet_formula <- as.formula(paste(unname(grouping_vars)[-1], collapse = " ~ "))
      } else {
        facet_formula <- NA
      }
      facet_formula
      
      
      # Plot
      p <-
        ggplot(x, aes(x=!!sym(grouping_vars[[1]]), y=!!sym(dependent_vars[[d]]))) +
          { if (is.formula(facet_formula)) 
            facet_grid(facet_formula, scales = "free_x", space = "free_x")
          } + # https://stackoverflow.com/questions/22915337/if-else-condition-in-ggplot-to-add-an-extra-layer
          geom_violin(trim = TRUE) +
          geom_boxplot(outlier.shape = NA, width = 0.1) +
          geom_jitter(width = 0.1, height = 0, size = 3) +
          { if (is.data.frame(ttest)) 
            stat_pvalue_manual(ttest,
                               label = ifelse("p.adj.signif" %in% names(ttest), "p.adj.signif", "p"),
                               tip.length = 0.01, hide.ns = hide_ns)
          } +
          stat_summary(fun = mean, colour = "darkred", geom = "point", 
                       shape = 18, size = 4, show.legend = FALSE) + 
          stat_summary(fun = mean, colour = "darkred", geom = "text", show.legend = FALSE, 
                       vjust = -0.7,
                       aes(label = after_stat(ifelse(y > 0.01,
                                                     round(y, 2),
                                                     formatC(y, format = "e", digits = 1))))) +
          # geom_hline(yintercept = 1, linetype = "dashed") +
          { if (y_log) 
            scale_y_continuous(trans = "log",
                               labels = function(x) { 
                                 ifelse(x > 1000,
                                        formatC(x, format = "e", digits = 1),
                                        round(x, digits = 2)) })
          } +
          labs(title = t,
               x = names(grouping_vars)[1],
               y = d) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      
      
      # Output
      # return(p)
      return(list(ttest = ttest,
                  plot = p))
      
    })
    
    # Output
    return(ttest_plot)
    
  })
      
  
  
  
  
  
  # Extract ttests
  ttests <- lapply(ttest_plot, function(x) {
    y <- lapply(x, extract_sublist, sub = "ttest")
    return(y)
  })
  
  # Extract plots
  plots <- lapply(ttest_plot, function(x) {
    y <- lapply(x, extract_sublist, sub = "plot")
    return(y)
  })
  
  # Print the plots
  for (d in names(plots)) {
    
    # Print to pdf
    print(
      wrap_plots(plots[[d]][!is.na(plots[[d]])]) +
        plot_annotation(title = d,
                        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    )
    
    # Save as png
    ggsave(file.path(out_dir, paste0("plot_", dependent_vars[[d]], ".png")))
    
  }
  
  
  
  # Output ttest tables
  return(ttests)
  
}
  
  
  
  
