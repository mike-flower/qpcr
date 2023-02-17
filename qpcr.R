#===============================================================================
# Set working directory
#===============================================================================
home <- "/Users/michaelflower/Library/Mobile Documents/com~apple~CloudDocs/Documents/bin/qpcr"
setwd(home)


#===============================================================================
# Useful commands
#===============================================================================
#rm(list = ls())
#.rs.restartR()
`%notin%` <- Negate(`%in%`)

# Extract sublist
extract_sublist <- function(x, sub) {
  y <- x[[sub]]
  return(y)
}

# Directory popup
choose_directory = function(caption = 'Select directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
  } else {
    tk_choose.dir(caption = caption)
  }
}


#===============================================================================
# Load programs
#===============================================================================
packages <- c("readxl", "ggplot2", "dplyr", "tidyr", "timevis", "tidyverse",
              "janitor", "patchwork", "pbapply", "seqinr", "shiny", "plotly",
              "ggnewscale", "ggdark", "gridExtra", "Fragman", "lubridate",
              "xlsx", "tcltk", "svDialogs", "ggstatsplot", "MASS", "relaimpo", 
              "stringr", "rstatix", "evaluate", "broom", "ggiraphExtra", 
              "xlsx", "ggpubr", "nlme", "ggpmisc", "emmeans", "data.table",
              "modelr", "plyr", "scales")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)



#===============================================================================
# Import qPCR data
#===============================================================================
datafiles <- tk_choose.files(caption = "Select qPCR data file/s",
                             multi = TRUE,
                             filter = matrix(c("All files", "*",
                                               "Text", ".txt"),
                                             2, 2, byrow = TRUE))
names(datafiles) <- basename(datafiles)

qpcr_import <- lapply(datafiles, function(x) {
  y <- read_excel(x, sheet = "Results", col_names = F) %>%
    janitor::clean_names()
  return(y)
})
names(qpcr_import)



#===============================================================================
# Import settings
#===============================================================================
settingsfiles <- tk_choose.files(caption = "Select settings file/s",
                                 multi = TRUE,
                                 filter = matrix(c("All files", "*",
                                                   "Excel", ".xlsx"),
                                                 2, 2, byrow = TRUE))
names(settingsfiles) <- basename(settingsfiles)

settings <- lapply(settingsfiles, function(x) {
  y <- read_excel(x) %>%
    mutate(row = substr(well, 1, 1),
           col = as.numeric(sub(".", "", well))) %>%
    mutate(well_short = paste0(row, col))  %>%
    dplyr::rename(well_long = well) %>%
    relocate(filename, well_long, well_short, row, col)
  return(y)
})
settings <- as.data.frame(data.table::rbindlist(settings))




#===============================================================================
# Variables
#===============================================================================
# Unique ID
unique_id = tk_select.list(choices = names(settings),
                           multiple = FALSE,
                           title = "Select unique sample identifier column")

unique_id_label = dlgInput(message = "Enter unique ID label",
                           default = unique_id)$res

# target source
target_source <- tk_select.list(choices = c("qPCR data file", "Settings file"),
                                preselect = "qPCR data file",
                                multiple = FALSE,
                                title = "Select source to take target information from")


# Calibrator samples
calibration_method = tk_select.list(choices = c("Default",
                                               "Callibrator sample/s",
                                               "Lowest GOI ct",
                                               "Highest GOI ct"),
                                   preselect = "Default",
                                   multiple = FALSE,
                                   title = paste0("Select calibrator method for comparative ct analysis \n
                                                  (Default uses mean of callibrator sample/s from settings, or if not available then the sample with lowest GOI ct)"))

calibrator_var = tk_select.list(choices = c(NA, names(settings)),
                                preselect = NA,
                                multiple = FALSE,
                                title = "Select calibrator sample column")

# Outlier threshold
outlier_threshold = as.numeric(dlgInput(message = "Enter outlier threshold [abs(value - median(value)) > threshold * sd]",
                                        default = 1.5)$res)

remove_outliers = tk_select.list(choices = c("Yes", "No"),
                                 preselect = "Yes",
                                 multiple = FALSE,
                                 title = "Remove outliers?")

# Grouping variables
xvar = tk_select.list(choices = c(NA, names(settings)),
                      preselect = NA,
                      multiple = FALSE,
                      title = "Select first grouping variable (x-axis groups)")

x_label <- dlgInput(message = "Enter x-axis label",
                    default = xvar)$res

facet_col = tk_select.list(choices = c(NA, names(settings)),
                           preselect = NA,
                           multiple = FALSE,
                           title = "Select second grouping variable (facet columns)")

facet_row = tk_select.list(choices = c(NA, names(settings)),
                           preselect = NA,
                           multiple = FALSE,
                           title = "Select third grouping variable (facet rows)")

# Exclusions
exclusion_var = tk_select.list(choices = c(NA, names(settings)),
                               preselect = NA,
                               multiple = FALSE,
                               title = "Select exclusion column (NA if not included)")
exclusion_var <- if(is.na(exclusion_var)) {"."} else {exclusion_var}

# Output directory
out.dir = choose_directory(caption = "Select output folder")



#===============================================================================
# Open a pdf
#===============================================================================
pdf(file.path(out.dir, "qpcr.pdf"),
    height = 10, width = 13, onefile = T)



#===============================================================================
# Extract information from qPCR data
#===============================================================================
# Header metrics to extract
terms <- c("Experiment File Name", "Experiment Name", "Experiment Run End Time",
           "Experiment Type", "Instrument Name", "Instrument Serial Number",
           "Instrument Type", "Passive Reference", "Quantification Cycle Method")

# Well short vector
well_short_vars <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", 
                     "A11", "A12", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", 
                     "B9", "B10", "B11", "B12", "C1", "C2", "C3", "C4", "C5", "C6", 
                     "C7", "C8", "C9", "C10", "C11", "C12", "D1", "D2", "D3", "D4", 
                     "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "E1", "E2", 
                     "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", 
                     "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
                     "F11", "F12", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", 
                     "G9", "G10", "G11", "G12", "H1", "H2", "H3", "H4", "H5", "H6", 
                     "H7", "H8", "H9", "H10", "H11", "H12")

# Extract header metrics and ct data
qpcr_data <- lapply(setNames(names(qpcr_import), names(qpcr_import)), function(n) {
  
  # Extract the plate
  q = qpcr_import[[n]]
  
  # Extract header terms
  headers <- lapply(terms, function(t) {
    term_rename <- t %>% janitor::make_clean_names()
    term_location <- which(q == t, arr.ind = TRUE)
    value <- as.character(q[term_location[1], term_location[2]+1])
    out <- data.frame(term = term_rename,
                      value = value)
    return(out)
  })
  headers <- as.data.frame(data.table::rbindlist(headers))
  
  # Transpose header terms to make them rbind friendly 
  headers <- setNames(data.frame(t(headers[,-1])), headers[,1]) %>%
    mutate(filename = n) %>%
    relocate(filename)
  
  # Extract Ct data
  data_start = which(q[,1] == "Well")
  ct <- q[data_start:nrow(q),] %>%
    row_to_names(row_number = 1) %>%
    janitor::clean_names() %>%
    mutate(filename = n) %>%
    relocate(filename) %>%
    dplyr::mutate(across(c("ct"), ~as.numeric(.))) %>%
    dplyr::rename(well_n = well,
                  well_short = well_position) %>%
    filter(well_short %in% well_short_vars)
  
  # Output
  return(list(headers = headers,
              ct = ct))
  
})

# Consolidate data from different plates
headers <- rbind.fill(lapply(qpcr_data, extract_sublist, sub = "headers"))
ct <- rbind.fill(lapply(qpcr_data, extract_sublist, sub = "ct"))
  

# Export header information
unlink(file.path(out.dir, "qpcr.xlsx"))
write.xlsx(headers,
           file = file.path(out.dir, "qpcr.xlsx"),
           sheetName = "headers", append = T, row.names = F)





#===============================================================================
# Annotate ct data with target information
#===============================================================================
if (target_source == "qPCR data file") {
  
  ct <- ct %>%
    select(filename, well_short, reporter, target_name, ct) %>%
    dplyr::rename(target = target_name)
  
} else if (target_source == "Settings file") {
  
  target_map <- settings %>%
    select(filename, well_short, starts_with(c("reporter", "target"))) %>%
    pivot_longer(cols = -c(filename, well_short),
                 names_pattern = "(.*)\\d$",
                 names_to = ".value") %>%
    filter(complete.cases(.)) # https://stackoverflow.com/questions/75465900/using-pivot-longer-on-two-sets-of-columns/75465941#75465941
  
  ct <- ct %>%
    select(filename, well_short, reporter, ct) %>%
    left_join(target_map,
              by = c("filename", "well_short", "reporter")) %>%
    relocate(target, .after = "reporter")
  
}


#===============================================================================
# Select GOI and HK genes
#===============================================================================
goi <- tk_select.list(choices = sort(unique(ct$target)),
                      multiple = TRUE,
                      title = "Select target/s of interest")

hk <- tk_select.list(choices = sort(unique(ct$target)),
                     preselect = setdiff(sort(unique(ct$target)), goi),
                      multiple = TRUE,
                      title = "Select housekeeping target/s")

allow_failed_hk = tk_select.list(choices = c("Yes", "No"),
                                   preselect = "No",
                                   multiple = FALSE,
                                   title = "Permit geomean calculation when one or more housekeeper fails")



#===============================================================================
# Annotate ct data with settings information
#===============================================================================
# Target levels
levels_target <- sort(unique(ct[["target"]]))

# Unique ID levels
levels_id <- unique(settings[[unique_id]])

# Grouping variable levels
levels_xvar <- unique(settings[[xvar]])
levels_facetcol <- unique(settings[[facet_col]])
levels_facetrow <- unique(settings[[facet_row]])

# Annotate with settings data and set variable levels
ct <- ct %>%
  left_join(settings %>%
              select(-c(well_long, row, col),
                     -starts_with(c("reporter", "target"))),
            by = c("filename", "well_short")) %>%
  relocate(c(reporter, target, ct), .after = last_col()) %>%
  dplyr::mutate(target = factor(target, levels = levels_target),
         !!sym(unique_id) := factor(!!sym(unique_id), levels = levels_id),
         !!sym(xvar) := factor(!!sym(xvar), levels = levels_xvar)) %>%
  droplevels()

if (!is.na(facet_col)) {
  ct <- ct %>%
    dplyr::mutate(!!sym(facet_col) := factor(!!sym(facet_col), levels = levels_facetcol))
}

if (!is.na(facet_row)) {
  ct <- ct %>%
    dplyr::mutate(!!sym(facet_row) := factor(!!sym(facet_row), levels = levels_facetrow))
}




#===============================================================================
# Split ct data into a list by target
#===============================================================================
ct <- split(ct, f = ct$target)




#===============================================================================
# Find outliers
#===============================================================================
# Find outliers for each target separately
ct <- lapply(ct, function(x) {
  
  # Target
  t = x %>%
    slice(1) %>%
    pull(target)
  
  # Find outliers
  y <- x %>%
    group_by(!!sym(unique_id)) %>%
    dplyr::mutate(tech_rep = row_number()) %>%
    dplyr::mutate(outlier = abs(ct - median(ct, na.rm = T)) > outlier_threshold * sd(ct, na.rm = T)) %>%
    mutate(outlier_threshold = outlier_threshold) %>%
    relocate(c(tech_rep, outlier_threshold, outlier), .after = ct) %>%
    mutate(outlier = ifelse(is.na(outlier), FALSE, outlier)) %>%
    mutate(outlier = factor(outlier, levels = c(TRUE, FALSE))) %>%
    ungroup()
  
  # Plot
  plot <-
    ggplot(y, aes(x = !!sym(unique_id), y = ct)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width=0.1) +
    geom_jitter(aes(colour = outlier), width = 0.1,height = 0, size = 3) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs(title = t,
         x = unique_id_label,
         y = "Threshold cycle (ct)",
         colour = "Outlier") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          # text = element_text(size = 15)
          )
  
  # Output
  return(list(data = y,
              plot_outliers = plot))
})

print(
  wrap_plots(lapply(ct, extract_sublist, sub = "plot_outliers")) +
    plot_annotation(title = "Threshold cycle",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
)



#===============================================================================
# Average tech reps for each sample (removing outliers and exclusions)
#===============================================================================
ct <- lapply(ct, function(x) {
  
  # Count outliers and exclusions
  oe <- x[["data"]] %>%
    group_by(!!sym(unique_id)) %>%
    dplyr::summarise(n_total = n(),
                     n_outliers = sum(outlier == TRUE),
                     n_excluded = sum(!is.na(!!sym(exclusion_var))))
  
  # Remove outliers and exclusions
  y <- x[["data"]] %>%
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
    mutate(target = x[["data"]] %>%
             slice(1) %>%
             pull(target)) %>%
    relocate(target)
  
  # Summarise number of tech reps included (has to be done separately)
  n_included <- y %>%
    group_by(!!sym(unique_id)) %>%
    dplyr::summarise(n_included = n())
  
  # Merge mean, sd and n
  z <- z %>%
    left_join(oe, by = unique_id) %>%
    left_join(n_included, by = unique_id) %>%
    relocate(c(ct, sd), .after = "n_included")
  
  # Output
  out <- append(x, list(summary = z))
  return(out)
  
})

# Export summaries
summaries <- rbindlist(lapply(ct, extract_sublist, sub = "summary"))
wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
removeSheet(wb, sheetName = "summary")
saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
write.xlsx(summaries, file = file.path(out.dir, "qpcr.xlsx"),
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
    rbindlist(lapply(ct[names(ct) %in% goi], extract_sublist, sub = "summary")) %>%
    slice(which.min(ct_mean)) %>%
    slice(1) %>%
    pull(!!sym(unique_id)))
    
} else if (calibration_method == "Callibrator sample/s") {
  
  calibrators <- as.character(
    settings %>%
    filter(!is.na(!!sym(calibrator_var))) %>%
    distinct(!!sym(unique_id)) %>%
    pull(!!sym(unique_id)))
  
} else if (calibration_method == "Lowest GOI ct") {
  
  calibrators <- as.character(
    rbindlist(lapply(ct[names(ct) %in% goi], extract_sublist, sub = "summary")) %>%
    slice(which.min(ct_mean)) %>%
    slice(1) %>%
    pull(!!sym(unique_id)))
  
} else if (calibration_method == "Highest GOI ct") {
  
  calibrators <- as.character(
    rbindlist(lapply(ct[names(ct) %in% goi], extract_sublist, sub = "summary")) %>%
    slice(which.max(ct_mean)) %>%
    slice(1) %>%
    pull(!!sym(unique_id)))
  
}

calibrators


#===============================================================================
# Calculate control mean ct for each target
#===============================================================================
control_mean_ct <- lapply(ct, function(x) {

  return(
    x[["summary"]] %>%
    filter(!!sym(unique_id) %in% calibrators) %>%
    summarise(ct = mean(ct, na.rm = TRUE)) %>%
    pull(ct)
  )
  
})

control_mean_ct




#===============================================================================
# Comparative ct calculations
#===============================================================================
# Calculate 2^dct (relative quantity, RQ) for each target
comparative_ct <- lapply(ct, function(x) {
  
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

# rbind targets, pivot wider and calculate geomean
comparative_ct <- rbindlist(comparative_ct) %>%
  pivot_wider(names_from = target,
              values_from = c(ct, cmct, dct, twodct)) %>%
  dplyr::mutate(calibration_method = calibration_method,
                calibrator = !!sym(unique_id) %in% calibrators) %>%
  relocate(all_of(c("calibration_method", "calibrator")), .after = unique_id) %>%
  rowwise() %>%
  dplyr::mutate(geomean_hk = case_when(allow_failed_hk == "Yes" ~ 
                                         exp(mean(log(c_across(c(paste0("twodct_", hk)))), na.rm = TRUE)),
                                       TRUE ~
                                         exp(mean(log(c_across(c(paste0("twodct_", hk)))))))) %>%
  ungroup()
  

# For each GOI calculate relative expression
for (g in goi) {
  
  comparative_ct <- comparative_ct %>%
    dplyr::mutate(!!sym(paste0("relative_expression_", g)) :=
                    !!sym(paste0("twodct_", g)) / geomean_hk)
  
}


#===============================================================================
# Add annotation from settings
#===============================================================================
# Select annotation from settings
s <- settings %>%
  select(-c(filename, well_long, well_short, row, col, !!sym(exclusion_var), !!sym(calibrator_var)),
         -starts_with(c("reporter", "target"))) %>%
  distinct() %>%
  filter(!is.na(!!sym(unique_id)))

# Add settings annotation and set levels
comparative_ct <- comparative_ct %>%
  left_join(s, by = unique_id) %>%
  relocate(all_of(names(s)))



#===============================================================================
# Add levels
#===============================================================================
# Unique ID and first grouping variable
comparative_ct <- comparative_ct %>%
  dplyr::mutate(!!sym(unique_id) := factor(!!sym(unique_id), levels = levels_id),
                !!sym(xvar) := factor(!!sym(xvar), levels = levels_xvar)) %>%
  droplevels()

if (!is.na(facet_col)) {
  comparative_ct <- comparative_ct %>%
    dplyr::mutate(!!sym(facet_col) := factor(!!sym(facet_col), levels = levels_facetcol))
}

if (!is.na(facet_row)) {
  comparative_ct <- comparative_ct %>%
    dplyr::mutate(!!sym(facet_row) := factor(!!sym(facet_row), levels = levels_facetrow))
}


#===============================================================================
# Export comparative ct table
#===============================================================================
wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
removeSheet(wb, sheetName = "comparative_ct")
saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
write.xlsx(as.data.frame(comparative_ct), file = file.path(out.dir, "qpcr.xlsx"),
           sheetName = "comparative_ct", append = T, row.names = F)


                  
#===============================================================================
# Plot relative expression
#===============================================================================
plots <- lapply(setNames(goi, goi), function(g) {
  
  # t-test
  ttest <- comparative_ct %>%
    t_test(as.formula(paste0("relative_expression_", g, " ~ ", xvar))) %>%
    add_y_position()
  ttest
  
  # Plot ttest
  plot_ttest <-
    ggplot(comparative_ct, aes(x = !!sym(xvar), y = !!sym(paste0("relative_expression_", g)))) +
    geom_boxplot(outlier.size = 0, position = position_dodge(width = 0.9)) +
    geom_jitter(width = 0.1,height = 0, size = 3) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 4, show.legend = FALSE) + 
    stat_summary(fun = mean, colour = "darkred", geom = "text", show.legend = FALSE, 
                 vjust = -0.7,
                 aes(label = after_stat(round(y, 2)))) +
    stat_pvalue_manual(ttest, 
                       label = ifelse("p.adj.signif" %in% names(ttest), "p.adj.signif", "p"),
                       tip.length = 0.01, hide.ns = T) +
    scale_y_continuous(expand = expansion(mult = c(NA, 0.2)),
                       limits = c(0, NA)) +
    labs(title = g,
         x = x_label,
         y = "Relative expression") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          text = element_text(size = 20))
  
  
  # ggstatsplot
  plot_ggstatsplot <-
    tryCatch(
      
      ggbetweenstats(data = comparative_ct,
                     x = !!sym(xvar), xlab = x_label,
                     y = !!sym(paste0("relative_expression_", g)), ylab = "Relative expression",
                     title = g) +
                     # package = "rcartocolor", palette = "Vivid") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
              text = element_text(size = 15))
      
      , error = function(e)
        
        ggbetweenstats(data = comparative_ct,
                       x = !!sym(xvar), xlab = x_label,
                       y = !!sym(paste0("relative_expression_", g)), ylab = "Relative expression",
                       title = g,
                       pairwise.comparisons = FALSE) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
              text = element_text(size = 15))
      
        )
        
  # Output
  return(list(plot_ttest = plot_ttest,
              plot_ggstatsplot = plot_ggstatsplot))
    
})

print(
  wrap_plots(lapply(plots, extract_sublist, sub = "plot_ttest"))
)

print(
  wrap_plots(lapply(plots, extract_sublist, sub = "plot_ggstatsplot"))
)




#===============================================================================
# Close pdf
#===============================================================================
dev.off()


#===============================================================================
# Save data
#===============================================================================
# Save data
save(list = ls(),
     file = file.path(out.dir, "qpcr.RData"))

# Load data
# out.dir = "/Users/michaelflower/Library/Mobile Documents/com~apple~CloudDocs/Documents/bin/qpcr/demo_results/Freja"
# lnames = load(file = file.path(out.dir, "qpcr.RData"))
# lnames