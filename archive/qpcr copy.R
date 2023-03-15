#===============================================================================
# Set working directory
#===============================================================================
setwd("/Users/michaelflower/Library/Mobile Documents/com~apple~CloudDocs/Documents/bin/qpcr")



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
              "modelr", "plyr", "scales", "ggrepel")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)
rm(packages)



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
rm(datafiles)






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
rm(settingsfiles)






#===============================================================================
# Variables common to all analyses
#===============================================================================
# Analysis method
analysis_method = tk_select.list(choices = c("Comparative ct (multiple housekeepers)",
                                             "Comparative ct (2^-ddct)",
                                             "Standard curve"),
                                 preselect = "Comparative ct (multiple housekeepers)",
                                 multiple = FALSE,
                                 title = paste0("Select qPCR analysis method\n
                                       The 'Multiple housekeepers' method uses the formula '2^dctGOI / geomean(2^dctREFS)'.\n
                                       The '2^-ddct' formula has been bodged to handle multiple hk targets too."))


# Unique ID
unique_id = tk_select.list(choices = names(settings),
                           multiple = FALSE,
                           title = "Select unique sample identifier column")

# Unique ID label
unique_id_label = dlgInput(message = "Enter unique ID label",
                           default = unique_id)$res


# Target source
target_source <- tk_select.list(choices = c("qPCR data file", "Settings file"),
                                preselect = "qPCR data file",
                                multiple = FALSE,
                                title = "Select source to take target (gene) information from")


# NC samples (used to exclude HK genes if they fail in a real sample)
empty_samples = tk_select.list(choices = c(NA,
                                           unique(settings[[unique_id]])[!is.na(unique(settings[[unique_id]]))]),
                               preselect = NA,
                               multiple = TRUE,
                               title = "Select empty sample/s\n
                              (samples you expect to be empty.\n
                              Used to exclude housekeeping genes if amplification fails in a real sample)")


# Outlier threshold
outlier_threshold = as.numeric(dlgInput(message = "Enter outlier threshold\n
                                        [abs(value - median(value)) > threshold * sd]",
                                        default = 1.5)$res)

# Remove outliers
remove_outliers = tk_select.list(choices = c("Yes", "No"),
                                 preselect = "Yes",
                                 multiple = FALSE,
                                 title = "Remove outliers from analysis?")


# Exclusions
exclusion_var = tk_select.list(choices = c(NA, names(settings)),
                               preselect = NA,
                               multiple = FALSE,
                               title = "Select column of excluded wells (NA if no exclusion column)")
exclusion_var <- if(is.na(exclusion_var)) {"."} else {exclusion_var}

# Output directory
out.dir = choose_directory(caption = "Select output folder")



#===============================================================================
# Variables for comparative ct analyses
#===============================================================================
if (grepl("Comparative", analysis_method, fixed = TRUE)) {
  
  # Calibrator samples
  calibration_method = tk_select.list(choices = c("Default",
                                                  "Callibrator sample/s",
                                                  "Lowest GOI ct",
                                                  "Highest GOI ct",
                                                  "Select sample/s"),
                                      preselect = "Default",
                                      multiple = FALSE,
                                      title = paste0("Select calibrator method for comparative ct analysis \n
                                                  (Default uses mean Ct of callibrator sample/s from Settings file, or if not available then the Ct of the sample with lowest gene of interest Ct)"))
  
  if(calibration_method == "Default" |
     calibration_method == "Callibrator sample/s") {
    
    calibrator_var = tk_select.list(choices = c(NA, names(settings)),
                                    preselect = NA,
                                    multiple = FALSE,
                                    title = "Select calibrator sample column")
    
  } else { calibrator_var = "" }
  
  
  if(calibration_method == "Select sample/s") {
    
    calibrator_samples = tk_select.list(choices = unique(settings[[unique_id]]),
                                        multiple = TRUE,
                                        title = "Select calibrator sample/s")
    
  } else { calibrator_samples = NA }
  
  
  # Exclude a housekeeping targets if they don't work for all samples
  exclude_failed_hk = tk_select.list(choices = c("Yes", "No"),
                                     preselect = "Yes",
                                     multiple = FALSE,
                                     title = "Exclude housekeeping targets (genes) if they don't work for all samples")
  
}







#===============================================================================
# Standard curve variables
#===============================================================================
if (analysis_method == "Standard curve") {
  
  # Standards
  standard_var = tk_select.list(choices = names(settings),
                                multiple = FALSE,
                                title = "Select standard concentration column")
  
  # Standard concentrations
  conc_unit <- dlgInput(message = "Enter concentration unit",
                        default = "ng/µL")$res
  
  # standard_label = dlgInput(message = "Enter standard concentration label",
  #                           default = standard_var)$res
  
  # Source volume
  vol_var = tk_select.list(choices = names(settings),
                           multiple = FALSE,
                           title = "Select sample volume column")
  
  # Is this a quantifiler experiment?
  quantifiler = tk_select.list(choices = c("No", "Yes"),
                               preselect = "No",
                               multiple = FALSE,
                               title = "Quantifiler experiment?")
  
}






#===============================================================================
# Open a pdf
#===============================================================================
pdf(file.path(out.dir, "qpcr.pdf"),
    height = 10, width = 13, onefile = T)







#===============================================================================
# Extract information from qPCR data
#===============================================================================
# Header metrics to extract
terms_local <- c("Experiment File Name", "Experiment Name", "Experiment Run End Time",
                 "Experiment Type", "Instrument Name", "Instrument Serial Number",
                 "Instrument Type", "Passive Reference", "Quantification Cycle Method")

terms_cloud <- c("File Name", "Operator", "Barcode", "Instrument Type", "Block Type", 
                 "Instrument Name", "Instrument Serial Number", "Heated Cover Serial Number",
                 "Block Serial Number","Run Start Date/Time", "Run End Date/Time",
                 "Run Duration", "Sample Volume", "Cover Temperature", "Passive Reference",
                 "PCR Stage/Step Number", "Quantification Cycle Method", "Analysis Date/Time",
                 "Software Name and Version", "Plugin Name and Version", "Exported By",
                 "Exported On")

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
  
  # Select term set
  if(length(which(q[,1] == terms_local[2])) > 0) {
    terms <- terms_local
  } else {
    terms <- terms_cloud
  }
  terms
  
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
  
  # Extract ct data
  data_start = which(q[,1] == "Well")
  ct <- q[data_start:nrow(q),] %>%
    row_to_names(row_number = 1) %>%
    janitor::clean_names() 
  
  # If data is from the cloud, convert cq to ct
  names(ct) <- gsub("cq", "ct", names(ct))
  
  # Format ct table
  ct <- ct %>%
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


# Export raw ct data
wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
removeSheet(wb, sheetName = "ct")
saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
rm(wb)
write.xlsx(as.data.frame(ct),
           file = file.path(out.dir, "qpcr.xlsx"),
           sheetName = "ct", append = T, row.names = F)





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

# Make target names machine friendly (for later t-testing)
ct$target <- gsub(" ", "_", ct$target)







#===============================================================================
# Set target purposes
#===============================================================================
if (grepl("Comparative", analysis_method, fixed = TRUE)) {
  
  # Select gene/s of interest
  goi <- tk_select.list(choices = unique(ct$target)[!is.na(unique(ct$target))],
                        multiple = TRUE,
                        title = "Select target/s of interest")
  goi <- goi[order(match(goi, unique(ct$target)))] 
  
  # Select housekeeping genes
  hk <- tk_select.list(choices = unique(ct$target)[!is.na(unique(ct$target))],
                       preselect = setdiff(sort(unique(ct$target)), goi),
                       multiple = TRUE,
                       title = "Select housekeeping target/s")
  hk <- hk[order(match(hk, unique(ct$target)))]
  
}


# IPC target
ipc <- tk_select.list(choices = c(NA, unique(ct$target)[!is.na(unique(ct$target))]),
                      preselect = "NA",
                      multiple = FALSE,
                      title = "Select internal positive control (IPC) target\n
                      (included in PCR mix for accuracy of pipetting)")


# Standard curve target variables
if (analysis_method == "Standard curve" &&
    quantifiler == "Yes") {
  
  quality_numerator <- tk_select.list(choices = c(NA, unique(ct$target)[!is.na(unique(ct$target))]),
                                      preselect = "NA",
                                      multiple = FALSE,
                                      title = "Select quality score numerator target (e.g. Small)")
  
  quality_denominator <- tk_select.list(choices = c(NA, unique(ct$target)[!is.na(unique(ct$target))]),
                                        preselect = "NA",
                                        multiple = FALSE,
                                        title = "Select quality score denominator target (e.g. Large)")
  
}






#===============================================================================
# Annotate ct data with settings information
#===============================================================================
# Target levels (in the order they were run on the plate)
levels_target <- (unique(ct[["target"]]))

# Unique ID levels (in the order they appear in the settings file)
levels_id <- unique(settings[[unique_id]])

# Annotate with settings data and set variable levels
ct <- ct %>%
  left_join(settings %>%
              select(-c(well_long, row, col),
                     -starts_with(c("reporter", "target"))),
            by = c("filename", "well_short")) %>%
  relocate(c(reporter, target, ct), .after = last_col()) %>%
  dplyr::mutate(target = factor(target, levels = levels_target),
                !!sym(unique_id) := factor(!!sym(unique_id), levels = levels_id)) %>%
  droplevels()




#===============================================================================
# Plot ct for targets
#===============================================================================
print(
  ggplot(ct, aes(x = target, y = ct)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1) +
    geom_jitter(aes(colour = !!sym(unique_id)), width = 0.1,height = 0, size = 3) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 4, show.legend = FALSE) + 
    stat_summary(fun = mean, colour = "darkred", geom = "text", show.legend = FALSE, 
                 vjust = -0.7,
                 aes(label = after_stat(round(y, 2)))) +
    labs(title = "Targets",
         x = "Target",
         y = "Threshold cycle (ct)",
         colour = unique_id_label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          # text = element_text(size = 15)
    )
)


#===============================================================================
# Plot ct for samples
#===============================================================================
print(
  ggplot(ct, aes(x = !!sym(unique_id), y = ct)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1) +
    geom_jitter(aes(colour = target), width = 0.1,height = 0, size = 3) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 4, show.legend = FALSE) + 
    stat_summary(fun = mean, colour = "darkred", geom = "text", show.legend = FALSE, 
                 vjust = -0.7,
                 aes(label = after_stat(round(y, 2)))) +
    labs(title = "Samples",
         x = unique_id_label,
         y = "Threshold cycle (ct)",
         colour = "Target") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          # text = element_text(size = 15)
    )
)





#===============================================================================
# Split ct data into a list by target
#===============================================================================
ct <- split(ct, f = ct$target)




#===============================================================================
# Find outlier internal positive control (IPC) wells (IF SPECIFIED)
#===============================================================================
# Synthetic IPC template is included in the primer mix
# Measures accuracy of pipetting master mix

if (!is.na(ipc)) {

  # Find IPC outliers
  ipc_outliers <- ct[[ipc]] %>%
    dplyr::mutate(ipc_outlier = abs(ct - median(ct, na.rm = T)) > outlier_threshold * sd(ct, na.rm = T))
  
  # Plot ipc variability
  print(
    ggplot(ipc_outliers, aes(x = target, y = ct)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(colour = ipc_outlier), width=0.1, height=0, size = 3) +
      scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      geom_text_repel(data = ipc_outliers %>% filter(ipc_outlier == T),
                      aes(label = well_short), colour = "red", size = 5) +
      labs(title = "IPC outlier wells",
           x = "Target",
           y = "Threshold cycle (ct)",
           colour = "IPC outlier") +
      theme(text = element_text(size = 20))
  )
  
  
  # Export ipc outlier table
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "ipc_outlier")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(as.data.frame(ipc_outliers),
             file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "ipc_outliers", append = T, row.names = F)

}




#===============================================================================
# Find outlier technical replicates
#===============================================================================
# Find outliers for each target separately
ct <- lapply(ct, function(x) {
  
  # Target
  t = as.character(
    x %>%
      slice(1) %>%
      pull(target)
  )
  
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
  
  # Output
  return(y)
  
})


# Export ct data
ct_export <- rbindlist(ct)
wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
removeSheet(wb, sheetName = "outliers")
saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
rm(wb)
write.xlsx(ct_export,
           file = file.path(out.dir, "qpcr.xlsx"),
           sheetName = "outliers", append = T, row.names = F)


# Plot ct data with outliers
plot <- lapply(ct, function(x) {
  
  # Target
  t = as.character(
    x %>%
      slice(1) %>%
      pull(target)
  )
  
  # Plot
  plot <-
    ggplot(x, aes(x = !!sym(unique_id), y = ct)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(outlier.shape = NA, width = 0.1) +
    geom_jitter(aes(fill = outlier), width = 0.1,height = 0, shape = 21, size = 3, alpha = 0.5) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs(title = t,
         x = unique_id_label,
         y = "Threshold cycle (ct)",
         fill = "Outlier") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          # text = element_text(size = 15)
    )
  
  # Output
  return(plot)
  
})

# Wrap plots
print(
  wrap_plots(plot) +
    plot_annotation(title = "Outliers",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
)






#===============================================================================
# Comparative ct (multiple housekeepers)
#===============================================================================
if (analysis_method == "Comparative ct (multiple housekeepers)") {
  
  # Prepare targets for comparative ct analysis
  source("./functions/comparative_prep.R")
  comparative_vars <- comparative_prep(ct = ct,
                                       settings = settings)
  
  # Assign variables
  for (n in names(comparative_vars)) {
    assign(n, comparative_vars[[n]])
  }
  
  # Comparative ct analysis
  # (N.B. use summary, as each unique sample has 1 value for goi and all hk)
  source("./functions/multiple_hk.R")
  comparative_ct <- multiple_hk(summary = summary)
  
  # Export relative expression data
  wb = loadWorkbook(file.path(out.dir, "qpcr.xlsx"))
  removeSheet(wb, sheetName = "relative_expression")
  saveWorkbook(wb, file.path(out.dir, "qpcr.xlsx"))
  rm(wb)
  write.xlsx(as.data.frame(comparative_ct),
             file = file.path(out.dir, "qpcr.xlsx"),
             sheetName = "relative_expression", append = T, row.names = F)
  
  # Plot relative expression
  plots <- lapply(setNames(goi, goi), function(g) {
    
    # Bar plot
    plot <-
      ggplot(comparative_ct, aes(x = !!sym(unique_id), y = !!sym(g))) +
      geom_bar(aes(fill = calibrator), stat = "identity", colour = "black") +
      scale_fill_manual(values = c("TRUE" = "darkgrey", "FALSE" = "blue")) +
      geom_text(aes(label = round(!!sym(g), 2)), vjust=-0.3) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_y_continuous(expand = expansion(mult = c(NA, 0.1)),
                         limits = c(0, NA)) +
      labs(title = g,
           x = unique_id_label,
           y = paste0(g, " relative expression"),
           fill = "Calibrator") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            # text = element_text(size = 15)
      )
    
    # Output
    return(plot)
    
  })
  
  # Wrap the plots
  print(
    wrap_plots(plots) +
      plot_annotation(title = "Relative expression",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
  
  
}



#===============================================================================
# Comparative ct (2^-ddct)
#===============================================================================
if (analysis_method == "Comparative ct (2^-ddct)") {
  
  # Prepare targets for comparative ct analysis
  source("./functions/comparative_prep.R")
  comparative_vars <- comparative_prep(ct = ct,
                                       settings = settings)
  
  # Assign variables
  for (n in names(comparative_vars)) {
    assign(n, comparative_vars[[n]])
  }
  
  # Comparative ct analysis
  # (N.B. use summary, as each unique sample has 1 value for goi and all hk)
  source("./functions/twoddct.R")
  comparative_ct <- twoddct(summary = summary)
  
  # Plot relative expression
  plots <- lapply(setNames(goi, goi), function(g) {
    
    # Bar plot
    plot <-
      ggplot(comparative_ct, aes(x = !!sym(unique_id), y = !!sym(g))) +
      geom_bar(aes(fill = calibrator), stat = "identity", colour = "black") +
      scale_fill_manual(values = c("TRUE" = "darkgrey", "FALSE" = "blue")) +
      geom_text(aes(label = round(!!sym(g), 2)), vjust=-0.3) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_y_continuous(expand = expansion(mult = c(NA, 0.1)),
                         limits = c(0, NA)) +
      labs(title = g,
           x = unique_id_label,
           y = paste0(g, " relative expression"),
           fill = "Calibrator") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            # text = element_text(size = 15)
      )
    
    # Output
    return(plot)
    
  })
  
  # Wrap the plots
  print(
    wrap_plots(plots) +
      plot_annotation(title = "Relative expression",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
  
}



#===============================================================================
# Standard curve
#===============================================================================
if (analysis_method == "Standard curve") {
  
  source("./functions/standard_curve.R")
  sc <- standard_curve(ct = ct)

}






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

# # Load data
# lnames = load(file = file.path(out.dir, "qpcr.RData"))
# lnames



