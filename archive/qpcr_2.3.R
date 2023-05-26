#===============================================================================
# Set working directory
#===============================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio
getwd()


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
choose_directory = function(caption = "Select directory",
                            default = "") {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption, default = default) 
  } else {
    tk_choose.dir(caption = caption, default = default)
  }
}

# Remove excel sheet
rm_excel_sheet <- function(file, sheetname) {
  wb = loadWorkbook(file)
  removeSheet(wb, sheetName = sheetname)
  saveWorkbook(wb, file)
  rm(wb)
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
  y <- read_excel(x) # %>%
    # mutate(row = substr(well, 1, 1),
           # col = as.numeric(sub(".", "", well))) %>%
    # mutate(well_short = paste0(row, col))  %>%
    # dplyr::rename(well_long = well) %>%
    # relocate(filename, well_long, well_short, row, col)
  return(y)
})
settings <- as.data.frame(data.table::rbindlist(settings))
# rm(settingsfiles)




#===============================================================================
# Variables common to all analyses that require input
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

# Output directory
out_dir = choose_directory(caption = "Select output folder",
                           default = dirname(settingsfiles)[[1]])



#===============================================================================
# Grouping variables
#===============================================================================
# Number of grouping variables
group_n <- tk_select.list(choices = c(NA, as.character(seq(3))),
                          preselect = NA,
                          multiple = FALSE,
                          title = "Select number of grouping variables")

# Remove pre-existing grouping vector
rm(grouping_vector)

# Create named vector of grouping variables (if no groupings selected, do stats on the unique sample IDs)
if (is.na(group_n)) {
  
  grouping_vector <- setNames(unique_id, unique_id_label)
  
} else {
  
  for (n in seq(as.numeric(group_n))) {
    
    # Select column
    var <- tk_select.list(choices = names(settings),
                          multiple = FALSE,
                          title = paste0("Select grouping variable ", n, "\n
                                       1 = x-axis, 2 = facet rows, 3 = facet columns"))
    
    # Create label
    label = dlgInput(message = paste0("Enter label for grouping variable ", n),
                     default = str_to_sentence(var))$res
    
    # Make vector
    vec <- setNames(var, label)
    
    # Add to the grouping vector
    if (exists("grouping_vector")) {
      grouping_vector <- c(grouping_vector, vec)
    } else {
      grouping_vector <- vec
    }
    
  }
  
}

grouping_vector




#===============================================================================
# Variables common to all analyses that can be set by DEFAULT
#===============================================================================
# Use default variables?
default <- tk_select.list(choices = c("Yes", "No"),
                          preselect = "No",
                          multiple = FALSE,
                          title = "Use default variables?\n
                          Unique sample ID = sample_id
                          Unique sample label = Sample ID
                          Target source = Settings file
                          Outlier threshold = 1.5
                          Remove outliers = Yes
                          Exclusion column = exclude
                          Hide non-significant p-values on plots = Yes
                          PDF dimensions = 10:13
                          Internal positive control target = NA\n
                          Comparative Ct variables:
                          Calibrator column = calibrator
                          Exclude HK if doesn't work for all samples = Yes\n
                          Standard curve variables:
                          Standard concentration column = standard_conc_ngul
                          Concentration unit = ng/µL
                          Sample volume column = vol_ul
                          Quantifiler experiment = No
                          Plot dependent variable log scale = No")
default <- ifelse(default == "Yes", TRUE, FALSE)


# Set variables
if (default) {
  
  unique_id = "sample_id"
  unique_id_label = "Sample ID"
  target_source = "Settings file"
  outlier_threshold = 1.5
  remove_outliers = TRUE
  exclusion_var = "exclude"
  hide_ns = TRUE
  y_log = FALSE
  pdf_height = 10
  pdf_width = 13
  
} else {

  # Unique ID
  unique_id = tk_select.list(choices = names(settings),
                             preselect = "sample_id",
                             multiple = FALSE,
                             title = "Select unique sample identifier column")
  
  # Unique ID label
  unique_id_label = dlgInput(message = "Enter unique ID label",
                             default = unique_id)$res
  
  
  # Target source
  target_source <- tk_select.list(choices = c("qPCR data file", "Settings file"),
                                  preselect = "Settings file",
                                  multiple = FALSE,
                                  title = "Select source to take target (gene) information from")
  
  # Tech rep outliers
  outlier_threshold = as.numeric(dlgInput(message = "Enter outlier threshold
                                        [abs(value - median(value)) > threshold * sd]",
                                          default = 1.5)$res)
  
  remove_outliers = tk_select.list(choices = c("Yes", "No"),
                                   preselect = "Yes",
                                   multiple = FALSE,
                                   title = "Remove outlier technical replicates from analysis?")
  remove_outliers <- ifelse(remove_outliers == "Yes", TRUE, FALSE)
  
  
  # Exclusions
  exclusion_var = tk_select.list(choices = c(NA, names(settings)),
                                 preselect = "exclude",
                                 multiple = FALSE,
                                 title = "Select column of excluded samples (NA if no exclusion column)")
  
  # Hide non-signficant p values?
  hide_ns <- tk_select.list(choices = c("Yes", "No"),
                            preselect = "Yes",
                            multiple = FALSE,
                            title = "Hide non-significant p-values on statistics plots?")
  hide_ns <- ifelse(hide_ns == "Yes", TRUE, FALSE)
  
  # Log scale y-axis
  y_log <- tk_select.list(choices = c("Yes", "No"),
                          preselect = "No",
                          multiple = FALSE,
                          title = "Log scale the dependent variable on plots (y-axis)?")
  y_log <- ifelse(y_log == "Yes", TRUE, FALSE)
  
  # PDF dimensions
  pdf_dimensions <- dlgInput(message = "Enter PDF dimensions (height:width)",
                             default = "10:13")$res
  
  pdf_height <- as.numeric(sub("\\:.*", "", pdf_dimensions))
  pdf_width <- as.numeric(sub(".*\\:", "", pdf_dimensions))
  
}



#===============================================================================
# Decide how to handle failed samples in housekeeping genes
#===============================================================================
# NC samples (used to exclude HK genes if they fail in a real sample)
options <- settings %>%
  unite("option", c(unique_id, unname(grouping_vector)), remove = FALSE, sep = " - ") %>%
  distinct(across(all_of(c(unique_id, unname(grouping_vector), "option"))))
  # dplyr::mutate(option = paste0(!!sym(unique_id), " (", !!sym(group_var), ")")) %>%
  # distinct(!!sym(unique_id), !!sym(group_var), option)

empty_samples = tk_select.list(choices = options$option,
                               multiple = TRUE,
                               title = "Select empty sample/s\n
                              Samples you expect NOT to amplify.
                              Conversely this is used to identify samples you expect to work
                              so that housekeeping genes can be excluded if a real sample fails")

if (identical(empty_samples, character(0))) {
  empty_samples <- NA
}

empty_samples <- options$sample_id[match(empty_samples, options$option)]




#===============================================================================
# Variables for comparative ct analyses
#===============================================================================
if (grepl("Comparative", analysis_method, fixed = TRUE)) {
  
  # Default method of calibrator selection
  calibration_method = tk_select.list(choices = c("Default",
                                                  "Calibrator samples",
                                                  "Lowest GOI ct",
                                                  "Highest GOI ct",
                                                  "Select samples"),
                                      preselect = "Default",
                                      multiple = FALSE,
                                      title = paste0("Choose a method to select calibrator samples \n
                                                  (Default uses mean Ct of callibrator samples as specified in the Settings file,
                                                     or if not available then the Ct of the sample with lowest gene of interest Ct)"))
  
  if(calibration_method == "Default" |
     calibration_method == "Calibrator samples") {
    
    
    if (default) {
      calibrator_var = "calibrator"
    } else {
      calibrator_var = tk_select.list(choices = names(settings),
                                      preselect = "calibrator",
                                      multiple = FALSE,
                                      title = "Select calibrator sample column")
    }
  } else { calibrator_var = NA }
  
  
  
  
  # Select calibrator samples
  if(calibration_method == "Select samples") {
    
    calibrator_samples = tk_select.list(choices = options$option,
                                        multiple = TRUE,
                                        title = "Select calibrator samples")
    calibrator_samples <- options$sample_id[match(calibrator_samples, options$option)]
    
  } else { calibrator_samples = NA }
  
  
  # Exclude a housekeeping target if it doesn't work for all samples
  if (default) {
    exclude_failed_hk <- TRUE
  } else {
    exclude_failed_hk = tk_select.list(choices = c("Yes", "No"),
                                       preselect = "Yes",
                                       multiple = FALSE,
                                       title = "Exclude housekeeping targets (genes) if they don't work for all samples")
    exclude_failed_hk <- ifelse(exclude_failed_hk == "Yes", TRUE, FALSE)
  }
  
}





#===============================================================================
# Standard curve variables
#===============================================================================
if (analysis_method == "Standard curve") {
  
  if (default) {
    
    standard_var = "standard_conc_ngul"
    conc_unit = "ng/µL"
    vol_var = "vol_ul"
    quantifiler = FALSE
    y_log = FALSE
    
  } else {
    
    # Standards
    standard_var = tk_select.list(choices = names(settings),
                                  multiple = FALSE,
                                  title = "Select standard concentration column")
    
    # Standard concentrations
    conc_unit <- dlgInput(message = "Enter concentration unit",
                          default = "ng/µL")$res
    
    # Source volume
    vol_var = tk_select.list(choices = names(settings),
                             multiple = FALSE,
                             title = "Select sample volume column")
    
    # Is this a quantifiler experiment?
    quantifiler = tk_select.list(choices = c("Yes", "No"),
                                 preselect = "No",
                                 multiple = FALSE,
                                 title = "Quantifiler experiment?")
    
    # Log scale y-axis
    y_log <- tk_select.list(choices = c("Yes", "No"),
                            preselect = "No",
                            multiple = FALSE,
                            title = "Log scale the dependent variable on plots (y-axis)?")
    y_log <- ifelse(y_log == "Yes", TRUE, FALSE)
    
  }
}






#===============================================================================
# Open a pdf
#===============================================================================
pdf(file.path(out_dir, "qpcr.pdf"),
    height = pdf_height, width = pdf_width, onefile = T)






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
    filter(well_short %in% well_short_vars) %>%
    dplyr::mutate(row = substr(well_short, 1, 1),
                  col = as.numeric(sub(".", "", well_short))) %>%
    dplyr::mutate(well = paste0(row, str_pad(col, 2, pad = "0"))) %>%
    relocate(c(well, well_short, well_n, row, col), .after = "filename")
    
  # Output
  return(list(headers = headers,
              ct = ct))
  
})


# Consolidate data from different plates
headers <- rbind.fill(lapply(qpcr_data, extract_sublist, sub = "headers"))
ct <- rbind.fill(lapply(qpcr_data, extract_sublist, sub = "ct"))

# Export header information
unlink(file.path(out_dir, "qpcr.xlsx"))
write.xlsx(headers,
           file = file.path(out_dir, "qpcr.xlsx"),
           sheetName = "headers", append = T, row.names = F)


# Export raw ct data
rm_excel_sheet(file = file.path(out_dir, "qpcr.xlsx"),
               sheetname = "ct")
write.xlsx(as.data.frame(ct),
           file = file.path(out_dir, "qpcr.xlsx"),
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
    select(filename, well, starts_with(c("reporter", "target"))) %>%
    pivot_longer(cols = -c(filename, well),
                 names_pattern = "(.*)\\d$",
                 names_to = ".value") %>%
    filter(complete.cases(.)) # https://stackoverflow.com/questions/75465900/using-pivot-longer-on-two-sets-of-columns/75465941#75465941
  
  ct <- ct %>%
    select(filename, well, reporter, ct) %>%
    left_join(target_map,
              by = c("filename", "well", "reporter")) %>%
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
                        title = "Select targets of interest")
  goi <- goi[order(match(goi, unique(ct$target)))]
  
  # Select housekeeping genes
  hk <- tk_select.list(choices = unique(ct$target)[!is.na(unique(ct$target))],
                       preselect = setdiff(sort(unique(ct$target)), goi),
                       multiple = TRUE,
                       title = "Select housekeeping target/s")
  hk <- hk[order(match(hk, unique(ct$target)))]
  
}


# IPC target
if (default) {
  ipc = NA
} else {
  ipc <- tk_select.list(choices = c(NA, unique(ct$target)[!is.na(unique(ct$target))]),
                        preselect = "NA",
                        multiple = FALSE,
                        title = "Select internal positive control (IPC) target\n
                      (included in PCR mix for accuracy of pipetting)")
}


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
              select(-starts_with(c("reporter", "target"))),
            by = c("filename", "well")) %>%
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
# Add group nomenclature
ct_option <- ct %>%
  left_join(options %>%
              select(!!sym(unique_id), option),
            by = unique_id) %>%
  relocate(option, .after = "well")

# Plot
print(
  ggplot(ct_option, aes(x = option, y = ct)) +
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
          # text = element_text(size = 15),
          plot.margin = unit(c(5.5, 5.5, 5.5, 50), "pt"))
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
  rm_excel_sheet(file = file.path(out_dir, "qpcr.xlsx"),
                 sheetname = "ipc_outlier")
  write.xlsx(as.data.frame(ipc_outliers),
             file = file.path(out_dir, "qpcr.xlsx"),
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
    dplyr::mutate(tech_outlier = abs(ct - median(ct, na.rm = T)) > outlier_threshold * sd(ct, na.rm = T)) %>%
    mutate(outlier_threshold = outlier_threshold) %>%
    relocate(c(tech_rep, outlier_threshold, tech_outlier), .after = ct) %>%
    mutate(tech_outlier = ifelse(is.na(tech_outlier), FALSE, tech_outlier)) %>%
    mutate(tech_outlier = factor(tech_outlier, levels = c(TRUE, FALSE))) %>%
    ungroup()
  
  # Output
  return(y)
  
})


# Export ct data
ct_export <- rbindlist(ct)
rm_excel_sheet(file = file.path(out_dir, "qpcr.xlsx"),
               sheetname = "outliers")
write.xlsx(ct_export,
           file = file.path(out_dir, "qpcr.xlsx"),
           sheetName = "outliers", append = T, row.names = F)


# Plot ct data with outliers
plot <- lapply(ct, function(x) {
  
  # # Manual
  # x = ct[[2]]
  
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
    # geom_jitter(aes(fill = tech_outlier), width = 0.1,height = 0, shape = 21, size = 3, alpha = 0.5) +
    geom_point(aes(fill = tech_outlier),
               position = position_jitter(seed = 42, width = 0.1, height = 0),
               shape = 21, size = 3, alpha = 0.5) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    labs(title = t,
         x = unique_id_label,
         y = "Threshold cycle (ct)",
         fill = "Outlier") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          # text = element_text(size = 15)
    ) +
    ggnewscale::new_scale_color() +
    geom_point(data = x %>%
                 dplyr::mutate(!!sym(exclusion_var) := 
                                 ifelse(is.na(!!sym(exclusion_var)), NA, TRUE)),
               aes(colour = !!sym(exclusion_var)),
               position = position_jitter(seed = 42, width = 0.1, height = 0),
               shape = 21, fill = NA, size = 5, stroke = 2, alpha = 0.75) +
    scale_colour_manual(values = c("TRUE" = "red"),
                        na.value = NA) +
    labs(colour = "Excluded")
  
  # Output
  return(plot)
  
})

# Wrap plots
print(
  wrap_plots(plot) +
    plot_annotation(title = "Technical replicate outliers",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
)







#===============================================================================
# Comparative ct (multiple housekeepers)
#===============================================================================
if (analysis_method == "Comparative ct (multiple housekeepers)") {
  
  # Prepare targets for comparative ct analysis
  source("./functions/comparative_prep_2.0.R")
  comparative_vars <- comparative_prep(ct = ct,
                                       settings = settings)
  
  # Assign variables
  for (n in names(comparative_vars)) {
    assign(n, comparative_vars[[n]])
  }
  
  # Comparative ct analysis
  # (N.B. use summary, as each unique sample has 1 value for goi and all hk)
  source("./functions/multiple_hk_2.0.R")
  comparative_ct <- multiple_hk(summary = summary)
  
  # Export relative expression data
  rm_excel_sheet(file = file.path(out_dir, "qpcr.xlsx"),
                 sheetname = "relative_expression")
  write.xlsx(as.data.frame(comparative_ct),
             file = file.path(out_dir, "qpcr.xlsx"),
             sheetName = "relative_expression", append = T, row.names = F)
  
  # Add group nomenclature
  comparative_ct_option <- comparative_ct %>%
    left_join(options %>%
                select(!!sym(unique_id), option),
              by = unique_id) %>%
    relocate(option, .after = unique_id)
  
  # Plot relative expression
  plots <- lapply(setNames(goi, goi), function(g) {
    
    # Bar plot
    plot <-
      ggplot(comparative_ct_option, aes(x = option, y = !!sym(g))) +
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
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            # text = element_text(size = 15),
            plot.margin = unit(c(5.5, 5.5, 5.5, 50), "pt"))
    
    # Output
    return(plot)
    
  })
  
  # Wrap the plots
  print(
    wrap_plots(plots) +
      plot_annotation(title = "Relative expression",
                      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  )
  
  # Prepare expression level data for stats
  quantification <- comparative_ct %>%
    select(any_of(c(unique_id, goi))) %>%
    pivot_longer(cols = any_of(goi),
                 names_to = "target",
                 values_to = "relative_expression") %>%
    split(f = as.factor(.$target))
  
  # Dependent variables vector for stats
  dependent_vector <- c("Relative expression" = "relative_expression")
  
}






#===============================================================================
# Comparative ct (2^-ddct)
#===============================================================================
if (analysis_method == "Comparative ct (2^-ddct)") {
  
  # Prepare targets for comparative ct analysis
  source("./functions/comparative_prep_2.0.R")
  comparative_vars <- comparative_prep(ct = ct,
                                       settings = settings)
  
  # Assign variables
  for (n in names(comparative_vars)) {
    assign(n, comparative_vars[[n]])
  }
  
  # Comparative ct analysis
  # (N.B. use summary, as each unique sample has 1 value for goi and all hk)
  source("./functions/twoddct_2.0.R")
  comparative_ct <- twoddct(summary = summary)
  
  # Plot relative expression
  plots <- lapply(setNames(goi, goi), function(g) {
    
    # # Manual
    # g = setNames(goi, goi)[[1]]
    
    # Bar plot
    plot <-
      ggplot(comparative_ct, aes(x = !!sym(unique_id), y = !!sym(g))) +
      geom_bar(aes(fill = calibrator), stat = "identity", colour = "black") +
      scale_fill_manual(values = c("TRUE" = "darkgrey", "FALSE" = "blue")) +
      geom_text(aes(label = round(!!sym(g), 2)), vjust=-0.3) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      { if (y_log) {
        scale_y_continuous(trans = "log",
                           labels = function(x) { round(x, digits = 0) })
      } else {
        scale_y_continuous(expand = expansion(mult = c(NA, 0.1)),
                           limits = c(0, NA))
      }} +
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
  
  # Prepare expression level data for stats
  quantification <- comparative_ct %>%
    select(any_of(c(unique_id, goi))) %>%
    pivot_longer(cols = any_of(goi),
                 names_to = "target",
                 values_to = "relative_expression") %>%
    split(f = as.factor(.$target))
  
  # Dependent variables vector for stats
  dependent_vector <- c("Relative expression" = "relative_expression")
  
}






#===============================================================================
# Standard curve
#===============================================================================
if (analysis_method == "Standard curve") {
  
  # Standard curve calculations
  source("./functions/standard_curve_2.0.R")
  sc <- standard_curve(ct = ct)
  
  # Prepare expression level data for stats
  quantification <- lapply(sc[["summary"]], function(x) {
    
    # Select columns
    y <- x %>%
      select(any_of(c(unique_id, "target", "predict_conc", "total_DNA")))
    
    # Output
    return(y)
    
  })
  
  # Dependent variables vector for stats
  dependent_vector <- c("predict_conc", "total_DNA")
  names(dependent_vector) <- c(paste0("Concentration (", conc_unit, ")"),
                               paste0("DNA yield (", gsub("/.*$", "", conc_unit), ")"))
  dependent_vector
  
}






#===============================================================================
# qPCR stats
#===============================================================================
source("./functions/qpcr_stats_2.2.R")
stats <- qpcr_stats(quantification = quantification,
                    grouping_vars = grouping_vector,
                    dependent_vars = dependent_vector,
                    y_log = y_log,
                    hide_ns = hide_ns)



#===============================================================================
# Close pdf
#===============================================================================
dev.off()


#===============================================================================
# Save data
#===============================================================================
# Save data
save(list = ls(),
     file = file.path(out_dir, "qpcr.RData"))

# # Load data
# lnames = load(file = file.path(out_dir, "qpcr.RData"))
# lnames


