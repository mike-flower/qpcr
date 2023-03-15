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
# Import qPCR data and remove what's not needed
#===============================================================================
# Select qPCR RData file
qpcrdata <- tk_choose.files(caption = "Select qPCR Rdata file",
                            multi = FALSE,
                            filter = matrix(c("All files", "*",
                                              "RData", ".RData"),
                                            2, 2, byrow = TRUE))

# Load instability data
lnames = load(file = qpcrdata)
lnames

# Remove unwanted objects
if (grepl("Comparative", analysis_method, fixed = TRUE)) {
  rm(list = lnames[!lnames %in% c("comparative_ct", "ct", "goi", "calibrators", "analysis_method",
                                  "unique_id", "unique_id_label", "levels_id", "levels_target",
                                  "out.dir",
                                  "`%notin%`", "extract_sublist", "choose_directory")])
}

if (analysis_method == "Standard curve") {
  
  rm(list = lnames[!lnames %in% c("sc", "ct", "analysis_method",
                                  "unique_id", "unique_id_label", "levels_id", "levels_target",
                                  "out.dir",
                                  "`%notin%`", "extract_sublist", "choose_directory")])
  
}

rm(qpcrdata)

# Replace %notin% function
`%notin%` <- Negate(`%in%`)



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
# Variables
#===============================================================================
# Number of grouping variables
group_n <- tk_select.list(choices = seq(1:3),
                                preselect = 1,
                                multiple = FALSE,
                                title = "Select number of grouping variables")

# Grouping variable 1 (x-axis)
group1 <- tk_select.list(choices = names(settings),
                         multiple = FALSE,
                         title = "Select first grouping variable (x-axis)")

group1_label = dlgInput(message = "Enter label for first grouping variable",
                        default = group1)$res

# Grouping variable 2 (facet columns)
if (group_n > 1) {
  group2 <- tk_select.list(choices = names(settings),
                           multiple = FALSE,
                           title = "Select second grouping variable (facet columns)")
  group2_label = dlgInput(message = "Enter label for second grouping variable",
                          default = group2)$res
  
} else { group2 <- "" }

# Grouping variable 2 (facet columns)
if (group_n > 2) {
  group3 <- tk_select.list(choices = names(settings),
                           multiple = FALSE,
                           title = "Select third grouping variable (facet rows)")
  group3_label = dlgInput(message = "Enter label for third grouping variable",
                          default = group3)$res
  
} else { group3 <- "" }



#===============================================================================
# Open a pdf
#===============================================================================
# pdf(file.path(out.dir, "qpcr_stats.pdf"),
#     height = 10, width = 13, onefile = T)



#===============================================================================
# Format expression data AND SPLIT INTO A LIST BY TARGET
#===============================================================================
if (grepl("Comparative", analysis_method, fixed = TRUE)) {
  
  expression <- comparative_ct %>%
    select(any_of(c(unique_id, goi))) %>%
    pivot_longer(cols = any_of(goi),
                 names_to = "target",
                 values_to = "relative_expression") %>%
    split(f = as.factor(.$target))
  
}

if (analysis_method == "Standard curve") {
  
  expression <- sc[["summary"]]
  
}



#===============================================================================
# Annotate with required information from settings
#===============================================================================
# Select settings data
s <- settings %>%
  select(any_of(c(unique_id, group1, group2, group3))) %>%
  distinct()

# Group1 variable levels (in the order they appear in the settings file)
levels_group1 <- unique(settings[[group1]])

# Group2 variable levels
if (group_n > 1) {
  levels_group2 <- unique(settings[[group2]])
} else { levels_group2 <- "" }

# Group3 variable levels
if (group_n > 2) {
  levels_group3 <- unique(settings[[group3]])
} else { levels_group3 <- "" }



# Annotate relative expression data
expression <- lapply(expression, function(x) {
  
  # Annotate with settings
  y <- x %>%
    left_join(s, by = unique_id) %>%
    relocate(names(s)) %>%
    dplyr::mutate(!!sym(group1) := factor(!!sym(group1), levels = levels_group1)) %>%
    droplevels()
  
  if (group_n > 1) {
    y <- y %>%
      dplyr::mutate(!!sym(group2) := factor(!!sym(group2), levels = levels_group2)) %>%
      droplevels()
  }
  
  if (group_n > 2) {
    y <- y %>%
      dplyr::mutate(!!sym(group2) := factor(!!sym(group2), levels = levels_group2),
                    !!sym(group3) := factor(!!sym(group3), levels = levels_group3)) %>%
      droplevels()
  }
  
  # Output
  return(y)
  
})






#===============================================================================
# t-test
#===============================================================================
# Clear export spreadsheet before starting
unlink(file.path(out.dir, "qpcr_stats.xlsx"))

# T-test each target
plot <- lapply(expression, function(x) {
 
  # Remove NAs
  y <- x %>%
    filter(!is.na(!!sym(unique_id)),
           !is.na(!!sym(group1)))
  
  # Extract target
  t <- as.character(
    y %>%
      slice(1) %>%
      pull(target)
  )
  t
  
  # t-test
  ttest <- y %>%
    t_test(as.formula(paste0("relative_expression", " ~ ", group1))) %>%
    add_y_position()
  ttest
  
  # Export t-test
  if (file.exists(file.path(out.dir, "qpcr_stats.xlsx"))) {
    wb = loadWorkbook(file.path(out.dir, "qpcr_stats.xlsx"))
    removeSheet(wb, sheetName = t)
    saveWorkbook(wb, file.path(out.dir, "qpcr_stats.xlsx"))
    rm(wb)
  }
  write.xlsx(data.frame(ttest) %>%
               mutate(groups = as.character(groups)),
             file = file.path(out.dir, "qpcr_stats.xlsx"),
             sheetName = t, append = T, row.names = F)
  
  # Plot
  p <-
    ggplot(y, aes(x=!!sym(group1), y=relative_expression)) +
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
         x = group1_label,
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
# T-test each target
plot <- lapply(expression, function(x) {
  
  # Remove NAs
  y <- x %>%
    filter(!is.na(!!sym(unique_id)),
           !is.na(!!sym(group1)))
  
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
                     x = !!sym(group1), xlab = group1_label,
                     y = relative_expression, ylab = "Relative expression",
                     title = t,
                     package = "rcartocolor", palette = "Vivid") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
      error = function(e)
        ggbetweenstats(data = x,
                       x = !!sym(group1), xlab = group1_label,
                       y = relative_expression, ylab = "Relative expression",
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



#===============================================================================
# Close pdf
#===============================================================================
dev.off()


#===============================================================================
# Save data
#===============================================================================
# Save data
save(list = ls(),
     file = file.path(out.dir, "qpcr_stats.RData"))

# # Load data
# lnames = load(file = file.path(out.dir, "qpcr_stats.RData"))
# lnames










