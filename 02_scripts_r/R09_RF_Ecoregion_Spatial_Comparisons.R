# ==============================================================================
# SCRIPT 09: Regional Comparison (Great Basin vs. Northern Great Plains)
# 
# DESCRIPTION:
# This script executes Random Forest models for two specific ecoregions
# (Great Basin and Northern Great Plains) using the Spatial Gradient (Mean) 
# approach. It aggregates the results to produce a direct comparison of 
# variable importance between these two key regions.
#
# KEY METHODOLOGY:
# 1. Spatial Mean: Temporal data is averaged for each unique location (.geo).
# 2. Targeted Filtering: Runs only for "greatPlains" and "greatBasin".
# 3. Comparative Visualization: Generates a "dumbbell" or connected boxplot 
#    to visualize how drivers shift between the two regions.
#
# INPUTS:
# - ../03_data_inputs/GEE06_Covariates_Points_[YearRange].csv
#
# OUTPUTS:
# - ../04_data_outputs/plots/R09_Comparison_GreatBasin_NorthernGreatPlains.png
# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_dir   <- "../03_data_inputs"
output_dir  <- "../04_data_outputs/plots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

set.seed(123)
run_date <- format(Sys.Date(), "%Y%m%d")

library(data.table)
library(dplyr)
library(randomForestSRC)
library(ggplot2)
library(tools)
library(cowplot)
library(tidyr)
library(stringr)

# 2. DATA INGESTION & PRE-PROCESSING -------------------------------------------

data1 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_1986_1995.csv'))
data2 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_1996_2004.csv'))
data3 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_2005_2014.csv'))
data4 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_2015_2023.csv'))
data <- rbind(data1, data2, data3, data4)

data$agDensity_05km <- data$agDensity_05km / 10000
data$seasonalWetlands_density_05km <- data$seasonalWetlands_density_05km / 10000
data$semiPermWetlands_density_05km <- data$semiPermWetlands_density_05km / 10000
data$tempWetlands_density_05km <- data$tempWetlands_density_05km / 10000

rm(data1, data2, data3, data4); gc()

# 3. COMPARISON CONFIGURATION --------------------------------------------------

# Only run the Great Plains and Great Basin models
unique_regions <- c("greatPlains", "greatBasin")

plot_labels <- c(
  "agDensity_05km" = "Irrigated agriculture density",
  "annual_huc8_runoff" = "Annual runoff (mm)",
  "period_huc8_runoff" = "late-season runoff runoff (mm)",
  "elevation" = "Elevation (m)",
  "global_CO2" = "Annual global CO2 (ppm)",
  "mean_Spring_precip" = "Spring precip. (mm)",
  "mean_annual_huc8_swe" = "Annual SWE (mm)",
  "mean_lateSeason_frostFree" = "Late-season frost free days",
  "mean_lateSeason_maxTemp" = "Late-season max temperature (C)",
  "mean_lateSeason_pdsi" = "Late-season PDSI",
  "mean_lateSeason_precip" = "Late-season precip. (mm)",
  "mean_lateSeason_z" = "Late-season Z",
  "mean_prevWinter_maxTemp" = "Winter max temp. (C)",
  "percent_annuals" = "Annual forbs/grass cover (%)",
  "percent_perennial" = "Perennial forbs/grass cover (%)",
  "percent_shrub" = "Shrub cover (%)",
  "percent_tree" = "Tree cover (%)",
  "seasonalWetlands_density_05km" = "Seasonal wetlands density",
  "semiPermWetlands_density_05km" = "Semi-Permanent wetlands density",
  "slope" = "Slope (%)",
  "tempWetlands_density_05km" = "Temporary wetlands density"
)

comparison_data_list <- list()

# 4. LOOP & MODELING -----------------------------------------------------------

for (region in unique_regions) {
  
  cat("\n--- Processing Region:", region, "---\n")
  
  region_title <- str_to_title(str_replace_all(region, "(?<!^)([A-Z])", " \\1"))
  region_title <- gsub("Great Plains", "Northern Great Plains", region_title)
  
  filtered_region_data <- data %>%
    filter(region == !!region) %>%
    group_by(.geo) %>%
    summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = 'drop') %>%
    select(-c('percent_perennial', 'percent_annuals', 'percent_shrub', 'percent_tree', '.geo', 'year')) %>%
    as.data.table()
  
  set.seed(123)
  model <- rfsrc(mean_lateSeason_ndvi~., data = filtered_region_data, ntree = 1000, importance = 'permute', mc.cores = 1)
  
  set.seed(123)
  subsampled_model <- subsample(model)
  plot_data <- as.data.frame(t(plot.subsample(subsampled_model, alpha = 0.05, jknife = TRUE, show.plots = FALSE)$stats))
  plot_data$variable_raw <- rownames(plot_data)
  
  plot_data_sorted <- plot_data %>% arrange(`50%`)
  replace_with_labels <- function(var_name, plot_labels) {
    if(var_name %in% names(plot_labels)) return(plot_labels[var_name])
    return(var_name)
  }
  
  plot_data_sorted$variable_labeled <- sapply(plot_data_sorted$variable_raw, replace_with_labels, plot_labels)
  plot_data_sorted$variable_labeled <- factor(plot_data_sorted$variable_labeled, levels = plot_data_sorted$variable_labeled)
  
  if(region %in% c("greatBasin", "greatPlains")) {
    temp_comp_data <- plot_data_sorted
    temp_comp_data$variable_labeled <- as.character(temp_comp_data$variable_labeled) 
    temp_comp_data$Region <- region
    comparison_data_list[[region]] <- temp_comp_data
  }
  
  # Cleanup
  objects_to_keep <- c("data", "unique_regions", "run_date", "plot_labels", "comparison_data_list", "input_dir", "output_dir")
  rm(list = setdiff(ls(), objects_to_keep))
  gc()
}

# 5. VISUALIZATION & EXPORT ----------------------------------------------------

if(length(comparison_data_list) > 1) {
  
  print("Generating Great Basin vs Northern Great Plains comparison plot...")
  comp_df <- do.call(rbind, comparison_data_list)
  
  comp_df$Region_Clean <- dplyr::recode(comp_df$Region, "greatBasin" = "Great Basin", "greatPlains" = "Northern Great Plains")
  
  order_stats <- comp_df %>%
    select(variable_labeled, Region_Clean, `50%`) %>%
    pivot_wider(names_from = Region_Clean, values_from = `50%`) %>%
    mutate(gap = abs(`Great Basin` - `Northern Great Plains`)) %>% 
    arrange(gap)
  
  comp_df$variable_labeled <- factor(comp_df$variable_labeled, levels = order_stats$variable_labeled)
  
  p_comparison <- ggplot() +
    geom_boxplot(data = comp_df,
                 aes(y = variable_labeled, xmin = `2.5%`, xlower = `25%`, xmiddle = `50%`, xupper = `75%`, xmax = `97.5%`,
                     fill = Region_Clean, color = Region_Clean),
                 stat = "identity", position = position_dodge(width = 0.25), width = 0.6, alpha = 0.7) + 
    scale_fill_manual(values = c("Great Basin" = "#C4318C", "Northern Great Plains" = "#A0D534")) + 
    scale_color_manual(values = c("Great Basin" = "#8a2262", "Northern Great Plains" = "#709524")) + 
    theme_minimal_hgrid(line_size = 0.8) +
    theme(legend.position = "top", legend.title = element_blank(), 
          legend.text = element_text(size = 14, face="bold"),
          axis.text = element_text(size = 14)) +
    labs(x = "Variable importance", y = "Covariates")
  
  # print(p_comparison)
  
  ggsave(file.path(output_dir, paste0(run_date, "_R09_SpatialGradient_Comparison_GreatBasin_NorthernGreatPlains.png")),
         plot = p_comparison, width = 12, height = 10, dpi = 300, bg = 'white')
}