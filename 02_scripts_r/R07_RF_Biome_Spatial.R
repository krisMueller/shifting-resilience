# ==============================================================================
# SCRIPT 07: Random Forest - Biome-Level Spatial Gradient Approach
# 
# DESCRIPTION:
# This script executes a Random Forest model for the entire Biome (all regions 
# combined). 
#
# KEY METHODOLOGY:
# 1. Stratified Weighted Sampling: To prevent sampling bias, points are sampled 
#    proportionally based on the area size of each ecoregion (Total N = 12,000).
# 2. Spatial Mean: Temporal data is averaged for each unique location (.geo) 
#    to analyze spatial gradients across the entire biome (1986-2023).
#
# INPUTS:
# - ../03_data_inputs/GEE06_Covariates_Points_[YearRange].csv
#
# OUTPUTS:
# - ../04_data_outputs/stats/R07_SpatialGradient_Biome_ModelStats.txt
# - ../04_data_outputs/plots/R07_SpatialGradient_Biome_VIMP.png
# - ../04_data_outputs/plots/R07_SpatialGradient_Biome_PartialPlots.png
# - ../04_data_outputs/plots/R07_SpatialGradient_Biome_Error.png

# NOTE: Ensure your working directory is set to the location of this script.
# In RStudio: Session -> Set Working Directory -> To Source File Location

# NOTE: This script requires very long processing times. We recommend running
# this script on a server within a tmux session overnight, or on a always on PC 
# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_dir        <- "../03_data_inputs"
output_dir_plots <- "../04_data_outputs/plots"
output_dir_stats <- "../04_data_outputs/stats"

if (!dir.exists(output_dir_plots)) dir.create(output_dir_plots, recursive = TRUE)
if (!dir.exists(output_dir_stats)) dir.create(output_dir_stats, recursive = TRUE)

set.seed(123)
run_date    <- format(Sys.Date(), "%Y%m%d")
file_prefix <- paste0(run_date, "_R07_SpatialGradient_Biome_")

library(data.table)
library(dplyr)
library(randomForestSRC)
library(ggplot2)
library(tools)
library(cowplot)
library(car)

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

# 3. WEIGHTED STRATIFIED SAMPLING ----------------------------------------------

weights <- data.frame(
  region_name = c("coloradoPlateau", "columbianPlateau", "greatBasin",
                  "wyomingBasin", "cascadeSierraMountains", "northernRockies",
                  "southernRockies", "greatPlains"),
  region_area = c(238707900640, 74875496572, 529363746173, 132467254363,
                  119568996639, 177919269233, 191922951890, 292582156256)
)

total_points_per_year <- 1500 * length(unique(weights$region_name))

weights <- weights %>%
  mutate(
    area_proportion = region_area / sum(region_area),
    points_to_sample = round(area_proportion * total_points_per_year)
  )

diff <- sum(weights$points_to_sample) - total_points_per_year
if (diff != 0) {
  max_idx <- which.max(weights$points_to_sample)
  weights$points_to_sample[max_idx] <- weights$points_to_sample[max_idx] - diff
}

sampled_data <- list()
sorted_regions <- sort(unique(weights$region_name))

set.seed(123)

for (current_region in sorted_regions) {
  region_points <- weights %>% filter(region_name == current_region) %>% pull(points_to_sample)
  
  region_geo <- data %>%
    filter(region == current_region) %>%
    distinct(.geo) %>%
    slice_sample(n = region_points)
  
  for (current_year in 1986:2023) {
    temp_df <- data %>%
      filter(region == current_region) %>%
      filter(year == current_year) %>%
      filter(.geo %in% region_geo$.geo)
    
    sampled_data[[paste(current_region, current_year)]] <- temp_df
  }
}

# 4. SPATIAL AVERAGING & CONFIGURATION -----------------------------------------

filtered_biome <- bind_rows(sampled_data) %>%
  group_by(.geo) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = 'drop') %>%
  select(-c('percent_perennial', 'percent_annuals', 'percent_shrub', 'percent_tree', '.geo', 'year')) %>%
  as.data.table()

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

# 5. RANDOM FOREST MODEL -------------------------------------------------------

print(paste('processing biome model:'))

set.seed(123)
model <- rfsrc(mean_lateSeason_ndvi~.,
               data = filtered_biome,
               ntree = 1000,
               importance = 'permute',
               mc.cores = 1)

# 6. VARIABLE IMPORTANCE (VIMP) ------------------------------------------------

set.seed(123) 
subsampled_model <- subsample(model)

plot_data <- as.data.frame(t(plot.subsample(subsampled_model, alpha = 0.05, jknife = TRUE, show.plots = FALSE)$stats))
plot_data$variable_raw <- rownames(plot_data) 

plot_data_sorted <- plot_data %>% arrange(`50%`)
vimp_order_raw_ascending <- plot_data_sorted$variable_raw
partial_plot_order_raw <- rev(vimp_order_raw_ascending)

replace_with_labels <- function(var_name, plot_labels) {
  if(var_name %in% names(plot_labels)) return(plot_labels[var_name])
  return(var_name)
}

plot_data_sorted$variable_labeled <- sapply(plot_data_sorted$variable_raw, replace_with_labels, plot_labels)
plot_data_sorted$variable_labeled <- factor(plot_data_sorted$variable_labeled, levels = plot_data_sorted$variable_labeled)

p_top_vimps <- ggplot(plot_data_sorted, aes(x = `50%`, y = variable_labeled)) +
  geom_boxplot(aes(xmin = `2.5%`, xlower = `25%`, xmiddle = `50%`, xupper = `75%`, xmax = `97.5%`),
               stat = "identity", position = position_dodge(width = -0.8), width = 0.8, fill = "grey") +
  theme_minimal_hgrid(line_size = 0.3) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"), legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), panel.grid.minor = element_blank()) +
  labs(x = "Variable importance", y = "Covariates")

# 7. PARTIAL DEPENDENCE PLOTS --------------------------------------------------

{
  print('running partials')
  create_combined_partial_plots <- function(model, vars_in_order, plot_labels, scales, plot_title) {
    extract_partial_data <- function(partials_obj, sorted_vars, plot_labels) {
      partial_dfs <- list()
      for (var_name in sorted_vars) {
        var_index <- which(sapply(partials_obj$pData, function(p) p$xvar.names) == var_name)
        if (length(var_index) == 0) next
        p_data <- partials_obj$pData[[var_index]]
        if (length(unique(c(length(p_data$x.uniq), length(p_data$yhat), length(p_data$yhat.se)))) > 1) next
        df <- data.frame(Variable = var_name, X = p_data$x.uniq, Yhat = p_data$yhat, Yhat_SE = p_data$yhat.se)
        partial_dfs[[var_name]] <- df
      }
      final_df <- do.call(rbind, partial_dfs)
      vimp_sorted_labels <- sapply(sorted_vars, function(var) if (var %in% names(plot_labels)) plot_labels[[var]] else var, USE.NAMES = FALSE)
      final_df$Variable <- factor(final_df$Variable, levels = sorted_vars, labels = vimp_sorted_labels)
      return(final_df)
    }
    partials <- plot.variable(model, vars_in_order, partial = TRUE, show.plots = FALSE, sorted = FALSE, npts = 25)
    data <- extract_partial_data(partials, vars_in_order, plot_labels)
    data$period <- "1986-2023"
    
    p <- ggplot(data, aes(x = X, y = Yhat, color = period, fill = period)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = Yhat - Yhat_SE, ymax = Yhat + Yhat_SE), alpha = 0.2, linetype = 0) +
      facet_wrap(~Variable, scales = scales, ncol = 4) +
      theme_bw() +
      labs(title = plot_title, x = "", y = "Predicted mean late-season NDVI", color = "Period", fill = "Period") +
      scale_color_manual(values = c("1986-2023" = "black")) +
      scale_fill_manual(values = c("1986-2023" = "black")) +
      theme(strip.text = element_text(size = 10, face = "bold"), axis.text = element_text(size = 13),
            axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.position = "none", legend.text = element_text(size = 13), legend.title = element_text(size = 14))
    return(p)
  }
  print('Generating combined partial plot for all variables...')
  plot_all <- create_combined_partial_plots(model, partial_plot_order_raw, plot_labels, 'free_x', NULL)
}

# 8. EXPORT --------------------------------------------------------------------

writeLines(capture.output({print('Model output for biome vimp model with all variables:'); print(model)}),
           file.path(output_dir_stats, paste0(file_prefix, "ModelStats.txt")))

ggsave(file.path(output_dir_plots, paste0(file_prefix, "Error.png")), plot = plot(model), width = 11, height = 5, dpi = 300, bg = 'white')
ggsave(file.path(output_dir_plots, paste0(file_prefix, "VIMP.png")), plot = p_top_vimps, width = 9.5, height = 8, dpi = 300, bg = 'white')
ggsave(file.path(output_dir_plots, paste0(file_prefix, "PartialPlots.png")), plot = plot_all, width = 11, height = 9.5, dpi = 300, bg = 'white')