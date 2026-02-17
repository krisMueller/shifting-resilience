# ==============================================================================
# SCRIPT 08: Random Forest - Ecoregion-Level Spatial Gradient Approach
# 
# DESCRIPTION:
# This script executes Random Forest models (via randomForestSRC) to identify 
# long-term environmental drivers of NDVI across multiple ecoregions.
#
# KEY METHODOLOGY:
# 1. Spatial Mean Approach: The script collapses the time-series by averaging 
#    temporal covariate data for each unique geographic location (.geo). This 
#    isolates spatial gradients across the landscape (1986-2023).
# 2. Variable Importance (VIMP): Uses permutation-based importance with 
#    subsampling to calculate 95% confidence intervals for driver significance.
# 3. Partial Dependence: Generates marginal effect plots for all covariates.
#
# INPUTS:
# - ../03_data_inputs/GEE06_Covariates_Points_[YearRange].csv
#
# OUTPUTS (Generated for each ecoregion):
# - ../04_data_outputs/stats/R08_SpatialGradient_[Region]_ModelStats.txt
# - ../04_data_outputs/plots/R08_SpatialGradient_[Region]_VIMP.png
# - ../04_data_outputs/plots/R08_SpatialGradient_[Region]_PartialPlots.png
# - ../04_data_outputs/plots/R08_SpatialGradient_[Region]_Error.png

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
file_prefix <- paste0(run_date, "_R08_SpatialGradient_")

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

# 3. REGIONAL ANALYSIS LOOP ----------------------------------------------------

unique_regions <- unique(data$region)

plot_labels <- c(
  "agDensity_05km" = "Irrigated agriculture density",
  "annual_huc8_runoff" = "Annual runoff (mm)",
  "period_huc8_runoff" = "late-season runoff (mm)",
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

for (region in unique_regions) {
  
  cat("\n--- Processing Region:", region, "---\n")
  
  # Ensure NGP naming consistency for plots
  region_label <- if(region == "greatPlains") "Northern Great Plains" else if (region == "columbianPlateau") "Columbia Plateau" else str_to_title(str_replace_all(region, "(?<!^)([A-Z])", " \\1"))
  
  # Filter & Average
  filtered_region_data <- data %>%
    filter(region == !!region) %>%
    group_by(.geo) %>%
    summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = 'drop') %>%
    select(-c('percent_perennial', 'percent_annuals', 'percent_shrub', 'percent_tree', '.geo', 'year')) %>%
    as.data.table()
  
  # Model
  set.seed(123)
  model <- rfsrc(mean_lateSeason_ndvi ~ ., data = filtered_region_data, ntree = 1000, importance = 'permute', mc.cores = 1)
  
  # VIMP logic
  set.seed(123); subsampled_model <- subsample(model)
  plot_data <- as.data.frame(t(plot.subsample(subsampled_model, alpha = 0.05, jknife = TRUE, show.plots = FALSE)$stats))
  plot_data$variable_raw <- rownames(plot_data)
  plot_data_sorted <- plot_data %>% arrange(`50%`)
  vimp_order_raw_ascending <- plot_data_sorted$variable_raw
  partial_plot_order_raw <- rev(vimp_order_raw_ascending)
  
  replace_with_labels <- function(var_name, plot_labels) { if(var_name %in% names(plot_labels)) return(plot_labels[var_name]); return(var_name) }
  plot_data_sorted$variable_labeled <- sapply(plot_data_sorted$variable_raw, replace_with_labels, plot_labels)
  plot_data_sorted$variable_labeled <- factor(plot_data_sorted$variable_labeled, levels = plot_data_sorted$variable_labeled)
  
  p_top_vimps <- ggplot(plot_data_sorted, aes(x = `50%`, y = variable_labeled)) +
    geom_boxplot(aes(xmin = `2.5%`, xlower = `25%`, xmiddle = `50%`, xupper = `75%`, xmax = `97.5%`), stat = "identity", position = position_dodge(width = -0.8), width = 0.8, fill = "grey") +
    theme_minimal_hgrid(line_size = 0.3) +
    theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title = element_text(size = 18, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
    labs(x = "Variable importance", y = "Covariates", title = region_label)
  
  # Partials
  {
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
      data_partial <- extract_partial_data(partials, vars_in_order, plot_labels); data_partial$period <- "1986-2023"
      ggplot(data_partial, aes(x = X, y = Yhat, color = period, fill = period)) +
        geom_line(linewidth = 1) + geom_ribbon(aes(ymin = Yhat - Yhat_SE, ymax = Yhat + Yhat_SE), alpha = 0.2, linetype = 0) +
        facet_wrap(~Variable, scales = scales, ncol = 4) + theme_bw() + labs(title = plot_title, x = "", y = "Predicted mean late-season NDVI") +
        scale_color_manual(values = c("1986-2023" = "black")) + scale_fill_manual(values = c("1986-2023" = "black")) +
        theme(strip.text = element_text(size = 10, face = "bold"), axis.text = element_text(size = 13), axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.position = "none")
    }
    plot_all <- create_combined_partial_plots(model, partial_plot_order_raw, plot_labels, 'free_x', region_label)
    }
  
  # Diagnostics & Export
  
  writeLines(capture.output(print(model)), file.path(output_dir_stats, paste0(file_prefix, region, "_ModelStats.txt")))
  
  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_Error.png")), plot = plot(model), width = 7, height = 5)
  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_VIMP.png")), plot = p_top_vimps, width = 9.5, height = 8)
  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_PartialPlots.png")), plot = plot_all, width = 11, height = 9.5)
  
  # Cleanup
  graphics.off()
  objects_to_keep <- c("data", "unique_regions", "run_date", "plot_labels", "input_dir", "output_dir_plots", "output_dir_stats", "file_prefix")
  rm(list = setdiff(ls(), objects_to_keep)); gc()
}