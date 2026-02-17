# ==============================================================================
# SCRIPT 06: Random Forest - Ecoregion-Level Temporal Anomaly Approach
# 
# DESCRIPTION:
# This script executes Random Forest models to identify the drivers of 
# year-to-year NDVI variability across individual ecoregions.
#
# KEY METHODOLOGY:
# 1. Temporal Centering: For each location (.geo), the period mean is subtracted 
#    from each observation. This isolates anomalies (deviations from average).
# 2. Period Comparison: Data is split into Period 1 (1986–2004) and Period 2 
#    (2005–2023) to detect shifts in driver sensitivity over time.
#
# INPUTS:
# - ../03_data_inputs/GEE06_Covariates_Points_[YearRange].csv
#
# OUTPUTS (Generated for each ecoregion):
# - ../04_data_outputs/stats/R06_TemporalAnomaly_[Region]_P1/P2_ModelStats.txt
# - ../04_data_outputs/plots/R06_TemporalAnomaly_[Region]_VIMP_Comparison.png
# - ../04_data_outputs/plots/R06_TemporalAnomaly_[Region]_PartialPlots.png
# - ../04_data_outputs/plots/R06_TemporalAnomaly_[Region]_P1/P2_Error.png

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
file_prefix <- paste0(run_date, "_R06_TemporalAnomaly_")

library(data.table); library(dplyr); library(randomForestSRC); library(ggplot2);
library(tools); library(cowplot); library(tidyr); library(stringr)

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
plot_labels <- c("agDensity_05km" = "Irrigated Agriculture Density", "annual_huc8_runoff" = "Annual Runoff (mm)", "period_huc8_runoff" = "Late-Season Runoff (mm)", "elevation" = "Elevation (m)", "global_CO2" = "Annual Global CO2 (ppm)", "mean_Spring_precip" = "Spring Precip. (mm)", "mean_annual_huc8_swe" = "Annual SWE (mm)", "mean_lateSeason_frostFree" = "Late-Season Frost Free Days", "mean_lateSeason_maxTemp" = "Late-Season Max Temp. (°C)", "mean_lateSeason_pdsi" = "Late-Season PDSI", "mean_lateSeason_precip" = "Late-Season Precip. (mm)", "mean_lateSeason_z" = "Late-Season Z", "mean_prevWinter_maxTemp" = "Winter Max Temp. (°C)", "percent_annuals" = "Annual Forb/Grass Cover (%)", "percent_perennial" = "Perennial Forb/Grass Cover (%)", "percent_shrub" = "Shrub Cover (%)", "percent_tree" = "Tree Cover (%)", "seasonalWetlands_density_05km" = "Seasonal Wetland Density", "semiPermWetlands_density_05km" = "Semi-Permanent Wetland Density", "slope" = "Slope (%)", "tempWetlands_density_05km" = "Temporary Wetland Density")

for (region in unique_regions) {
  
  cat("\n--- Processing Region:", region, "---\n")
  
  # Naming logic: Use 'greatPlains' for filtering, but 'Northern Great Plains' for display
  region_label <- if(region == "greatPlains") "Northern Great Plains" else if (region == "columbianPlateau") "Columbia Plateau" else str_to_title(str_replace_all(region, "(?<!^)([A-Z])", " \\1"))
  
  region_data <- data %>% filter(region == !!region)
  
  # P1 Centering
  filtered_p1 <- region_data %>% filter(year <= 2004) %>% group_by(.geo) %>%
    mutate(across(where(is.numeric), ~ round(. - mean(., na.rm = TRUE), 8))) %>% ungroup() %>%
    select(-c('system:index', '.geo', 'year', 'region', 'percent_perennial', 'percent_annuals', 'percent_shrub', 'percent_tree', 'elevation', 'slope', 'seasonalWetlands_density_05km', 'semiPermWetlands_density_05km', 'tempWetlands_density_05km')) %>% as.data.table()
  
  # P2 Centering
  filtered_p2 <- region_data %>% filter(year > 2004, year != 2012) %>% group_by(.geo) %>%
    mutate(across(where(is.numeric), ~ round(. - mean(., na.rm = TRUE), 8))) %>% ungroup() %>%
    select(-c('system:index', '.geo', 'year', 'region', 'percent_perennial', 'percent_annuals', 'percent_shrub', 'percent_tree', 'elevation', 'slope', 'seasonalWetlands_density_05km', 'semiPermWetlands_density_05km','tempWetlands_density_05km')) %>% as.data.table()
  
  # Models
  set.seed(123); model_p1 <- rfsrc(mean_lateSeason_ndvi ~ ., data = filtered_p1, ntree = 1000, importance = 'permute')
  set.seed(123); model_p2 <- rfsrc(mean_lateSeason_ndvi ~ ., data = filtered_p2, ntree = 1000, importance = 'permute')
  
  # VIMP Comparison
  set.seed(123); subsampled_model_p1 <- subsample(model_p1)
  set.seed(123); subsampled_model_p2 <- subsample(model_p2)
  
  plot_data_p1 <- as.data.frame(t(plot.subsample(subsampled_model_p1, alpha = 0.05, jknife = TRUE, show.plots = FALSE)$stats))
  plot_data_p1$variable_raw <- rownames(plot_data_p1); plot_data_p1$period <- "p1"
  
  plot_data_p2 <- as.data.frame(t(plot.subsample(subsampled_model_p2, alpha = 0.05, jknife = TRUE, show.plots = FALSE)$stats))
  plot_data_p2$variable_raw <- rownames(plot_data_p2); plot_data_p2$period <- "p2"
  
  replace_with_labels <- function(var_name, plot_labels) { if(var_name %in% names(plot_labels)) return(plot_labels[var_name]); return(var_name) }
  plot_data_p1$variable <- sapply(plot_data_p1$variable_raw, replace_with_labels, plot_labels)
  plot_data_p2$variable <- sapply(plot_data_p2$variable_raw, replace_with_labels, plot_labels)
  
  vimp_order_data <- rbind(plot_data_p1, plot_data_p2) %>% group_by(variable) %>% summarise(mean_vimp = mean(`50%`)) %>% arrange(mean_vimp)
  vimp_plot_order_labeled <- vimp_order_data$variable
  partial_plot_order_raw <- rev((rbind(plot_data_p1, plot_data_p2) %>% mutate(variable = factor(variable, levels = vimp_plot_order_labeled)) %>% arrange(variable) %>% distinct(variable, variable_raw))$variable_raw)
  
  plot_data <- rbind(plot_data_p1, plot_data_p2) %>% mutate(variable = factor(variable, levels = vimp_plot_order_labeled), period = factor(period, levels = c("p2", "p1")))
  
  p_top_vimps <- ggplot(plot_data, aes(x = `50%`, y = variable, fill = period)) +
    geom_boxplot(aes(xmin = `2.5%`, xlower = `25%`, xmiddle = `50%`, xupper = `75%`, xmax = `97.5%`), stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    theme_minimal_hgrid(line_size = 0.3) +
    theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title = element_text(size = 18, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
    labs(title = region_label, x = "Variable importance", y = "Covariates", fill = "Period") +
    scale_fill_manual(values = c("p1" = "tan", "p2" = "lightblue"), labels = c("p1" = "1986-2004", "p2" = "2005-2023"), breaks = c("p1", "p2"))
  
  # Combined Partials
  {
    create_combined_partial_plots <- function(model_p1, model_p2, vars_in_order, plot_labels, scales, plot_title) {
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
        vimp_sorted_labels <- sapply(sorted_vars, function(var) if (var %in% names(plot_labels)) plot_labels[var] else var, USE.NAMES = FALSE)
        final_df$Variable <- factor(final_df$Variable, levels = sorted_vars, labels = vimp_sorted_labels)
        return(final_df)
      }
      partials_p1 <- plot.variable(model_p1, vars_in_order, partial = TRUE, show.plots = FALSE, sorted = FALSE, npts = 25)
      data_p1 <- extract_partial_data(partials_p1, vars_in_order, plot_labels); data_p1$period <- "1986-2004 (P1)"
      partials_p2 <- plot.variable(model_p2, vars_in_order, partial = TRUE, show.plots = FALSE, sorted = FALSE, npts = 25)
      data_p2 <- extract_partial_data(partials_p2, vars_in_order, plot_labels); data_p2$period <- "2005-2023 (P2)"
      combined_data <- rbind(data_p1, data_p2)
      ggplot(combined_data, aes(x = X, y = Yhat, color = period, fill = period)) +
        geom_line(linewidth = 1) + geom_ribbon(aes(ymin = Yhat - Yhat_SE, ymax = Yhat + Yhat_SE), alpha = 0.2, linetype = 0) +
        facet_wrap(~Variable, scales = scales, ncol = 4) + theme_bw() +
        labs(title = plot_title, x = "Covariate value deviation from period mean", y = "Predicted deviation in late-season NDVI from period mean", color = "Period", fill = "Period") +
        scale_color_manual(values = c("1986-2004 (P1)" = "sienna", "2005-2023 (P2)" = "cyan4")) +
        scale_fill_manual(values = c("1986-2004 (P1)" = "sienna", "2005-2023 (P2)" = "cyan4")) +
        theme(strip.text = element_text(size = 10, face = "bold"), axis.text = element_text(size = 13), axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.position = "bottom")
    }
    combined_plot_all <- create_combined_partial_plots(model_p1, model_p2, partial_plot_order_raw, plot_labels, 'free_x', region_label)
    }

  # EXPORT
  writeLines(capture.output(print(model_p1)), file.path(output_dir_stats, paste0(file_prefix, region, "_P1_ModelStats.txt")))
  writeLines(capture.output(print(model_p2)), file.path(output_dir_stats, paste0(file_prefix, region, "_P2_ModelStats.txt")))

  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_VIMP_Comparison.png")), plot = p_top_vimps, width = 10, height = 8)
  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_PartialPlots_Comparison.png")), plot = combined_plot_all, width = 11, height = 9.5)
  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_P1_Error.png")), plot = plot(model_p1), width = 11, height = 5)
  ggsave(file.path(output_dir_plots, paste0(file_prefix, region, "_P2_Error.png")), plot = plot(model_p2), width = 11, height = 5)
  
  # Cleanup
  objects_to_keep <- c("data", "unique_regions", "run_date", "plot_labels", "input_dir", "output_dir_plots", "output_dir_stats", "file_prefix")
  rm(list = setdiff(ls(), objects_to_keep)); gc()
}
