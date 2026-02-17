# ==============================================================================
# SCRIPT 03: Ecoregion-Level NDVI & PDSI Time Series Analysis
# 
# DESCRIPTION:
# This script generates individual time-series and correlation plots for 
# each of the eight distinct ecoregions.
#
# KEY METHODOLOGY:
# 1. Iterative Processing: Loops through all input CSVs to process regions independently.
# 2. Correlation Analysis: Calculates R-squared and p-values for the relationship
#    between NDVI and PDSI for Period 1 (1984-2004) and Period 2 (2005-2024).
#
# INPUTS:
# - ../03_data_inputs/GEE03_[Region]_NDVI_PDSI_TimeSeries.csv
#
# OUTPUTS:
# - ../04_data_outputs/plots/R03_[Region]_NDVI_PDSI_Plot.png

# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_dir   <- "../03_data_inputs"
output_dir  <- "../04_data_outputs/plots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(tidyr)

# 2. DATA INGESTION & PROCESSING FUNCTIONS -------------------------------------

csv_files <- list.files(input_dir, pattern = "GEE03_.*NDVI_PDSI.*\\.csv", full.names = TRUE)
decade_palette <- c("1984-2004" = "#66c2a5", "2005-2024" = "#8da0cb")

extract_region_name <- function(file_name) {
  sub(".*GEE03_(.*)_NDVI_PDSI.*\\.csv", "\\1", basename(file_name))
}

create_plots_for_region <- function(file_name) {
  data <- read.csv(file_name)
  raw_region_name <- extract_region_name(file_name)
  
  # Standardize display name
  if(raw_region_name == "greatPlains" || raw_region_name == "northern_great_Plains") {
    region_title <- "Northern Great Plains"
  } else if(raw_region_name == "columbianPlateau") {
      region_title <- "Columbia Plateau"
  } else {
    region_title <- str_replace_all(raw_region_name, "([a-z])([A-Z])", "\\1 \\2")
    region_title <- str_to_title(region_title)
  }
  
  means <- data %>%
    group_by(year) %>%
    summarise(
      mean_ndvi = mean(mean_ndvi),
      mean_pdsi = mean(mean_period_pdsi),
      .groups = 'drop'
    ) %>%
    filter(year != 2012)
  
  df_long <- means %>%
    pivot_longer(cols = c(mean_ndvi, mean_pdsi), names_to = "variable", values_to = "value") %>%
    mutate(decade = case_when(
      year >= 1984 & year < 2005 ~ "1984-2004",
      year >= 2005 & year < 2024 ~ "2005-2024",
      TRUE ~ NA_character_
    ))
  
  # Plots
  ndvi_plot <- ggplot(df_long %>% filter(variable == "mean_ndvi"), aes(x = year, y = value)) +
    geom_line(color = "#66c2a5", linewidth = 0.8) +
    geom_smooth(method = "lm", se = FALSE, aes(group = decade), color = "black") +
    theme_classic() +
    labs(y = "Mean NDVI", x = "", title = region_title) +
    theme(legend.position = "bottom", axis.title = element_text(size = 12),
          axis.text = element_text(size = 12), legend.text = element_text(size = 12))
  
  period_plots <- df_long %>% filter(variable == "mean_pdsi") %>%
    ggplot(aes(x = year, y = value, fill = value > 0)) +
    geom_bar(stat = "identity") +
    scale_color_manual(values = decade_palette) +
    scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "tan")) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1, labeller = labeller(variable = NULL)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "grey") +  
    labs(x = "Year", y = "Late Season PDSI") +
    geom_smooth(method = "lm", se = FALSE, aes(group = decade), color = "black") +
    theme(legend.position = "none", axis.title = element_text(size = 12),
          axis.text = element_text(size = 12), strip.text = element_blank())
  
  # Regression Logic
  first_period_data <- means %>% filter(year >= 1984 & year <= 2004)
  second_period_data <- means %>% filter(year >= 2005 & year <= 2024)
  
  lm_model_first <- lm(mean_ndvi ~ mean_pdsi, data = first_period_data)
  lm_model_second <- lm(mean_ndvi ~ mean_pdsi, data = second_period_data)
  
  rsq_first <- summary(lm_model_first)$r.squared
  pval_first <- summary(lm_model_first)$coefficients[2, 4]
  formatted_pval1 <- ifelse(pval_first < 0.001, "<0.001", format(round(pval_first, 2), nsmall = 2))
  
  rsq_second <- summary(lm_model_second)$r.squared
  pval_second <- summary(lm_model_second)$coefficients[2, 4]
  formatted_pval2 <- ifelse(pval_second < 0.001, "<0.001", format(round(pval_second, 2), nsmall = 2))
  
  y_limits <- range(c(first_period_data$mean_ndvi, second_period_data$mean_ndvi))
  
  scatter_plot_first <- ggplot(first_period_data, aes(x = mean_pdsi, y = mean_ndvi)) +
    geom_point(alpha = 0.5) + theme_classic() +
    labs(x = "Mean late-season PDSI", y = "Mean late-season NDVI", title = "1984-2004") +
    annotate("text", x = -Inf, y = Inf, label = paste("R² =", round(rsq_first, 2), "\np =", formatted_pval1), 
             hjust = -0.1, vjust = 1.1, size = 4, color = "black") +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
          strip.text = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold")) + ylim(y_limits)
  
  scatter_plot_second <- ggplot(second_period_data, aes(x = mean_pdsi, y = mean_ndvi)) +
    geom_point(alpha = 0.5) + theme_classic() +
    labs(x = "Mean late-season PDSI", y = "Mean late-season NDVI",  title = "2005-2024") +
    annotate("text", x = -Inf, y = Inf, label = paste("R² =", round(rsq_second, 2), "\np =", formatted_pval2), 
             hjust = -0.1, vjust = 1.1, size = 4, color = "black") +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
          strip.text = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold")) + ylim(y_limits)
  
  return((ndvi_plot / period_plots) / (scatter_plot_first + scatter_plot_second))
}

# 3. GENERATE & EXPORT LOOP ----------------------------------------------------

for (file in csv_files) {
  region_name <- extract_region_name(file)
  
  # Generate Plot
  region_plot <- create_plots_for_region(file)
  print(paste("Generated plot for:", region_name))
  
  # Format filename
  if(region_name == "greatPlains") region_name <- "northernGreatPlains"
  if(region_name == "columbianPlateau") region_name <- "columbiaPlateau"
  save_name <- paste0("R03_", region_name, "_NDVI_PDSI_Plot.png")
  # print(region_plot)
  ggsave(file.path(output_dir, save_name), region_plot, width = 6.5, height = 9, dpi = 300)
}