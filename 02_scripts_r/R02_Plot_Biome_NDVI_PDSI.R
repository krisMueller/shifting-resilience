# ==============================================================================
# SCRIPT 02: Biome-Level NDVI & PDSI Time Series Analysis
# 
# DESCRIPTION:
# This script aggregates ecoregion-level time series data to create a 
# biome-wide analysis of vegetation productivity (NDVI) and drought (PDSI).
#
# KEY METHODOLOGY:
# 1. Weighted Aggregation: Calculates biome-level means by weighting each 
#    ecoregion's data by its total area.
# 2. Trend Analysis: Linear regression is applied to two distinct periods 
#    (1984–2004 vs 2005–2024) to assess shifts in climate sensitivity.
#
# INPUTS:
# - ../03_data_inputs/GEE03_[Region]_NDVI_PDSI_TimeSeries.csv
#
# OUTPUTS:
# - ../04_data_outputs/plots/R02_Biome_NDVI_TimeSeries_Plot.png

# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_dir   <- "../03_data_inputs"
output_dir  <- "../04_data_outputs/plots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 2. DATA CONFIGURATION --------------------------------------------------------

# Define region names and areas (Updated Feb 2026 with exact GEE Metadata)
region_df <- data.frame(
  region_name = c("coloradoPlateau", "columbianPlateau", "greatBasin", 
                  "wyomingBasin", "cascadeSierraMountains", "northernRockies", 
                  "southernRockies", "greatPlains"),
  # Exact totalArea from GEE FeatureCollection
  region_area = c(238707900639.61792, 
                  74875496572.15283, 
                  529363746173.2681, 
                  132467254362.48683,
                  119568996639.30588, 
                  177919269233.08115, 
                  191922951890.01584, 
                  292582156256.39624)
)

# Calculate total area and region weights
total_area <- sum(region_df$region_area)
region_df$weight <- region_df$region_area / total_area

# 3. DATA INGESTION ------------------------------------------------------------

# List CSV files matching the pattern
csv_files <- list.files(input_dir, pattern = "GEE03_.*NDVI_PDSI.*\\.csv", full.names = TRUE)

# Extract ecoregion name from the file name
extract_region_name <- function(file_name) {
  sub(".*GEE03_(.*)_NDVI_PDSI.*\\.csv", "\\1", basename(file_name))
}

# Read CSV and add region name and weight
read_and_add_region <- function(file_name, region_df) {
  data <- read.csv(file_name)
  region_name <- extract_region_name(file_name)
  
  # Standardize northern great plains naming if necessary
  if(region_name == "northern_great_Plains") region_name <- "greatPlains"
  
  data$region_name <- region_name
  data$weight <- region_df$weight[region_df$region_name == region_name]
  return(data)
}

# Define a custom color palette for decades
decade_palette <- c("1984-2004" = "#66c2a5", "2005-2024" = "#8da0cb")

# Combine all CSV files into a single dataset
all_data <- do.call(rbind, lapply(csv_files, read_and_add_region, region_df = region_df)) %>%
  select(-c(.geo, system.index))

# 4. ANALYSIS & PLOTTING -------------------------------------------------------

# Calculate the weighted means for all years
weighted_means <- all_data %>%
  group_by(year) %>%
  summarise(
    weighted_mean_ndvi = sum(mean_ndvi * weight) / sum(weight),
    weighted_mean_pdsi = sum(mean_period_pdsi * weight) / sum(weight),
    .groups = 'drop'
  ) %>%
  filter(year != 2012)  # Remove the year 2012

# Transform data to long format
df_long <- weighted_means %>%
  pivot_longer(cols = c(weighted_mean_ndvi, weighted_mean_pdsi),
               names_to = "variable",
               values_to = "value") %>%
  mutate(decade = case_when(
    year >= 1984 & year < 2005 ~ "1984-2004",
    year >= 2005 & year < 2024 ~ "2005-2024",
    TRUE ~ NA_character_
  ))

# Create the NDVI line plot
ndvi_plot <- ggplot(df_long %>% filter(variable == "weighted_mean_ndvi"), 
                    aes(x = year, y = value)) +
  geom_line(color = "#66c2a5", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(group = decade), color = "black") +
  theme_classic() +
  labs(y = "Mean late-season NDVI", x = "", title = "") +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))

# Create the bar plots for PDSI
period_plots <- df_long %>% 
  filter(variable == "weighted_mean_pdsi") %>%
  ggplot(aes(x = year, y = value, fill = value > 0)) +
  geom_bar(stat = "identity") +
  scale_color_manual(values = decade_palette) +
  scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "tan")) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1, 
             labeller = labeller(variable = NULL)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "grey") +
  labs(x = "Year", y = "Mean late-season PDSI") +
  geom_smooth(method = "lm", se = FALSE, aes(group = decade), color = "black") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_blank())

# Combine vertically
combined_plot <- ndvi_plot / period_plots 

# --- Regression Analysis ---

first_period_data <- weighted_means %>% filter(year >= 1984 & year <= 2004)
second_period_data <- weighted_means %>% filter(year >= 2005 & year <= 2024)

lm_model_first <- lm(weighted_mean_ndvi ~ weighted_mean_pdsi, data = first_period_data)
lm_model_second <- lm(weighted_mean_ndvi ~ weighted_mean_pdsi, data = second_period_data)

# Extract R-squared and p-value
rsq_first <- summary(lm_model_first)$r.squared
pval_first <- summary(lm_model_first)$coefficients[2, 4]
formatted_pval1 <- ifelse(pval_first < 0.001, "<0.001", format(round(pval_first, 2), nsmall = 2))

rsq_second <- summary(lm_model_second)$r.squared
pval_second <- summary(lm_model_second)$coefficients[2, 4]
formatted_pval2 <- ifelse(pval_second < 0.001, "<0.001", format(round(pval_second, 2), nsmall = 2))

y_limits <- range(c(first_period_data$weighted_mean_ndvi, second_period_data$weighted_mean_ndvi))

scatter_plot_first <- ggplot(first_period_data, aes(x = weighted_mean_pdsi, y = weighted_mean_ndvi)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "Mean late-season PDSI", y = "Mean late-season NDVI", title = "1984-2004") +
  annotate("text", x = -Inf, y = Inf, label = paste("R² =", round(rsq_first, 2), "\np =", formatted_pval1), 
           hjust = -0.1, vjust = 1.1, size = 5, color = "black") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
        strip.text = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylim(y_limits)

scatter_plot_second <- ggplot(second_period_data, aes(x = weighted_mean_pdsi, y = weighted_mean_ndvi)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "Mean late-season PDSI", y = "Mean late-season NDVI",  title = "2005-2024") +
  annotate("text", x = -Inf, y = Inf, label = paste("R² =", round(rsq_second, 2), "\np =", formatted_pval2), 
           hjust = -0.1, vjust = 1.1, size = 5, color = "black") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
        strip.text = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylim(y_limits)

scatter_combined <- scatter_plot_first + scatter_plot_second
final_plot <- combined_plot / scatter_combined

# 5. EXPORT --------------------------------------------------------------------

ggsave(file.path(output_dir, 'R02_Biome_NDVI_TimeSeries_Plot.png'),
       plot = final_plot, width = 6.5, height = 9, dpi = 300)