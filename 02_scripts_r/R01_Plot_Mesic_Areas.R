# ==============================================================================
# SCRIPT 01: Mesic Area and Proportion Visualization
# 
# DESCRIPTION:
# This script processes the ecoregion-level statistics exported from GEE to 
# visualize the total area (hectares) and proportion (%) of mesic resources 
# within each ecoregion.
#
# KEY METHODOLOGY:
# 1. Data Processing: Converts raw hectare values to million-hectares for readability.
# 2. Visualization: Generates a faceted bar chart ranking ecoregions by 
#    mesic abundance.
#
# INPUTS:
# - ../03_data_inputs/GEE02_Mesic_Area_Stats.csv
#
# OUTPUTS:
# - ../04_data_outputs/plots/R01_Mesic_Area_Plot.svg
# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_file  <- "../03_data_inputs/GEE02_Mesic_Area_Stats.csv"
output_dir  <- "../04_data_outputs/plots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

library(ggplot2)
library(dplyr)
library(scales)
library(stringr)
library(tidyr)

# 2. DATA PROCESSING -----------------------------------------------------------

data <- read.csv(input_file) %>%
  mutate(layer = str_to_title(gsub("_", " ", layer)),
         mesic_area_m = mesic_area_ha / 1000000) %>%
  pivot_longer(
    cols = c(mesic_area_m, mesic_percent),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = case_when(
      metric == "mesic_area_m" ~ "Mesic Area (million ha)",
      metric == "mesic_percent" ~ "Mesic Proportion (%)"
    )
  )

# 3. VISUALIZATION -------------------------------------------------------------

# Create the plot
p <- ggplot(data, aes(x = reorder(layer, value), y = value)) +
  geom_col(fill = "black", width = 0.7) +
  geom_text(aes(label = ifelse(metric == "Mesic Area (million ha)",
                               sprintf("%.1f", value),
                               sprintf("%.1f%%", value))),
            vjust = -0.5,
            size = 4.5) +  
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Ecoregion",
       y = NULL,
       title = "Mesic Area and Proportion by Ecoregion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 20))
  )

# 4. EXPORT --------------------------------------------------------------------

ggsave(file.path(output_dir, "R01_Mesic_Area_Plot.png"), plot = p, width = 6.5, height = 7)