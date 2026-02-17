# ==============================================================================
# SCRIPT 04: Kendall's Tau Trend Distribution (Boxplots)
# 
# DESCRIPTION:
# This script visualizes the distribution of long-term NDVI trends (Kendall's Tau)
# derived from 200,000 sample points per ecoregion.
#
# KEY METHODOLOGY:
# 1. Data Aggregation: Combines large trend tables exported from GEE.
# 2. Period Comparison: Visualizes the shift in trend direction and magnitude
#    between P1 (1984-2004) and P2 (2005-2024) using side-by-side boxplots.
#
# INPUTS:
# - ../03_data_inputs/GEE05_Tau_Extract_[Region].csv
#
# OUTPUTS:
# - ../04_data_outputs/plots/R04_Tau_Boxplots.png

# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_dir   <- "../03_data_inputs"
output_dir  <- "../04_data_outputs/plots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)

# 2. DATA PROCESSING -----------------------------------------------------------

# Get list of files matching the GEE05 export pattern
file_list <- list.files(input_dir, pattern = "GEE05_Tau_Extract_.*\\.csv$", full.names = FALSE)

# Define explicit abbreviations from original script (plus NGP for consistency)
region_abbreviations <- list(
  "Columbian Plateau" = "CP",
  "Colorado Plateau" = "CoP",
  "Northern Great Plains" = "NGP"
)

extract_title <- function(filename) {
  # 1. Extract the raw region name from filename (e.g., "cascadeSierraMountains")
  raw_name <- str_extract(filename, "(?<=Extract_).*(?=\\.csv)")
  
  # 2. Handle specific naming for Great Plains to ensure it becomes "Northern Great Plains"
  if(raw_name == "greatPlains" || raw_name == "northern_great_Plains") {
    title <- "Northern Great Plains"
  } else {
    # 3. Convert CamelCase to Title Case with spaces (e.g., "cascadeSierraMountains" -> "Cascade Sierra Mountains")
    # This ensures the initials generator works for all 8 regions
    title <- str_replace_all(raw_name, "([a-z])([A-Z])", "\\1 \\2")
    title <- str_replace_all(title, "_", " ") # Handle any leftover underscores
    title <- str_to_title(title)
  }
  
  # 4. Get Abbreviation
  abbreviation <- region_abbreviations[[title]]
  
  # 5. Fallback: Generate initials if not in list (e.g., "Great Basin" -> "GB")
  if (is.null(abbreviation)) {
    abbreviation <- str_extract_all(title, "\\b\\w") %>% unlist() %>% paste(collapse = "")
  }
  return(abbreviation)
}

# Read and combine all data using extracted titles
all_data <- rbindlist(lapply(file_list, function(file) {
  # Read File
  df <- fread(file.path(input_dir, file))
  
  # Pivot Long
  df_long <- tidyr::pivot_longer(df, 
                                 cols = c(ndviTau_1984_2004, ndviTau_2005_2024),
                                 names_to = "Period",
                                 values_to = "NDVI_Tau")
  
  # Add Region Name
  df_long$Region <- extract_title(file)
  
  return(df_long)
}))

# Set factor levels for Period
all_data$Period <- factor(all_data$Period, 
                          levels = c("ndviTau_1984_2004", "ndviTau_2005_2024"), 
                          labels = c("1984-2004", "2005-2024"))

# 3. VISUALIZATION -------------------------------------------------------------

# Define custom colors for regions
region_colors <- c("#fcccab", "#cbd5e8", "#b3e2cd", "#f1c7e1", 
                   "#e6f5c9", "#fff2ae", "#f1e2cc", "#cccccc")

# Create the combined plot
combined_plot <- ggplot(all_data, aes(x = Region, y = NDVI_Tau, fill = Region)) +
  geom_hline(yintercept = 0, size=0.8, color = 'black', linetype = 'dashed') +
  geom_boxplot(outliers = FALSE) +
  facet_wrap(~ Period, scales = "free_x") +
  scale_fill_manual(values = region_colors) +
  labs(y = "NDVI tau-b") +
  ylim(-1, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.position = "none",  # Remove legend
        strip.text = element_text(size = 12))

# 4. EXPORT --------------------------------------------------------------------

# print(combined_plot)

ggsave(file.path(output_dir, 'R04_Tau_Boxplots.png'),
       plot = combined_plot, width=9, height=5, dpi=300)