# ==============================================================================
# SCRIPT 10: Export Biome Model Sample Points (Shapefile)
# 
# DESCRIPTION:
# This script replicates the stratified weighted sampling logic used in 
# SCRIPT 05 (Temporal) and SCRIPT 07 (Spatial) to export the exact ~12,000 
# points used in the Random Forest models as a spatial shapefile.
#
# KEY METHODOLOGY:
# 1. Loads the raw covariate data (same inputs as analysis scripts).
# 2. Applies the exact same seed (set.seed(123)) and weighting logic.
# 3. Cleans and Parses the '.geo' column to extract Latitude/Longitude.
# 4. Exports an ESRI Shapefile for GIS visualization.
#
# INPUTS:
# - ../03_data_inputs/GEE06_Covariates_Points_[YearRange].csv
#
# OUTPUTS:
# - ../05_spatial_data/R10_Biome_Model_Points_12k.shp
# ==============================================================================

# 1. GLOBAL PATHS & SETTINGS ---------------------------------------------------

input_dir   <- "../03_data_inputs"
output_dir  <- "../05_spatial_data"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Global reproducibility seed (MUST MATCH R05 AND R07)
set.seed(123)

library(data.table)
library(dplyr)
library(sf)       
library(jsonlite) 
library(stringr)

# 2. DATA INGESTION ------------------------------------------------------------

message("Loading data...")
data1 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_1986_1995.csv'))
data2 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_1996_2004.csv'))
data3 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_2005_2014.csv'))
data4 <- fread(file.path(input_dir, 'GEE06_Covariates_Points_2015_2023.csv'))

# Combine just like the analysis scripts
data <- rbind(data1, data2, data3, data4)

# Cleanup to save RAM
rm(data1, data2, data3, data4); gc()

# 3. WEIGHTED STRATIFIED SAMPLING ----------------------------------------------
# NOTE: This section is identical to R05 and R07 to ensure the same points are selected.

message("Applying stratified sampling...")

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

# Adjust for rounding differences
diff <- sum(weights$points_to_sample) - total_points_per_year
if (diff != 0) {
  max_idx <- which.max(weights$points_to_sample)
  weights$points_to_sample[max_idx] <- weights$points_to_sample[max_idx] - diff
}

sampled_locations_list <- list()
sorted_regions <- sort(unique(weights$region_name))

# SAMPLING LOOP
for (current_region in sorted_regions) {
  region_points <- weights %>% filter(region_name == current_region) %>% pull(points_to_sample)
  
  # Identify the unique locations for this region using the Seed
  # This corresponds to 'region_geo' in R05/R07
  region_geo <- data %>%
    filter(region == current_region) %>%
    distinct(.geo, region) %>% # Keep region name for the shapefile attributes
    slice_sample(n = region_points)
  
  # Store these unique locations
  sampled_locations_list[[current_region]] <- region_geo
}

# 4. PARSE COORDINATES AND EXPORT ----------------------------------------------

message("Parsing coordinates and creating Shapefile...")

# Combine list into one dataframe of unique .geo strings
final_locations <- bind_rows(sampled_locations_list)

# Function to clean CSV quote escaping
clean_json <- function(geo_string) {
  # Replace double double-quotes ("") with single double-quotes (")
  cleaned <- gsub('""', '"', geo_string)
  # Sometimes strings start/end with a quote that needs removal if fread didn't catch it
  cleaned <- gsub('^"|"$', '', cleaned) 
  return(cleaned)
}

final_locations_parsed <- final_locations %>%
  rowwise() %>%
  mutate(
    # Clean the string before parsing
    geo_clean = clean_json(.geo),
    coords = list(fromJSON(geo_clean)$coordinates),
    lon = coords[1],
    lat = coords[2]
  ) %>%
  ungroup() %>%
  select(region, lon, lat) # Keep Region ID and coords

# Convert to SF object
points_sf <- st_as_sf(final_locations_parsed, coords = c("lon", "lat"), crs = 4326)

# ==============================================================================
# VERIFICATION: COUNT POINTS PER REGION
# ==============================================================================

# Calculate summary statistics
summary_stats <- points_sf %>%
  st_drop_geometry() %>%             # Remove spatial data for calculation
  group_by(region) %>%
  summarise(
    n_points = n(),
    percent_of_total = round((n() / nrow(points_sf)) * 100, 2)
  ) %>%
  arrange(desc(n_points))            # Sort by highest count

# Print to console
print("--- Point Counts per Ecoregion (Biome Model) ---")
print(summary_stats)

# Verify Total
print(paste("Total Points:", sum(summary_stats$n_points)))

# Save
output_path <- file.path(output_dir, "R10_Biome_Model_Points_12k.shp")
st_write(points_sf, output_path, delete_dsn = TRUE)

message(paste("Success! Shapefile exported to:", output_path))
message(paste("Total unique points exported:", nrow(points_sf)))
