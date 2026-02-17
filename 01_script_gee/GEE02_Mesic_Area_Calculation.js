/**
 * SCRIPT: GEE02_Mesic_Area_Calculation
 * DESCRIPTION: Calculates mesic area statistics per ecoregion.
 * INPUT: Asset from GEE01
 * OUTPUT: CSV (GEE02_Mesic_Area_Stats.csv)
 */

// Load datasets
var ecoregions = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/biome_clipped_reduced_eight_ecoregions");
var mesic = ee.Image("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE01_Final_Riparian_Mesic_Mask");
var biome = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome");

// Ensure the mesic mask has values of 1 for mesic areas and 0 elsewhere
var binaryMesic = mesic.gt(0);

// Function to calculate mesic area within a feature
var calculateMesicArea = function(feature) {
  // Clip the mesic layer to the current feature
  var clippedMesic = binaryMesic.clip(feature.geometry());
  
  // Calculate area of mesic within the feature
  var mesicArea = clippedMesic.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: feature.geometry(),
    scale: 30, // Adjust scale based on the resolution of your mesic data
    maxPixels: 1e13
  });
  
  // Get the area value (convert from square meters to hectares)
  var areaValue = ee.Number(mesicArea.values().get(0)).divide(10000); // Convert to hectares
  
  // Return the feature with the calculated area
  return feature.set({
    'mesic_area_ha': areaValue,
    'total_area_ha': feature.geometry().area().divide(10000),
    'mesic_percent': areaValue.divide(feature.geometry().area().divide(10000)).multiply(100)
  });
};

// Map the function over the ecoregions feature collection
var ecoregionsWithArea = ecoregions.map(calculateMesicArea);

// Calculate total mesic area for the entire biome
var biomeClippedMesic = binaryMesic.clip(biome.geometry());
var biomeMesicArea = biomeClippedMesic.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: biome.geometry(),
  scale: 30,
  maxPixels: 1e13
});

// Calculate total biome area
var biomeTotalArea = biome.geometry().area().divide(10000); // Convert to hectares

// Get the mesic area value for the entire biome
var biomeMesicAreaValue = ee.Number(biomeMesicArea.values().get(0)).divide(10000); // Convert to hectares

// Calculate proportion of mesic area in the biome
var biomeMesicProportion = biomeMesicAreaValue.divide(biomeTotalArea).multiply(100);

// Print results
print('Ecoregions with mesic area:', ecoregionsWithArea);
print('Total Biome Area (hectares):', biomeTotalArea);
// print('Total Mesic Area in Biome (hectares):', biomeMesicAreaValue);
// print('Proportion of Mesic Area in Biome (%):', biomeMesicProportion);

// Export the results to a CSV file in the 'EE_Exports' folder
Export.table.toDrive({
  collection: ecoregionsWithArea,
  description: 'GEE02_Mesic_Area_Stats',
  folder: 'YOUR_GOOGLE_DRIVE_FOLDER',
  fileFormat: 'CSV'
});

// Optional: You could also export the biome-level mesic information
// Uncomment if you want to create a separate export for biome-level data
Export.table.toDrive({
  collection: ee.FeatureCollection([
    ee.Feature(null, {
      'total_biome_area_ha': biomeTotalArea,
      'total_mesic_area_ha': biomeMesicAreaValue,
      'mesic_proportion_percent': biomeMesicProportion
    })
  ]),
  description: 'GEE02_Biome_Mesic_Total',
  folder: 'YOUR_GOOGLE_DRIVE_FOLDER',
  fileFormat: 'CSV'
});