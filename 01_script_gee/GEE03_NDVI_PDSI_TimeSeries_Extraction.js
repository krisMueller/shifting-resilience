/**
 * SCRIPT: GEE03_NDVI_PDSI_TimeSeries_Extraction
 * DESCRIPTION: Extracts annual (LATE-SEASON; July 15 - Septembr 30) mean NDVI and PDSI for ecoregions.
 * INPUT: Mesic Mask
 * OUTPUT: CSVs (GEE03_[Region]_NDVI_PDSI_TimeSeries.csv)
 */

// Load datasets
var riparian_mesic_mask = ee.Image("projects/ee-krismueller134/assets/Mesic/shiftingResilience_Assets/GEE01_Final_Riparian_Mesic_Mask");
var ecoregions = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/biome_clipped_reduced_eight_ecoregions");
var sagebrushBiome = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome")

// Load the NDVI image collection 
var ndviCollection = ee.ImageCollection('projects/ee-krismueller134/assets/Mesic/FinalLayers/PeriodMeanNdviCollection');

// Convert collection back to NDVI values
var divideBy10000 = function(image) {
  var dividedImage = image.divide(10000).toFloat();
  return dividedImage.setMulti(image.toDictionary(image.propertyNames()));
};

ndviCollection = ndviCollection.map(divideBy10000).map(function(img){
  return img;
});


// Load the PDSI image collection
var pdsiCollection = ee.ImageCollection("GRIDMET/DROUGHT").select('pdsi');

// Function to calculate annual mean
function calculateAnnualMean(collection, bandName) {
  return ee.ImageCollection(
    ee.List.sequence(1984, 2024).map(function(year) {
      var start = ee.Date.fromYMD(year, 1, 1);
      var end = ee.Date.fromYMD(year, 12, 31);
      var annualCollection = collection.filterDate(start, end);
      var annualMean = annualCollection.mean().select(bandName).set('year', year);
      return annualMean;
    })
  );
}

// Function to calculate period mean
function calculatePeriodMean(collection, bandName) {
  return ee.ImageCollection(
    ee.List.sequence(1984, 2024).map(function(year) {
      var start = ee.Date.fromYMD(year, 7, 15);
      var end = ee.Date.fromYMD(year, 9, 30);
      var periodCollection = collection.filterDate(start, end);
      var periodMean = periodCollection.mean().select(bandName).set('year', year);
      return periodMean;
    })
  );
}

// Calculate annual means
var annualNdvi = calculateAnnualMean(ndviCollection, 'mean_ndvi');
var periodPdsi = calculatePeriodMean(pdsiCollection, 'pdsi');

// Function to process data for a given region
function processRegion(region) {
  var combinedCollection = ee.FeatureCollection(                            
    ee.List.sequence(1984, 2024).map(function(year) {                         
      var ndviImage = annualNdvi.filter(ee.Filter.eq('year', year)).first().updateMask(riparian_mesic_mask) // <-------------- apply mask
      var pdsiPeriodImage = periodPdsi.filter(ee.Filter.eq('year', year)).first().updateMask(riparian_mesic_mask)

      var feature = ee.Feature(null, {
        'year': year,
        'mean_ndvi': ndviImage.reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: region.geometry(),
          scale: 30,
          maxPixels: 1e13,
          tileScale: 2
        }).get('mean_ndvi'),
        'mean_period_pdsi': pdsiPeriodImage.reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: region.geometry(),
          scale: 4600,
          maxPixels: 1e13,
          tileScale: 1
        }).get('pdsi')
      });
      return feature;
    })
  );
  
  return combinedCollection;
}

// Get unique nal3_name values
var names = ecoregions.aggregate_array('layer').distinct();
print(names)
// Function to format region names for export
function formatRegionName(name) {
  return name;
}

// Process each nal3_name region for valleys (riparian areas)
names.evaluate(function(names) {
  names.forEach(function(name) {
    var region = ecoregions.filter(ee.Filter.eq('layer', name));
    var formattedName = formatRegionName(name);
    
    // Process valleys (riparian areas)
    var valleysCollection = processRegion(region);
    Export.table.toDrive({
      collection: valleysCollection,
      description: 'GEE03_' + formattedName + '_NDVI_PDSI_TimeSeries',
      folder: "YOUR_GOOGLE_DRIVE_FOLDER",
      fileNamePrefix: "GEE03_" + formattedName + "_NDVI_PDSI_TimeSeries",
      fileFormat: "csv"
    });
  });
});