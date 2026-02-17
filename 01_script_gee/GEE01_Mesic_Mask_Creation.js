/**
 * SCRIPT: GEE01_Mesic_Mask_Creation
 * DESCRIPTION: Generates the foundational mesic mask used in subsequent analysis.
 * OUTPUT: Asset (Mesic Mask Image)
 */


////////////////////////////////////
/// at least one year mesic Mask ///
////////////////////////////////////
{
  
var ndviCollection = ee.ImageCollection("projects/ee-krismueller134/assets/Mesic/FinalLayers/PeriodMeanNdviCollection")
  
  
// -------------------------------- Convert collection back to NDVI values -------------------
// revert back to ndvi values
var divideBy10000 = function(image) {
  var dividedImage = image.divide(10000).toFloat();
  return dividedImage.setMulti(image.toDictionary(image.propertyNames()));
};

// Map the function over the entire collection
var ndviCollection = ndviCollection.map(divideBy10000)

//////// 40 year number of years mean ndvi above 0.3 
var number_productive_40 = ndviCollection.map(function(img){
  return img.gte(0.3); // <------------ only years 'productive'
}).sum().rename('productive_years');

var atLeastOneYearMesic_mask = number_productive_40.gt(0).selfMask()
Map.addLayer(atLeastOneYearMesic_mask, null, 'atLeastOneYearMesic_mask', false)
  
}


//////////////////////////
/// environmental Mask ///
//////////////////////////
{

/*

This mask assigns unwanted pixels to 0. 
Make sure during implementaiton this mask is used in the updateMask function to keep everything that is equal to 1.

this mask removes:
- pixels that were classified by NASS CDL as cultivated ag  at least 3 times in the past 10 years (dataset 1997 - 2023, annual)
- pixels classified as open water or developement in NASS CDL 2023 (dataset 1997 - 2023)
- RAP mean tree cover < 10% in the last 5 years
- any pixel that was ever > 10% tree cover, within a MTBS fire polygon since 1984, and had a sig. drop in NDVI between any two years (> 0.2) (dataset RAP Tree cover)
  
*/

// ----------------------------- Define Overhead variables -----------------------------
var RAPCollection = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3')
var SagebrushBiome = ee.FeatureCollection('projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome')
var nass = ee.ImageCollection("USDA/NASS/CDL")
var MTBSfires = ee.FeatureCollection("USFS/GTAC/MTBS/burned_area_boundaries/v1")
var biome = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome")


// projection
var RapProj = RAPCollection.first().projection().getInfo()

// biome
var roi = biome.geometry().bounds()


// -------------- CROPLAND COVER MASK (MORE THAN 2 YEARS IN THE PAST 11) (NASS) ---------------


var cropLandcoverNass = nass
          .sort('system:time_start', false)
          .filter(ee.Filter.neq('system:index', '2024'))
          .limit(11) // Last 10 years
          .map(function(image) {return image.select('cultivated').eq(2)})
          .sum() // number of years pixel was cultivated
          .gt(0) // must be at least 1 year cultivated in the last 10
          .eq(0)
           

cropLandcoverNass = cropLandcoverNass.focalMode(7, 'square', 'pixels', 1)
// Map.addLayer(cropLandcoverNass, null, 'cdl cultivated cover', false)


// -------------------  Water  (if ever water) --------------

var waterMask = nass
          .map(function(image) {return image.select('cropland').eq(111)}) // water
          .sum()
          .gt(0)
          .eq(0)
// Map.addLayer(waterMask, null, 'water mask', false)


// -------------------  DEVELOPED COVER MASK NASS CDL (if ever developed) --------------

var developedCover = nass
          .map(function(image) {return image.select('cropland').eq([82, 121, 122, 123, 124])}) // developement
          .sum()
          .gt(0) 
          .reduce(ee.Reducer.max())
          .eq(0)

// Map.addLayer(developedCover, null, 'developed cover', false)


// // ----------- RAP TREE COVER was > 10% at least once taken out --------------

var treeCover = RAPCollection
          .map(function(image) {return image.select('TRE').gt(10)})
          .sum() 
          .gt(0) 
          .neq(1)

// Map.addLayer(treeCover, null, 'tree cover < 10%', false)          


// -------------- Fire Mask -------------

var burnsPolysMask = ee.Image().paint({
  featureCollection: MTBSfires,
  color: 'Ig_Date'
})

var fire_mask = RAPCollection
                      .map(function(image) {return image.select('TRE').gt(10)})
                      .sum()
                      .gt(0)
                      .updateMask(burnsPolysMask)
fire_mask = fire_mask.unmask().eq(0).clip(biome.geometry().bounds())

// Map.addLayer(fire_mask, null, 'forest fire mask', false)

// ------------------------------------ Create and Export --------------------------------

var environmental_mask = ee.Image(1).updateMask(cropLandcoverNass) // <------------- CDL NASS - Cultivation over time
                      .updateMask(developedCover) // <------------------- CDL NASS - Current developement (2023)
                      .updateMask(treeCover) // <------------------------ RAP - mean Tree cover in last 10 years
                      .updateMask(fire_mask) // <------------------------ RAP tre > 10% - forest fires
                      .updateMask(waterMask) // <------------------------ CDL NASS - water
                      .eq(1)
                      
Map.addLayer(environmental_mask, null, 'environmental_mask', false)


// Export.image.toAsset({
// 	image:environmental_mask.clip(biome.geometry().bounds()),
// 	description:'environmental_mask',
// 	assetId:'projects/ee-krismueller134/assets/Mesic/FinalLayers/environmental_mask',
// 	pyramidingPolicy:'mean',
// 	region:roi,
// 	crs:RapProj.crs,
// 	crsTransform:RapProj.transform,
// 	maxPixels:1e13,
// })

}


//////////////////////////
/// elevation Mask ///////
//////////////////////////
{
/*

This script removes elevation values above 90th percentile in Rocky mountain regions
                                          97th percentile in great Basin and cascade regions
                                          and no elevation clipping elswhere

*/

// variables
// Define global variables
var dem = ee.Image("USGS/3DEP/10m");
var ecoregions = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/biome_clipped_reduced_eight_ecoregions");

// Function to clip DEM to percentile and mask based on region
function clipToPercentile(feature) {
  var clipped = dem.clip(feature);
  var regionName = feature.get('layer');
  
  // Use ee.String().match() for string comparison instead of getInfo()
  var isRockies = ee.String(regionName).match('.*Rockies').length().gt(0);
  var isBasinOrSierra = ee.String(regionName).match('greatBasin|cascadeSierra').length().gt(0);
  
  return ee.Algorithms.If({
    condition: isRockies,
    trueCase: maskByPercentile(clipped, feature, 90),
    falseCase: ee.Algorithms.If({
      condition: isBasinOrSierra,
      trueCase: maskByPercentile(clipped, feature, 97),
      falseCase: clipped
    })
  });
}

// Helper function to mask by percentile
function maskByPercentile(image, feature, percentileValue) {
  var percentile = image.reduceRegion({
    reducer: ee.Reducer.percentile([percentileValue]),
    geometry: feature.geometry(),
    scale: 10,
    maxPixels: 1e13
  });
  var threshold = ee.Number(percentile.values().get(0));
  return image.updateMask(image.lt(threshold));
}

// Create the mosaic without using evaluate()
var processedImages = ecoregions.map(function(feature) {
  return clipToPercentile(feature);
});

var mosaic = ee.ImageCollection(processedImages).mosaic();
var elevation_mask = ee.Image(0).where(mosaic.gt(0), 1).selfMask();

// Add the layer to the map
Map.addLayer(elevation_mask, null, 'elevation_mask', false);

}


/////////////////////////////
/// woody trends Mask ///////
/////////////////////////////
{


var tree_trend = RAPCollection.map(function(img) {
  return img.select('TRE')
}).reduce(ee.Reducer.kendallsCorrelation()).select('TRE_tau');

var shrub_trend = RAPCollection.map(function(img) {
  return img.select('SHR')
}).reduce(ee.Reducer.kendallsCorrelation()).select('SHR_tau');

var n = ee.Image(38);
var mk_z_tree = ee.Image.constant(3).multiply(tree_trend).multiply(n.multiply(n.subtract(1)).sqrt())
              .divide(ee.Image.constant(2).multiply(((ee.Image.constant(2).multiply(n)).add(5)).sqrt()))
              .rename('Z_tree');
var mk_z_shrub = ee.Image.constant(3).multiply(shrub_trend).multiply(n.multiply(n.subtract(1)).sqrt())
              .divide(ee.Image.constant(2).multiply(((ee.Image.constant(2).multiply(n)).add(5)).sqrt()))
              .rename('Z_shrub');

var woody_mask = ee.Image(1).updateMask(mk_z_tree.abs().lt(1.96))
                          .updateMask(mk_z_shrub.abs().lt(1.96)).eq(1)
                          
Map.addLayer(woody_mask, null, 'woody_mask (no sig woody change 1986-2023', false)

}


/////////////////////////////
/// riparian Mask ///////////
/////////////////////////////
{
  
var vbet = ee.Image("projects/ee-krismueller134/assets/shiftingResilience_Assets/MERGED_VBET_RASTER_20240716"),
    nass = ee.ImageCollection("USDA/NASS/CDL"),
    biome = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome"),
    bps = ee.Image("projects/ee-krismueller134/assets/shiftingResilience_Assets/LANDFIRE_BPS");

var riparian_mask = ee.Image(0).where(bps.eq(7), 1)
                      .where(vbet.gte(0.6), 1)
                      .where(nass.map(function(image) {return image.select('cropland').eq(195)}).sum().gt(0), 1) // ----> CLD layer
                      .selfMask()
Map.addLayer(riparian_mask, null, 'riparian_mask', false)   
  
}

//////////////////////////////////////////////////
/// putting it all together Final Mask ///////////
//////////////////////////////////////////////////
{
  
var trend_palette = {palette:['red', 'white', 'green'], min:-1, max:1}
  
  
var final_riparian_mesic = ee.Image(1).updateMask(atLeastOneYearMesic_mask)
                                      .updateMask(environmental_mask)
                                      .updateMask(elevation_mask)
                                      .updateMask(woody_mask)
                                      .updateMask(riparian_mask)

Map.addLayer(final_riparian_mesic, {palette:'red'}, 'final riparian mesic mask')


Map.setOptions("HYBRID")

// Export mask
Export.image.toAsset({
	image:final_riparian_mesic,
	description:'GEE01_Final_Riparian_Mesic_Mask',
	assetId:'YOUR_PROJECT_ASSET_FOLDER',
	region:biome.geometry().bounds(),
	scale:30,
	crs:'EPSG:4326',
	maxPixels:1e13,
})


}