/**
 * SCRIPT: GEE04_Sample_Points_Extraction
 * DESCRIPTION: Generates stratified random samples (200k points) per ecoregion.
 * OUTPUT: Assets (FeatureCollections of points)
 */

var ecoregions = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/biome_clipped_reduced_eight_ecoregions")
var biome = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome")
var riparian_mesic_mask = ee.Image("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE01_Final_Riparian_Mesic_Mask")

var mesicMask = ee.Image(1).updateMask(riparian_mesic_mask) // <--------------- apply mask **************
Map.addLayer(mesicMask)


var ecoregions = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/biome_clipped_reduced_eight_ecoregions");
var regions = ecoregions.aggregate_array('layer').getInfo(); 

function processRegion(regionName) {
  var regionToClip = ecoregions.filter(ee.Filter.eq('layer', regionName));
  var maskImg = ee.Image().paint({
    featureCollection: regionToClip, 
    color: 'fid'
  });

  var regionPoints = mesicMask.updateMask(maskImg).stratifiedSample({
    numPoints: 200000,
    classBand: 'constant',
    region: regionToClip.geometry(),
    scale: 30,
    geometries: true,
    seed: 123,
    tileScale: 2
  }).map(function(feature) {
    return feature.set('region', regionName);
  });

  Export.table.toAsset({
    collection: regionPoints,
    description: "GEE04_200k_Points_" + regionName,
    assetId: "YOUR_PROJECT_ASSET_FOLDER_" + regionName,
  });

  Map.addLayer(maskImg, null, 'masked_' + regionName + '_img', false);

  return regionPoints;
}

var allRegionsPoints = regions.map(function(regionName) {
  return processRegion(regionName);
});