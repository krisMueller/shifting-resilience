/**
 * SCRIPT: GEE05_Tau_Trend_Export
 * DESCRIPTION: Calculates Kendall's Tau trends for 200k points.
 * INPUT: Assets from GEE04
 * OUTPUT: CSVs (GEE05_Tau_Extract_[Region].csv)
 */

var basinOutline = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-121.84314234849965, 43.1432005595187],
          [-120.54675562974965, 42.181843724639236],
          [-120.53576930162465, 41.395476420373576],
          [-121.27185328599965, 41.11466032171771],
          [-121.22790797349965, 39.751670101501205],
          [-120.09631617662465, 38.12820283893192],
          [-118.47033961412465, 35.8557972405676],
          [-117.50354273912465, 35.552471943814005],
          [-115.53699000474965, 35.668585516560704],
          [-113.80115016099965, 34.98738630343024],
          [-112.62561305162465, 35.57928234029463],
          [-112.31799586412465, 37.346271595231784],
          [-111.04358180162465, 39.540183005430066],
          [-110.41736109849965, 43.89207164811022],
          [-111.68078883287465, 44.75660109129916],
          [-114.14172633287465, 43.89998838685847],
          [-115.30627711412465, 44.55340769831797],
          [-116.58069117662465, 45.01347504647374],
          [-121.30481227037465, 44.12911591938823]]]);

// Define tables
var tables = {
    cascadeSierraMountains: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_cascadeSierraMountains"),
    coloradoPlateau: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_coloradoPlateau"),
    columbianPlateau: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_columbianPlateau"),
    greatBasin: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_greatBasin"),
    greatPlains: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_greatPlains"), 
    northernRockies: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_northernRockies"),
    southernRockies: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_southernRockies"),
    wyomingBasin: ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE04_200k_Points_wyomingBasin")
};

Map.addLayer(tables.northernRockies)
// calculate period trends
var ndviCollection = ee.ImageCollection("projects/ee-krismueller134/assets/Mesic/FinalLayers/PeriodMeanNdviCollection");

// Function to divide by 10000
var divideBy10000 = function(image) {
  var dividedImage = image.divide(10000).toFloat();
  return dividedImage.setMulti(image.toDictionary(image.propertyNames()));
};

// Apply division and filter out 2012
var ndviCollection = ndviCollection.map(divideBy10000).filter(ee.Filter.neq('year', 2012));

// Function to add time band
var addTime = function(image) {
  return image.addBands(image.metadata('year'));
};

var ndviCollectionTime = ndviCollection.map(addTime);

// Function to calculate Kendall's tau trend
var calculateTrend = function(collection) {
  var kendallReducer = ee.Reducer.kendallsCorrelation(2);
  var kendallTrendCol = collection.select(['year', 'mean_ndvi']).reduce(kendallReducer);
  return kendallTrendCol.select('tau');
};

// Filter collections for the two periods
var p1Collection = ndviCollectionTime.filter(ee.Filter.rangeContains('year', 1984, 2004));
var p2Collection = ndviCollectionTime.filter(ee.Filter.rangeContains('year', 2005, 2024));

// Calculate trends for both periods
var p1Trend = calculateTrend(p1Collection).rename('ndviTau_1984_2004');
var p2Trend = calculateTrend(p2Collection).rename('ndviTau_2005_2024');


// Combine images into a single band image
var bandImage = p1Trend.addBands(p2Trend)

Map.addLayer(bandImage, null, 'image bands', false);
Map.addLayer(basinOutline.bounds().coveringGrid("EPSG:4326", 5e4), null, 'example grid', false)

// Function to create a covering grid and reduce points within each grid cell
function exportLargeTable(table, name) {
    var grid = table.geometry().bounds().coveringGrid("EPSG:4326", 5e4)

    var reducedPointsList = grid.toList(grid.size()).map(function(gridCell) {
        gridCell = ee.Feature(gridCell);
        var reducedPoints = bandImage.reduceRegions({
            collection: table.filterBounds(gridCell.geometry()),
            reducer: ee.Reducer.first(),
            scale: 10,
        });
        return reducedPoints;
    });

    var mergedPoints = ee.FeatureCollection(reducedPointsList).flatten();

    Export.table.toDrive({
        collection: mergedPoints,
        description: "GEE05_Tau_Extract_" + name,
        folder: " YOUR_DRIVE_FOLDER",
        fileNamePrefix: "GEE05_Tau_Extract_" + name,
    });
}

// Export the data
exportLargeTable(tables.cascadeSierraMountains, 'cascadeSierraMountains');
exportLargeTable(tables.coloradoPlateau, 'coloradoPlateau');
exportLargeTable(tables.columbianPlateau, 'columbianPlateau');
exportLargeTable(tables.greatBasin, 'greatBasin');
exportLargeTable(tables.greatPlains, 'greatPlains');
exportLargeTable(tables.northernRockies, 'northernRockies');
exportLargeTable(tables.southernRockies, 'southernRockies');
exportLargeTable(tables.wyomingBasin, 'wyomingBasin');