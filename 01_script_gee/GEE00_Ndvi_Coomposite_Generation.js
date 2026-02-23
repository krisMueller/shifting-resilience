/**
 * =================================================================================
 * SCRIPT: GEE00_NDVI_Composite_Generation.js
 * MANUSCRIPT: Shifting Resilience: Trends and Predictors of Mesic Resource Productivity
 * 
 * DESCRIPTION:
 * This script processes raw Landsat imagery (L5, L7, L8, L9) to generate annual, 
 * late-season (July 15 - Sept 30) mean NDVI composites. 
 * 
 * ROLE IN WORKFLOW:
 * This is the foundational data generation step. The output ImageCollection 
 * ("PeriodMeanNdviCollection") is the primary input used by:
 *   - GEE01 (to create the mesic mask)
 *   - GEE03 (to extract time-series NDVI values)
 * 
 * PROCESS:
 * 1. Filters Landsat collections by the study area and date range.
 * 2. Masks clouds, shadows, and snow; applies radiometric rescaling.
 * 3. Calculates NDVI for every available image.
 * 4. Reduces image collections to a single annual image using the Mean reducer.
 * 5. Scales values by 10,000 (Int16) for efficient storage.
 * 6. Exports images to a private Earth Engine Asset.
 * =================================================================================
 */

//============================================================================
//    helper functions
//============================================================================

// biome bounds
var biome = ee.FeatureCollection("projects/ee-krismueller134/assets/Misc/Sagebrush_Biome");


// region to filter to
var roi = biome.geometry().bounds();


// Rename landsat bands
var bandRenamel8l9 = function(image){
  var rename = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
  ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']); 
  return rename; 
};

// Mask clouds 
var cloudMask = function(image){
  var quality =image.select(['QA_PIXEL']); 
  // clear = no clouds, coud shadow, or snow
  var clear = quality.bitwiseAnd(1 << 1).eq(0) // dilated cloud
                .and(quality.bitwiseAnd(1 << 2).eq(0)) // cirrus
                .and(quality.bitwiseAnd(1 << 3).eq(0)) // cloud
                .and(quality.bitwiseAnd(1 << 4).eq(0)) // cloud shadow
                .and(quality.bitwiseAnd(1 << 0).eq(0)) // fill
                .and(quality.bitwiseAnd(1 << 7).eq(0)); // water
  image = image.mask(clear);
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2); // <-- apply radiometric rescaling coefficients
  
  return image.addBands(opticalBands, null, true);
};   

// function - calculate ndvi
var ndviCalc = function(image) {
  var ndvi = image.normalizedDifference(['SR_B4', 'SR_B3']);  
  ndvi = ndvi.select([0], ['ndvi']); // select ndvi band and name it ndvi
  var startTime = image.get('system:time_start'); // get image date
  var date =  ee.Date(startTime).format('YYYY-MM-DD');  // format to readable date
  var year = date.slice(0,4); // grab just the year
  year = ee.Number.parse(year); // convert year to number
  return ndvi.set('system:time_start', startTime)
    .set('date',date)
    .set('year',year);
};

//============================================================================
//    create imageCollection - filter landsat l5, l7, l8, and merge 
//                      *** Maybe need to add L4 ***
//============================================================================

var periodl5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(roi) 
    .filter((ee.Filter.date('1984-01-01', '2011-12-31')))
    .filter(ee.Filter.dayOfYear(196,274)) //Jul 15-Sep 30
    .map(cloudMask);

var periodl7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(roi) 
    .filter((ee.Filter.date('2012-01-01', '2012-12-31'))) // <------------- Maybe extend to 1999 - 2022
    .filter(ee.Filter.dayOfYear(196,274)) //Jul 15-Sep 30
    .map(cloudMask);

var periodl8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(roi)
    .filter((ee.Filter.date('2013-01-01', '2024-12-31')))
    .filter(ee.Filter.dayOfYear(196,274)) //Jul 15-Sep 30
    .map(cloudMask)
    .map(bandRenamel8l9);

var periodl9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(roi)
    .filter((ee.Filter.date('2020-01-01', '2024-12-31')))
    .filter(ee.Filter.dayOfYear(196,274)) //Jul 15-Sep 30
    .map(cloudMask)
    .map(bandRenamel8l9);

var landsatPeriod = ee.ImageCollection(periodl5.merge(periodl7).merge(periodl8).merge(periodl9)); 


//============================================================================
//    generate ndvi data by year (40 year)
//============================================================================


// ----------------------------------------- define years to export ---------------------------------------------------------------------------------
// define years of interest
var startYear = 1984;    // <------------ MAY HAVE TO CHUNK THIS UP INTO SMALLER YEAR SECTIONS
var endYear = 2024;
// --------------------------------------------------------------------------------------------------------------------------------------------------

// generate list of years
var years = ee.List.sequence(startYear, endYear);

// filter to years and calculate ndvi for all images 
//  ---------- *** ---------
var periodNdvi = landsatPeriod.filterDate(startYear.toString(), (endYear +1).toString())
  .map(ndviCalc);


// Rescale values to export as Int16 (multiply by 10,000)
// Set startTime as Jan 1 for given year.

var meanAnnualNdvi = ee.ImageCollection(
  years.map(function(year) {
  var ndvi_col = periodNdvi.filter(ee.Filter.eq('year',year)); // filter for each year
  
  // Get date of first image in the yearly series
  var firstImageMillis = ndvi_col.first().get('system:time_start'); // get the image date -------> WORK ON??
  var firstImageDate = ee.Date(firstImageMillis);
  
  // Set start time as Jan 1 for given year
  var yearString= ee.Number(year).format('%04d');
  var stringDate = ee.Date(yearString);
  var dateMillis = stringDate.millis();
  
  return ndvi_col.mean() // <----------------------- calculate mean pixels CHANGE TO MEDIAN IF NEEDED
    .multiply(10000)
    .toInt16()
    .rename('mean_ndvi') // rename the band
    .set('year',year) // set the year to property
    .set('system:time_start',dateMillis) // set the time stamp to January 1 of year.
    .set('firstImageDate',firstImageDate);
  })); 
  
  
// Map.addLayer(meanAnnualNdvi.filter(ee.Filter.eq('year', 2024)) , {}, 'mean ndvi')



//  ---------- *** ---------



// Map.addLayer(meanAnnualNdvi, {}, 'Mean')
// Map.addLayer(medianAnnualNdvi, {}, 'Median')
// Map.centerObject(roi, 11)




/////////////
// export
/////////////

var exportIndex = meanAnnualNdvi.aggregate_array("system:index").getInfo();
var exportYear = meanAnnualNdvi.aggregate_array("year").getInfo();
var exportSize = meanAnnualNdvi.size().getInfo(); 

var path2ImageCollection = 'PATH/TO/YOUR/ASSET/FOLDER';

// Use RAP CRS for export for easier comparison later
var rap = ee.ImageCollection("projects/rap-data-365417/assets/vegetation-cover-v3");
var projExport = rap.first().projection().getInfo();

for (var i = 0; i < exportSize; i++) {

  var exportImage = meanAnnualNdvi
    .filter(ee.Filter.eq('system:index',exportIndex[i]))
    .first()
    .set('description','composite NDVI - values scaled by 10000')
    .set('description2','composite NDVI - median values between DOY 196 to 274')
 
  // print(exportImage)
  
  Export.image.toAsset({
		image: exportImage,
	  description: 'ndviComposite_' + exportYear[i],
	  assetId: path2ImageCollection + '/ndviComposite_' + exportYear[i],
	  pyramidingPolicy: 'mean',
	  region: roi,
	  crs: projExport.crs,
	  crsTransform: projExport.transform,
	  maxPixels: 1e13,
  });
}