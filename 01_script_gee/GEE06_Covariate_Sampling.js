/**
 * SCRIPT: GEE06_Covariate_Sampling
 * DESCRIPTION: Samples environmental covariates at 10k points for Random Forest models.
 * OUTPUT: 4 CSVs (Time periods)
 */

// Load datasets
var ecoregions = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/biome_clipped_reduced_eight_ecoregions");
var biome = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/Sagebrush_Biome")
var RAPcol = ee.ImageCollection("projects/rap-data-365417/assets/vegetation-cover-v3")
var riparianValleys_mask = ee.Image("projects/ee-krismueller134/assets/shiftingResilience_Assets/mask_riparianValleyBottoms")

// Load PDSI and Z statistics collections
var pdsiCollection = ee.ImageCollection("GRIDMET/DROUGHT").select('pdsi');
var zCollection = ee.ImageCollection("GRIDMET/DROUGHT").select('z');
var daymet = ee.ImageCollection("NASA/ORNL/DAYMET_V4")

// Load and preprocess NDVI collection
var ndviCollection = ee.ImageCollection('projects/ee-krismueller134/assets/Mesic/FinalLayers/PeriodMeanNdviCollection')
  .map(function(image) {
    return image.divide(10000).toFloat()
      .setMulti(image.toDictionary(image.propertyNames()));
  });

// Function to calculate annual, spring, or period mean
function calculateMean(collection, bandName, periodType) {
  return ee.ImageCollection(
    ee.List.sequence(1986, 2023).map(function(year) {
      var start, end;
      
      // Set date ranges based on period type
      if (periodType === 'period') {
        start = ee.Date.fromYMD(year, 7, 15);
        end = ee.Date.fromYMD(year, 9, 30);
      } else if (periodType === 'winter') {
        start = ee.Date.fromYMD(ee.Number(year).subtract(1), 12, 1);
        end = ee.Date.fromYMD(year, 2, 28);
      } else if (periodType === 'spring') {
        start = ee.Date.fromYMD(year, 3, 1);
        end = ee.Date.fromYMD(year, 6, 30);
      } else { // annual
        start = ee.Date.fromYMD(year, 1, 1);
        end = ee.Date.fromYMD(year, 12, 31);
      }
      
      return collection.filterDate(start, end)
        .select(bandName)
        .mean()
        .set('year', year);
    })
  );
}


// var annualNdvi = calculateMean(ndviCollection, 'mean_ndvi', 'annual');
var periodPdsi = calculateMean(pdsiCollection, 'pdsi', 'period');
var periodZ = calculateMean(zCollection, 'z', 'period');
var periodMaxTemp = calculateMean(daymet, 'tmax', 'period');
var winterMaxTemp = calculateMean(daymet, 'tmax', 'winter');
var periodPrecip = calculateMean(daymet, 'prcp', 'period');
var springPrecip = calculateMean(daymet, 'prcp', 'spring');





//////////////////////// public lands
// var public = ee.FeatureCollection("USGS/GAP/PAD-US/v20/easement")
//       .merge(ee.FeatureCollection("USGS/GAP/PAD-US/v20/fee"))
//       .merge(ee.FeatureCollection("USGS/GAP/PAD-US/v20/designation"))
//       .filterBounds(biome)
// var publicLands = ee.Image().paint(public, 'ID', null)
//                             .reproject('EPSG:4326', null, 30)
//                             .eq(0)
//                             .unmask()
//                             .rename('isPublicLand')


//////////////////////// Calculate frost-free days
var gridmet = ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
  .select('tmmn')
  .filterDate('1984-01-01', '2024-12-31')
  .filter(ee.Filter.dayOfYear(196, 273));

var frostFree = ee.ImageCollection(ee.List.sequence(1984, 2024).map(function(year) {
  return gridmet.filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year, 12, 31))
    .map(function(image) {
      return image.gte(273.15).rename('frost_free');
    })
    .sum()
    .set('year', year);
}));



//////////////////////// Ag Density
var irrMapper = ee.ImageCollection("UMT/Climate/IrrMapper_RF/v1_2");

function calculateAgDensity(year) {
  
  // Filter collection to year and mosaic
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);
  var yearImage = irrMapper.filterDate(startDate, endDate).mosaic();
  
  // Calculate pixel area
  var areaRaster = yearImage.eq(0).multiply(ee.Image.pixelArea())
    .unmask(0) // convert masked areas to pixels with 0 value.
  
  // calculate ag density 
  var agDensity = areaRaster.reduceNeighborhood({
    reducer: ee.Reducer.sum(),
    kernel: ee.Kernel.circle({
  		radius:500,
  	  units:"meters",
  	  normalize:false 
    }) 
  })
    
  return agDensity 
    .set('year',year)
}


// //////////////////// flood Irr Densiity
// var floodIrrMask = ee.Image("projects/ee-krismueller134/assets/Mesic/DonnellyLayers/MaskForFloodIrrigatedWetlands")
// var floodIrrDist = floodIrrMask.fastDistanceTransform({neighborhood:1024}).sqrt().multiply(ee.Image.pixelArea().sqrt()).divide(1e3)


///////////////////////////// SWE within HUC8 watersheds
// Load HUC8 watersheds
var huc8 = ee.FeatureCollection("USGS/WBD/2017/HUC08").filterBounds(biome);

function calculateHUC8SWE(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year, 12, 31);
  
  // Calculate mean SWE for the year
  var meanSWE = daymet
    .select('swe')
    .filterDate(startDate, endDate)
    .mean();
  
  // Calculate mean SWE within each HUC8 watershed
  var huc8SWE = meanSWE.reduceRegions({
    collection: huc8,
    reducer: ee.Reducer.mean(),
    scale: 4600  // Adjust scale as needed
  });
  
  // Convert the FeatureCollection to an image
  var huc8SWEImage = huc8SWE.reduceToImage({
    properties: ['mean'],
    reducer: ee.Reducer.first()
  }).rename('mean_annual_watershed_swe');
  
  return huc8SWEImage.set('year', year);
}

// Calculate HUC8 SWE for all years
var annualHUC8SWE = ee.ImageCollection(
  ee.List.sequence(1984, 2023).map(calculateHUC8SWE)
);
annualHUC8SWE = annualHUC8SWE.merge(ee.Image(0).set('year', 2024))



// ///////////////// wetland densities
var col = ee.ImageCollection("projects/ee-4932539/assets/WET/wetAg/wetWest_2015-2022")

var semiPermWetlands = col.select('wet').count()
          .updateMask(col.select('wet').count().gt(8))
          .gt(0).multiply(ee.Image.pixelArea()).unmask().reduceNeighborhood({
          reducer: ee.Reducer.sum(),
          kernel: ee.Kernel.circle(500, 'meters', false)
        });


var seasonalWetlands = col.select('wet').count()
          .updateMask(col.select('wet').count().gt(2).and(col.select('wet').count().lt(9)))
          .gt(0).multiply(ee.Image.pixelArea()).unmask().reduceNeighborhood({
          reducer: ee.Reducer.sum(),
          kernel: ee.Kernel.circle(500, 'meters', false)
        }).unmask();

var tempWetlands = col.select('wet').count()
          .updateMask(col.select('wet').count().lt(2))
          .gt(0).multiply(ee.Image.pixelArea()).unmask().reduceNeighborhood({
          reducer: ee.Reducer.sum(),
          kernel: ee.Kernel.circle(500, 'meters', false)
        }).unmask();



/////////////////// elevation 
var elevation = ee.Image("USGS/3DEP/10m")


// slope
var slope = ee.Terrain.slope(elevation)


/////////////////// Modification Layer
// var modification = ee.Image('CSP/HM/GlobalHumanModification/2016');


///////////////////// Runoffs
// Initialize map variables
var startYear = 1986;
var endYear = 2023;

// Load feature collections
var runoffFeat = ee.FeatureCollection('projects/ee-krismueller134/assets/shiftingResilience_Assets/HUC8_avg_yearly_runoff');
var runoffFeatPeriod = ee.FeatureCollection('projects/ee-krismueller134/assets/shiftingResilience_Assets/HUC8_avg_period_runoff');

// Function to create annual runoff images
var runoff = function(year) {
  year = ee.Number(year);  // Ensure year is an ee.Number
  
  return ee.Algorithms.If(
    year.eq(2023),
    // For 2023
    ee.Image.constant(0)
      .rename('annual_huc8_runoff')
      .set({
        'year': year,
        'system:time_start': ee.Date.fromYMD(year, 1, 1).millis()
      }),
    // For other years
    ee.Image().paint(
      runoffFeat.filter(ee.Filter.eq('year', year)), 
      'value'
    ).set({
      'year': year,
      'system:time_start': ee.Date.fromYMD(year, 1, 1).millis()
    }).rename('annual_huc8_runoff')
  );
};

// Function to create period runoff images
var runoffPeriod = function(year) {
  year = ee.Number(year);  // Ensure year is an ee.Number
  
  return ee.Algorithms.If(
    year.eq(2023),
    // For 2023
    ee.Image.constant(0)
      .rename('period_huc8_runoff')
      .set({
        'year': year,
        'system:time_start': ee.Date.fromYMD(year, 1, 1).millis()
      }),
    // For other years 
    ee.Image().paint(
      runoffFeatPeriod.filter(ee.Filter.eq('year', year)),
      'periodMean'
    ).set({
      'year': year,
      'system:time_start': ee.Date.fromYMD(year, 1, 1).millis()
    }).rename('period_huc8_runoff')
  );
};

// Create image collections
var runoffCollection = ee.ImageCollection(
  ee.List.sequence(startYear, endYear).map(runoff)
).filterBounds(biome);

var runoffPeriodCollection = ee.ImageCollection(
  ee.List.sequence(startYear, endYear).map(runoffPeriod)
).filterBounds(biome);
// Map.addLayer(runoffCollection.filter(ee.Filter.eq('year', 2023)))


// GLOBAL CO2 CONCENTRATIONS 
// Define your data as an array of objects
var data = [
  {year: 1986, value: 346.97},
  {year: 1987, value: 348.68},
  {year: 1988, value: 351.16},
  {year: 1989, value: 352.79},
  {year: 1990, value: 354.05},
  {year: 1991, value: 355.39},
  {year: 1992, value: 356.09},
  {year: 1993, value: 356.83},
  {year: 1994, value: 358.33},
  {year: 1995, value: 360.17},
  {year: 1996, value: 361.93},
  {year: 1997, value: 363.05},
  {year: 1998, value: 365.70},
  {year: 1999, value: 367.80},
  {year: 2000, value: 368.96},
  {year: 2001, value: 370.57},
  {year: 2002, value: 372.58},
  {year: 2003, value: 375.15},
  {year: 2004, value: 376.95},
  {year: 2005, value: 378.98},
  {year: 2006, value: 381.15},
  {year: 2007, value: 382.90},
  {year: 2008, value: 385.02},
  {year: 2009, value: 386.50},
  {year: 2010, value: 388.75},
  {year: 2011, value: 390.62},
  {year: 2012, value: 392.65},
  {year: 2013, value: 395.40},
  {year: 2014, value: 397.34},
  {year: 2015, value: 399.65},
  {year: 2016, value: 403.06},
  {year: 2017, value: 405.22},
  {year: 2018, value: 407.61},
  {year: 2019, value: 410.07},
  {year: 2020, value: 412.44},
  {year: 2021, value: 414.70},
  {year: 2022, value: 417.08},
  {year: 2023, value: 419.32} 
];

var globalCo2Collection = ee.ImageCollection(data.map(function(item) {
  // Create a constant image with the value
  return ee.Image.constant(item.value).toFloat().rename('globalCO2')
                .set('year', item.year); // Set the year as a property
})).filterBounds(biome);


///////////////////////// attempts to sample

// Function to create a single multi-band image for a given year
function createCombinedImage(year) {
  return ee.Image.cat([
    elevation,
    slope,
    ndviCollection.filter(ee.Filter.eq('year', year)).first().rename('mean_lateSeason_ndvi').toFloat(),
    calculateAgDensity(year).rename('agDensity_05km').toFloat(),
    semiPermWetlands.rename('semiPermWetlands_density_05km'),
    seasonalWetlands.rename('seasonalWetlands_density_05km'),
    tempWetlands.rename('tempWetlands_density_05km'),
    RAPcol.select('TRE').filter(ee.Filter.eq('year', year)).first().rename('percent_tree').toInt(),
    RAPcol.select('SHR').filter(ee.Filter.eq('year', year)).first().rename('percent_shrub').toInt(),
    RAPcol.select('PFG').filter(ee.Filter.eq('year', year)).first().rename('percent_perennial').toInt(),
    RAPcol.select('AFG').filter(ee.Filter.eq('year', year)).first().rename('percent_annuals').toInt(),
    periodPdsi.filter(ee.Filter.eq('year', year)).first().rename('mean_lateSeason_pdsi').toFloat(),
    periodZ.filter(ee.Filter.eq('year', year)).first().rename('mean_lateSeason_z').toFloat(),
    periodPrecip.filter(ee.Filter.eq('year', year)).first().rename('mean_lateSeason_precip').toFloat(),
    springPrecip.filter(ee.Filter.eq('year', year)).first().rename('mean_Spring_precip').toFloat(),
    winterMaxTemp.filter(ee.Filter.eq('year', year)).first().rename('mean_prevWinter_maxTemp').toFloat(),
    periodMaxTemp.filter(ee.Filter.eq('year', year)).first().rename('mean_lateSeason_maxTemp').toFloat(),
    frostFree.filter(ee.Filter.eq('year', year)).first().rename('mean_lateSeason_frostFree').toInt(),
    runoffCollection.filter(ee.Filter.eq('year', year)).first(),
    runoffPeriodCollection.filter(ee.Filter.eq('year', year)).first(),
    annualHUC8SWE.filter(ee.Filter.eq('year', year)).first().rename('mean_annual_huc8_swe').toFloat(),
    globalCo2Collection.filter(ee.Filter.eq('year', year)).first().rename('global_CO2').toFloat(),
    ee.Image.constant(year).rename('year')
  ]).set({
    'year': year,
    'date_processed': ee.Date(Date.now()).format('YYYY-MM-dd')
  });
}
Map.addLayer(createCombinedImage(2002))

var points = ee.FeatureCollection("projects/ee-krismueller134/assets/shiftingResilience_Assets/GEE06_10k_points_per_region")
Map.addLayer(points, null, 'points')
print(points.size(), 'point size')


// Function to sample points for a single year
function sampleYearPoints(year) {
  var yearImage = createCombinedImage(year);
   
  // Sample the image at the point locations
  var sampledPoints = yearImage.reduceRegions({
    collection: points,
    reducer: ee.Reducer.first(),
    scale: 30,
  });
  
  // Add year as a property to each feature
  return sampledPoints.map(function(feature) {
    return feature.set('year', year);
  });
}

// Create four lists of years
var years1 = ee.List.sequence(1986, 1995);
var years2 = ee.List.sequence(1996, 2004);
var years3 = ee.List.sequence(2005, 2014);
var years4 = ee.List.sequence(2015, 2023);

// Map over each set of years
var yearCollections1 = years1.map(function(year) {
  return sampleYearPoints(year);
});
var yearCollections2 = years2.map(function(year) {
  return sampleYearPoints(year);
});
var yearCollections3 = years3.map(function(year) {
  return sampleYearPoints(year);
});
var yearCollections4 = years4.map(function(year) {
  return sampleYearPoints(year);
});

// Merge collections for each time period
var sampledPoints1 = ee.FeatureCollection(yearCollections1).flatten();
var sampledPoints2 = ee.FeatureCollection(yearCollections2).flatten();
var sampledPoints3 = ee.FeatureCollection(yearCollections3).flatten();
var sampledPoints4 = ee.FeatureCollection(yearCollections4).flatten();

// Export the first period (1986-1995)
Export.table.toDrive({
  collection: sampledPoints1,
  description: 'GEE06_Covariates_Points_1986_1995',
  folder: 'YOUR_DRIVE_FOLDER',
  fileFormat: 'CSV'
});

// Export the second period (1996-2004)
Export.table.toDrive({
  collection: sampledPoints2,
  description: 'GEE06_Covariates_Points_1996_2004',
  folder: 'YOUR_DRIVE_FOLDER',
  fileFormat: 'CSV'
});

// Export the third period (2005-2014)
Export.table.toDrive({
  collection: sampledPoints3,
  description: 'GEE06_Covariates_Points_2005_2014',
  folder: 'YOUR_DRIVE_FOLDER',
  fileFormat: 'CSV'
});

// Export the fourth period (2015-2023)
Export.table.toDrive({
  collection: sampledPoints4,
  description: 'GEE06_Covariates_Points_2015_2023',
  folder: 'YOUR_DRIVE_FOLDER',
  fileFormat: 'CSV'
});