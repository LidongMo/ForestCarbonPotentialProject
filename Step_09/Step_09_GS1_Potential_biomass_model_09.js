//load the plot points data for each biome
var biome_1_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_1_Seed_500_FinalModel");
var biome_2_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_2_Seed_500_FinalModel");
var biome_4_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_4_Seed_508_FinalModel");
var biome_5_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_5_Seed_508_FinalModel");
var biome_6_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_6_Seed_508_FinalModel");
var biome_7_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_7_Seed_500_FinalModel");
var biome_8_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_8_Seed_500_FinalModel");
var biome_10_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_10_Seed_500_FinalModel");
var biome_11_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_11_Seed_500_FinalModel");
var biome_12_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_12_Seed_508_FinalModel");
var biome_13_Seed_500 = ee.FeatureCollection("users/leonidmoore/ForestBiomass/PlotsShapefiles/Biome_13_Seed_500_FinalModel");

    
var nitrogen = ee.Image("users/haozhima95/total_nitrogen_soilgrids_0_to_200cm").rename('Nitrogen');
var cnRatio = ee.Image("users/haozhima95/cnratio_from_caesar").rename('cnRatio');

print(nitrogen);
print(cnRatio);
// // // add those points on to the map
// Map.addLayer(Type_1001_seed_500,{},'Type 1001');
// Map.addLayer(biome_1_Seed_500,{},'Type 1002');
// Map.addLayer(Type_1003_seed_500,{},'Type 1003');
// Map.addLayer(Type_1004_seed_500,{},'Type 1004');
// print the preperties in the data frame
print(biome_1_Seed_500.limit(5));

var forestCover = ee.Image("users/leonidmoore/ForestCover_Raster_30ArcSec");
print(forestCover);
// Map.addLayer(forestCover);
var WWFBiome = ee.Image("users/devinrouth/WWF_Biomes_30ArcSec");
var conservativeWDPA = ee.Image("users/leonidmoore/Conservative_WDPA_Raster_30ArcSec").rename('WDPA').gt(0);
// Map.addLayer(WWFBiome);
// load the full composite image with all bands
// var compositeImage = ee.Image("users/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec");
// // var gapFillAndExtendBounds = require('users/devinrouth/toolbox:GapFillAndExtendBounds.js');
// // var compositeImage = gapFillAndExtendBounds.gapFillAndExtendBounds(rawCompositeImage.select(rawCompositeImage.bandNames()),rawCompositeImage.bandNames(),1000);
// // print('Filled Image', compositeImage);
// print(compositeImage.bandNames())
var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);

// filter the image by band names which selected by us
var selectedBandNames = [ 'Aridity_Index',
                          'Abs_Lat',
                          'CHELSA_Annual_Mean_Temperature',
                          'CHELSA_Annual_Precipitation',
                          'CHELSA_Isothermality',
                          'CHELSA_Max_Temperature_of_Warmest_Month',
                          'CHELSA_Mean_Diurnal_Range',
                          'CHELSA_Mean_Temperature_of_Coldest_Quarter',
                          'CHELSA_Mean_Temperature_of_Driest_Quarter',
                          'CHELSA_Mean_Temperature_of_Warmest_Quarter',
                          'CHELSA_Mean_Temperature_of_Wettest_Quarter',
                          'CHELSA_Min_Temperature_of_Coldest_Month',
                          'CHELSA_Precipitation_Seasonality',
                          'CHELSA_Precipitation_of_Coldest_Quarter',
                          'CHELSA_Precipitation_of_Driest_Month',
                          'CHELSA_Precipitation_of_Driest_Quarter',
                          'CHELSA_Precipitation_of_Warmest_Quarter',
                          'CHELSA_Precipitation_of_Wettest_Month',
                          'CHELSA_Precipitation_of_Wettest_Quarter',
                          'CHELSA_Temperature_Annual_Range',
                          'CHELSA_Temperature_Seasonality',
                          'Depth_to_Water_Table',
                          'EarthEnvCloudCover_MODCF_interannualSD',
                          'EarthEnvCloudCover_MODCF_intraannualSD',
                          'EarthEnvCloudCover_MODCF_meanannual',
                          // 'EarthEnvCloudCover_MODCF_seasonality_concentration',
                          // 'EarthEnvCloudCover_MODCF_seasonality_theta',
                          'EarthEnvTopoMed_AspectCosine',
                          'EarthEnvTopoMed_AspectSine',
                          'EarthEnvTopoMed_Eastness',
                          'EarthEnvTopoMed_Elevation',
                          'EarthEnvTopoMed_Northness',
                          'EarthEnvTopoMed_ProfileCurvature',
                          'EarthEnvTopoMed_Roughness',
                          'EarthEnvTopoMed_Slope',
                          'EarthEnvTopoMed_TangentialCurvature',
                          'EarthEnvTopoMed_TerrainRuggednessIndex',
                          'EarthEnvTopoMed_TopoPositionIndex',
                          'EarthEnvTopoMed_VectorRuggednessMeasure',
                          'SG_Absolute_depth_to_bedrock',
                          'WorldClim2_SolarRadiation_AnnualMean',
                          'WorldClim2_WindSpeed_AnnualMean',
                          'NDVI',
                          'EVI',
                          'Lai',
                          'Fpar',
                          'GlobBiomass_GrowingStockVolume',
                          'Npp',
                          'Population_Density']
// 
var finalImage = ee.Image("users/leonidmoore/ForestBiomass/20200915_Forest_Biomass_Predictors_Image")
print('Final Image Band Names',finalImage);
// load the plots to sample

// Map.addLayer(finalImage);

// Map.addLayer(Type_1001_seed_500,{},'Type 1001');
// Map.addLayer(Type_1002_seed_500,{},'Type 1002');
// Map.addLayer(Type_1003_seed_500,{},'Type 1003');
// Map.addLayer(Type_1004_seed_500,{},'Type 1004');

// Sample the points using a mapped reduceRegions
var sampledBiome_1_Seed_500 = biome_1_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 1",sampledBiome_1_Seed_500.limit(5));

var sampledBiome_2_Seed_500 = biome_2_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 2",sampledBiome_2_Seed_500.limit(5));

var sampledBiome_4_Seed_500 = biome_4_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 4",sampledBiome_4_Seed_500.limit(5));

var sampledBiome_5_Seed_500 = biome_5_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 5",sampledBiome_5_Seed_500.limit(5));

var sampledBiome_6_Seed_500 = biome_6_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 6",sampledBiome_6_Seed_500.limit(5));

var sampledBiome_7_Seed_500 = biome_7_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 7",sampledBiome_7_Seed_500.limit(5));

var sampledBiome_8_Seed_500 = biome_8_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 8",sampledBiome_8_Seed_500.limit(5));

var sampledBiome_10_Seed_500 = biome_10_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 10",sampledBiome_10_Seed_500.limit(5));
var sampledBiome_11_Seed_500 = biome_11_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 11",sampledBiome_11_Seed_500.limit(5));
var sampledBiome_12_Seed_500 = biome_12_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 12",sampledBiome_12_Seed_500.limit(5));
var sampledBiome_13_Seed_500 = biome_13_Seed_500.select(['BmssD','x','y','FrstT']).map(function(f){
  return f.set(finalImage.reduceRegion('first',f.geometry()));
});
// show the sampled results
print("Sampled Points Biome 13",sampledBiome_13_Seed_500.limit(5));

// define the predictor vector
var predictorList = ['Aridity_Index',
                     'CHELSA_Annual_Mean_Temperature',
                     'CHELSA_Annual_Precipitation',
                     'CHELSA_Isothermality',
                     'CHELSA_Max_Temperature_of_Warmest_Month',
                     'CHELSA_Mean_Diurnal_Range',
                     'CHELSA_Mean_Temperature_of_Coldest_Quarter',
                     'CHELSA_Mean_Temperature_of_Driest_Quarter',
                     'CHELSA_Mean_Temperature_of_Warmest_Quarter',
                     'CHELSA_Mean_Temperature_of_Wettest_Quarter',
                     'CHELSA_Min_Temperature_of_Coldest_Month',
                     'CHELSA_Precipitation_Seasonality',
                     'CHELSA_Precipitation_of_Coldest_Quarter',
                     'CHELSA_Precipitation_of_Driest_Month',
                     'CHELSA_Precipitation_of_Driest_Quarter',
                     'CHELSA_Precipitation_of_Warmest_Quarter',
                     'CHELSA_Precipitation_of_Wettest_Month',
                     'CHELSA_Precipitation_of_Wettest_Quarter',
                     'CHELSA_Temperature_Annual_Range',
                     'CHELSA_Temperature_Seasonality',
                     'Depth_to_Water_Table',
                     'EarthEnvTopoMed_Eastness',
                     'EarthEnvTopoMed_Elevation',
                     'EarthEnvTopoMed_Northness',
                     'EarthEnvTopoMed_ProfileCurvature',
                     'EarthEnvTopoMed_Roughness',
                     'EarthEnvTopoMed_Slope',
                     'SG_Absolute_depth_to_bedrock',
                     'WorldClim2_SolarRadiation_AnnualMean',
                     'WorldClim2_WindSpeed_AnnualMean',
                     'EarthEnvCloudCover_MODCF_interannualSD',
                     'EarthEnvCloudCover_MODCF_intraannualSD',
                     'EarthEnvCloudCover_MODCF_meanannual',
                     'EarthEnvTopoMed_AspectCosine',
                     'EarthEnvTopoMed_AspectSine',
                     'LandCoverClass_Cultivated_and_Managed_Vegetation',
                     'Human_Disturbance',
                     'LandCoverClass_Urban_Builtup',
                     'SG_Clay_Content_0_100cm',
                     'SG_Coarse_fragments_0_100cm',
                     'SG_Sand_Content_0_100cm',
                     'SG_Silt_Content_0_100cm',
                     'SG_Soil_pH_H2O_0_100cm',
                     'WDPA',
                     'cropland',
                     'grazing',
                     'pasture',
                     'rangeland']
                      
print(predictorList);

// var Type_1001_seed_500_Filtered = sampledType_1001_seed_500.filter(ee.Filter.notNull(predictorVector));
// print(Type_1001_seed_500_Filtered.size());

// var Type_1002_seed_500_Filtered = sampledType_1002_seed_500.filter(ee.Filter.notNull(predictorVector));
// print(Type_1002_seed_500_Filtered.size());

// var Type_1003_seed_500_Filtered = sampledType_1003_seed_500.filter(ee.Filter.notNull(predictorVector));
// print(Type_1003_seed_500_Filtered.size());

// var Type_1004_seed_500_Filtered = sampledType_1004_seed_500.filter(ee.Filter.notNull(predictorVector));
// print(Type_1004_seed_500_Filtered.size());


var predictorImage = finalImage.select(predictorList);

// modifiy the variables we need control 
var optimizedCultivated = predictorImage.select('LandCoverClass_Cultivated_and_Managed_Vegetation').lt(0);
var optimizedUrban = predictorImage.select('LandCoverClass_Urban_Builtup').lt(0);
var optimizedDisturbance = predictorImage.select('Human_Disturbance').lt(0);
var optimizedCropland = predictorImage.select('cropland').lt(0);
var optimizedGrazing = predictorImage.select('grazing').lt(0);
var optimizedPasture = predictorImage.select('pasture').lt(0);
var optimizedRangeland = predictorImage.select('rangeland').lt(0);
var optimizedWDPA = predictorImage.select('WDPA').gte(0);
// Map.addLayer(optimizedRangeland);
// define the vector to get the retained variable names
var retainVector = ['Aridity_Index',
                     'CHELSA_Annual_Mean_Temperature',
                     'CHELSA_Annual_Precipitation',
                     'CHELSA_Isothermality',
                     'CHELSA_Max_Temperature_of_Warmest_Month',
                     'CHELSA_Mean_Diurnal_Range',
                     'CHELSA_Mean_Temperature_of_Coldest_Quarter',
                     'CHELSA_Mean_Temperature_of_Driest_Quarter',
                     'CHELSA_Mean_Temperature_of_Warmest_Quarter',
                     'CHELSA_Mean_Temperature_of_Wettest_Quarter',
                     'CHELSA_Min_Temperature_of_Coldest_Month',
                     'CHELSA_Precipitation_Seasonality',
                     'CHELSA_Precipitation_of_Coldest_Quarter',
                     'CHELSA_Precipitation_of_Driest_Month',
                     'CHELSA_Precipitation_of_Driest_Quarter',
                     'CHELSA_Precipitation_of_Warmest_Quarter',
                     'CHELSA_Precipitation_of_Wettest_Month',
                     'CHELSA_Precipitation_of_Wettest_Quarter',
                     'CHELSA_Temperature_Annual_Range',
                     'CHELSA_Temperature_Seasonality',
                     'Depth_to_Water_Table',
                     'EarthEnvTopoMed_Eastness',
                     'EarthEnvTopoMed_Elevation',
                     'EarthEnvTopoMed_Northness',
                     'EarthEnvTopoMed_ProfileCurvature',
                     'EarthEnvTopoMed_Roughness',
                     'EarthEnvTopoMed_Slope',
                     'SG_Absolute_depth_to_bedrock',
                     'WorldClim2_SolarRadiation_AnnualMean',
                     'WorldClim2_WindSpeed_AnnualMean',
                     'EarthEnvCloudCover_MODCF_interannualSD',
                     'EarthEnvCloudCover_MODCF_intraannualSD',
                     'EarthEnvCloudCover_MODCF_meanannual',
                     'EarthEnvTopoMed_AspectCosine',
                     'EarthEnvTopoMed_AspectSine',
                     'SG_Clay_Content_0_100cm',
                     'SG_Coarse_fragments_0_100cm',
                     'SG_Sand_Content_0_100cm',
                     'SG_Silt_Content_0_100cm',
                     'SG_Soil_pH_H2O_0_100cm']
var potentialPredictorimage = finalImage.select(retainVector)
                              .addBands(optimizedCultivated)
                              .addBands(optimizedUrban)
                              .addBands(optimizedDisturbance)
                              .addBands(optimizedCropland)
                              .addBands(optimizedGrazing)
                              .addBands(optimizedPasture)
                              .addBands(optimizedRangeland)
                              .addBands(optimizedWDPA);
                              
print(potentialPredictorimage)

// define the classifier paramters
var randomForestClassifier_01 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 13, 
	minLeafPopulation:5,
	maxNodes:10,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');

var randomForestClassifier_02 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 9, 
	minLeafPopulation:1,
	maxNodes:10,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');

var randomForestClassifier_04 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 9, 
	minLeafPopulation:5,
	maxNodes:50,
	bagFraction: 0.632,
	seed: 508
}).setOutputMode('REGRESSION');

var randomForestClassifier_05 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 9, 
	minLeafPopulation:1,
	maxNodes:10,
	bagFraction: 0.632,
	seed: 508
}).setOutputMode('REGRESSION');

var randomForestClassifier_06 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 6, 
	minLeafPopulation:5,
	maxNodes:60,
	bagFraction: 0.632,
	seed: 508
}).setOutputMode('REGRESSION');

var randomForestClassifier_07 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 9, 
	minLeafPopulation:3,
	maxNodes:60,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');

var randomForestClassifier_08 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 6, 
	minLeafPopulation:10,
	maxNodes:60,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');


var randomForestClassifier_10 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 6, 
	minLeafPopulation:3,
	maxNodes:50,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');

var randomForestClassifier_11 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 9, 
	minLeafPopulation:3,
	maxNodes:30,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');

var randomForestClassifier_12 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 6, 
	minLeafPopulation:5,
	maxNodes:60,
	bagFraction: 0.632,
	seed: 508
}).setOutputMode('REGRESSION');

var randomForestClassifier_13 = ee.Classifier.smileRandomForest({
	numberOfTrees:100,
	variablesPerSplit: 9, 
	minLeafPopulation:5,
	maxNodes:60,
	bagFraction: 0.632,
	seed: 500
}).setOutputMode('REGRESSION');
// define the classifier

var trainedClassifier_01 = randomForestClassifier_01.train({
  features:sampledBiome_1_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_02 = randomForestClassifier_02.train({
  features:sampledBiome_2_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_04 = randomForestClassifier_04.train({
  features:sampledBiome_4_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_05 = randomForestClassifier_05.train({
  features:sampledBiome_5_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_06 = randomForestClassifier_06.train({
  features:sampledBiome_6_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_07 = randomForestClassifier_07.train({
  features:sampledBiome_7_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_08 = randomForestClassifier_08.train({
  features:sampledBiome_8_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_10 = randomForestClassifier_10.train({
  features:sampledBiome_10_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_11 = randomForestClassifier_11.train({
  features:sampledBiome_11_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_12 = randomForestClassifier_12.train({
  features:sampledBiome_12_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

var trainedClassifier_13 = randomForestClassifier_13.train({
  features:sampledBiome_13_Seed_500,
  classProperty:'BmssD',
  inputProperties:predictorList
});

// generate the mask for each forest type 

var mask01 = WWFBiome.eq(1);
var mask02 = WWFBiome.eq(2);
var mask04 = WWFBiome.eq(4);
var mask05 = WWFBiome.eq(5);
var mask06 = WWFBiome.eq(6);
var mask07 = WWFBiome.eq(7);
var mask08 = WWFBiome.eq(8);
var mask10 = WWFBiome.eq(10);
var mask11 = WWFBiome.eq(11);
var mask12 = WWFBiome.eq(12);
var mask13 = WWFBiome.eq(13);

// Map it for check
// Map.addLayer(mask1001);
var predictorImage01 = potentialPredictorimage.mask(mask01);
var predictorImage02 = potentialPredictorimage.mask(mask02);
var predictorImage04 = potentialPredictorimage.mask(mask04);
var predictorImage05 = potentialPredictorimage.mask(mask05);
var predictorImage06 = potentialPredictorimage.mask(mask06);
var predictorImage07 = potentialPredictorimage.mask(mask07);
var predictorImage08 = potentialPredictorimage.mask(mask08);
var predictorImage10 = potentialPredictorimage.mask(mask10);
var predictorImage11 = potentialPredictorimage.mask(mask11);
var predictorImage12 = potentialPredictorimage.mask(mask12);
var predictorImage13 = potentialPredictorimage.mask(mask13);

// Map.addLayer(predictionImage1001);
// execute the prediction
var predictedMap_01 = predictorImage01.classify(trainedClassifier_01);
var predictedMap_02 = predictorImage02.classify(trainedClassifier_02);
var predictedMap_04 = predictorImage04.classify(trainedClassifier_04);
var predictedMap_05 = predictorImage05.classify(trainedClassifier_05);
var predictedMap_06 = predictorImage06.classify(trainedClassifier_06);
var predictedMap_07 = predictorImage07.classify(trainedClassifier_07);
var predictedMap_08 = predictorImage08.classify(trainedClassifier_08);
var predictedMap_10 = predictorImage10.classify(trainedClassifier_10);
var predictedMap_11 = predictorImage11.classify(trainedClassifier_11);
var predictedMap_12 = predictorImage12.classify(trainedClassifier_12);
var predictedMap_13 = predictorImage13.classify(trainedClassifier_13);
// export to asset
Export.image.toAsset({
	image: predictedMap_01,
	description: '20200702_Fores_Biomass_500_Asset_Biome_01_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_01_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: predictedMap_02,
	description: '20200702_Fores_Biomass_500_Asset_Biome_02_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_02_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: predictedMap_04,
	description: '20200702_Fores_Biomass_500_Asset_Biome_04_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_04_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: predictedMap_05,
	description: '20200702_Fores_Biomass_500_Asset_Biome_05_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_05_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: predictedMap_06,
	description: '20200702_Fores_Biomass_500_Asset_Biome_06_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_06_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});
Export.image.toAsset({
	image: predictedMap_07,
	description: '20200702_Fores_Biomass_500_Asset_Biome_07_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_07_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});
Export.image.toAsset({
	image: predictedMap_08,
	description: '20200702_Fores_Biomass_500_Asset_Biome_08_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_08_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});
Export.image.toAsset({
	image: predictedMap_10,
	description: '20200702_Fores_Biomass_500_Asset_Biome_10_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_10_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});
Export.image.toAsset({
	image: predictedMap_11,
	description: '20200702_Fores_Biomass_500_Asset_Biome_11_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_11_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: predictedMap_12,
	description: '20200702_Fores_Biomass_500_Asset_Biome_12_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_12_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: predictedMap_13,
	description: '20200702_Fores_Biomass_500_Asset_Biome_13_Potential',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_13_500_Final_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

// Map.addLayer(image);
var forestMask = forestCover.gte(10);
// Map.addLayer(forestMask,{},'mask');
var potentialForestCover = ee.Image("users/leonidmoore/PotentialForestCover_merged_raster").rename('PotentialForestCover');
var potentialMask= potentialForestCover.gte(0.1);
// get the canopy height  layer
var biomassMap01 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_01_500_Final_Ver09');
var biomassMap02 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_02_500_Final_Ver09');
var biomassMap04 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_04_500_Final_Ver09');
var biomassMap05 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_05_500_Final_Ver09');
var biomassMap06 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_06_500_Final_Ver09');
var biomassMap07 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_07_500_Final_Ver09');
var biomassMap08 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_08_500_Final_Ver09');
var biomassMap10 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_10_500_Final_Ver09');
var biomassMap11 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_11_500_Final_Ver09');
var biomassMap12 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_12_500_Final_Ver09');
var biomassMap13 = ee.Image('users/leonidmoore/ForestBiomass/PredictedMaps/Predicted_PotentialMap_Biome_13_500_Final_Ver09');

// Map.addLayer(biomassMap1002);
var vibgYOR = ['330044', '220066', '1133cc', '33dd00', 'ffda21', 'ff6622', 'd10000'];


var allTypesComposite = ee.ImageCollection.fromImages([biomassMap01,
                                                       biomassMap02,
                                                       biomassMap04,
                                                       biomassMap05,
                                                       biomassMap06,
                                                       biomassMap07,
                                                       biomassMap08,
                                                       biomassMap10,
                                                       biomassMap11,
                                                       biomassMap12,
                                                       biomassMap13]);
// merge into one image
var mergedPredictionPotential = allTypesComposite.mosaic().exp().subtract(1);
// export the raw potential density map
Export.image.toAsset({
	image: mergedPredictionPotential,
	description: '20200702_Fores_Biomass_500_Asset_Potential_Density',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver09',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

