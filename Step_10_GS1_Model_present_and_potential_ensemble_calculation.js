// define the working boundary
var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);
// load the present  and potential forest cover 
var presentForestCover = ee.Image("users/leonidmoore/ForestCover_Raster_30ArcSec").divide(100); // uniform with potential in the  0-1 scale
var potentialForestCover = ee.Image("users/leonidmoore/PotentialForestCover_merged_raster").rename('PotentialForestCover');
// generete the forest cover masks
var presentMask= presentForestCover.gt(0.01);
var potentialMask= potentialForestCover.gt(0.01);


//  load the 15 maps for potential biomass density 
var potentialDensityVer01 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver01");
var potentialDensityVer02 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver02");
var potentialDensityVer03 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver03");
var potentialDensityVer04 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver04");
var potentialDensityVer05 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver05");
var potentialDensityVer06 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver06");
var potentialDensityVer07 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver07");
var potentialDensityVer08 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver08");
var potentialDensityVer09 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver09");
var potentialDensityVer10 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver10");
var potentialDensityVer11 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver11");
var potentialDensityVer12 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver12");
var potentialDensityVer13 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver13");
var potentialDensityVer14 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver14");
var potentialDensityVer15 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ver15");
// compile all maps into a image collection
var potentialEnsambleComposite = ee.ImageCollection.fromImages([potentialDensityVer01,
                                                                 potentialDensityVer02,
                                                                 potentialDensityVer03,
                                                                 potentialDensityVer04,
                                                                 potentialDensityVer05,
                                                                 potentialDensityVer06,
                                                                 potentialDensityVer07,
                                                                 potentialDensityVer08,
                                                                 potentialDensityVer09,
                                                                 potentialDensityVer10,
                                                                 potentialDensityVer11,
                                                                 potentialDensityVer12,
                                                                 potentialDensityVer13,
                                                                 potentialDensityVer14,
                                                                 potentialDensityVer15]);
// generate the mean and sd of the ensamble maps
var meanPotentialDensity = potentialEnsambleComposite.mean();
var sdPotentialDensity = potentialEnsambleComposite.reduce(ee.Reducer.stdDev());
var variancePotentialDensity = sdPotentialDensity.divide(meanPotentialDensity);  

var vibgYOR = ['330044', '220066', '1133cc', '33dd00', 'ffda21', 'ff6622', 'd10000'];
Map.addLayer(meanPotentialDensity.mask(potentialMask),{palette:vibgYOR,min:0, max:200},'Potential Density Mean');
Map.addLayer(sdPotentialDensity.mask(potentialMask),{palette:vibgYOR,min:0, max:20},'Potential Density SD');
Map.addLayer(variancePotentialDensity.mask(potentialMask),{palette:vibgYOR,min:0, max:0.3},'Potential Density Variance');

Export.image.toAsset({
	image: meanPotentialDensity,
	description: '20200716_Fores_Biomass_Potential_Density_mean',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ensamble_Mean',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: sdPotentialDensity,
	description: '20200716_Fores_Biomass_Potential_Density_SD',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ensamble_SD',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: variancePotentialDensity,
	description: '20200716_Fores_Biomass_Potential_Density_Variance',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Potential_Biomass_Density_Ensamble_Variance',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});




//  load the 15 maps for present biomass density 
var presentDensityVer01 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver01");
var presentDensityVer02 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver02");
var presentDensityVer03 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver03");
var presentDensityVer04 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver04");
var presentDensityVer05 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver05");
var presentDensityVer06 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver06");
var presentDensityVer07 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver07");
var presentDensityVer08 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver08");
var presentDensityVer09 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver09");
var presentDensityVer10 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver10");
var presentDensityVer11 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver11");
var presentDensityVer12 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver12");
var presentDensityVer13 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver13");
var presentDensityVer14 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver14");
var presentDensityVer15 = ee.Image("users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ver15");
// // compile all maps into a image collection
var presentEnsambleComposite = ee.ImageCollection.fromImages([presentDensityVer01,
                                                                presentDensityVer02,
                                                                presentDensityVer03,
                                                                presentDensityVer04,
                                                                presentDensityVer05,
                                                                presentDensityVer06,
                                                                presentDensityVer07,
                                                                presentDensityVer08,
                                                                presentDensityVer09,
                                                                presentDensityVer10,
                                                                presentDensityVer11,
                                                                presentDensityVer12,
                                                                presentDensityVer13,
                                                                presentDensityVer14,
                                                                presentDensityVer15]);
// // generate the mean and sd of the ensamble maps
var meanPresentDensity = presentEnsambleComposite.mean();
var sdPresentDensity = presentEnsambleComposite.reduce(ee.Reducer.stdDev());
var variancePresentDensity = sdPresentDensity.divide(meanPresentDensity);  

var convexHullResult = ee.Image("users/leonidmoore/ForestBiomass/PresentIntExt");

var vibgYOR = ['330044', '220066', '1133cc', '33dd00', 'ffda21', 'ff6622', 'd10000'];
Map.addLayer(meanPresentDensity.mask(presentMask),{palette:vibgYOR,min:0, max:200},'Present Density Mean');
Map.addLayer(sdPresentDensity.mask(presentMask),{palette:vibgYOR,min:0, max:20},'Present Density SD');
Map.addLayer(variancePresentDensity.mask(presentMask),{palette:vibgYOR,min:0, max:0.3},'Present Density Variance');
Map.addLayer(convexHullResult,{palette:vibgYOR,min:0, max:0.3},'ConvexHull');

Export.image.toAsset({
	image: meanPresentDensity,
	description: '20200716_Fores_Biomass_Present_Density_mean',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ensamble_Mean',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: sdPresentDensity,
	description: '20200716_Fores_Biomass_Present_Density_SD',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ensamble_SD',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

Export.image.toAsset({
	image: variancePresentDensity,
	description: '20200716_Fores_Biomass_Present_Density_Variance',
	assetId: 'users/leonidmoore/ForestBiomass/PredictedMaps/PredictedMap_Present_Biomass_Density_Ensamble_Variance',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});


// ######################
Export.image.toDrive({
  image:variancePotentialDensity,
  description:'Variance_Potential_Density_To_Drive',
  fileNamePrefix: 'Variance_of_Potential_Density_Export',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});

Export.image.toDrive({
  image:variancePresentDensity,
  description:'Variance_Present_Density_To_Drive',
  fileNamePrefix: 'Variance_of_Present_Density_Export',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});

var convexHullImg = ee.Image("users/leonidmoore/ForestBiomass/PresentIntExt");

Export.image.toDrive({
  image:convexHullImg,
  description:'convexHullImg_To_Drive',
  fileNamePrefix: 'convexHull_Image_Export',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});


var hansenTreeCover = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec").select('HansenEtAl_TreeCover_Year2010').divide(100)
Export.image.toDrive({
  image:hansenTreeCover,
  description:'hansenTreeCover_Export',
  fileNamePrefix: 'Hansen_Tree_Cover_Map',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});