// load the composite
var compositeImageNew = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec");
// extract the biome layer
var biomeLayer= compositeImageNew.select("WWF_Biome").toInt();
// define the three quantile level scalers for median/50% quantile
// Tropical 1.22
// Temperate 1.33
// Boreal 1.8
// Dryland 1.18
Map.addLayer(biomeLayer);

var deadWoodLitterFactorBiome01 = ee.Image(biomeLayer).eq(1).multiply(1.22).rename("Biome01");
var deadWoodLitterFactorBiome02 = ee.Image(biomeLayer).eq(2).multiply(1.22).rename("Biome02");
var deadWoodLitterFactorBiome03 = ee.Image(biomeLayer).eq(3).multiply(1.22).rename("Biome03");
var deadWoodLitterFactorBiome04 = ee.Image(biomeLayer).eq(4).multiply(1.33).rename("Biome04");
var deadWoodLitterFactorBiome05 = ee.Image(biomeLayer).eq(5).multiply(1.33).rename("Biome05");
var deadWoodLitterFactorBiome06 = ee.Image(biomeLayer).eq(6).multiply(1.80).rename("Biome06");
var deadWoodLitterFactorBiome07 = ee.Image(biomeLayer).eq(7).multiply(1.22).rename("Biome07");
var deadWoodLitterFactorBiome08 = ee.Image(biomeLayer).eq(8).multiply(1.33).rename("Biome08");
var deadWoodLitterFactorBiome09 = ee.Image(biomeLayer).eq(9).multiply(1.22).rename("Biome09");
var deadWoodLitterFactorBiome10 = ee.Image(biomeLayer).eq(10).multiply(1.33).rename("Biome10");
var deadWoodLitterFactorBiome11 = ee.Image(biomeLayer).eq(11).multiply(1.80).rename("Biome11");
var deadWoodLitterFactorBiome12 = ee.Image(biomeLayer).eq(12).multiply(1.21).rename("Biome12");
var deadWoodLitterFactorBiome13 = ee.Image(biomeLayer).eq(13).multiply(1.21).rename("Biome13");
var deadWoodLitterFactorBiome14 = ee.Image(biomeLayer).eq(14).multiply(1.22).rename("Biome14");
var deadWoodLitterFactorBiome98 = ee.Image(biomeLayer).eq(98).multiply(1).rename("Biome98");
var deadWoodLitterFactorBiome99 = ee.Image(biomeLayer).eq(99).multiply(1).rename("Biome99");
var deadWoodLitterFactorBiome0 = ee.Image(biomeLayer).eq(0).multiply(1).rename("Biome0"); // after we did unmask, the masked regions was turned into 0, and we replaced it with 1

// merge these biome level factor maps into one
var allBiomeFactorsComposite = deadWoodLitterFactorBiome01.addBands(deadWoodLitterFactorBiome02)
                                                          .addBands(deadWoodLitterFactorBiome03)
                                                          .addBands(deadWoodLitterFactorBiome04)
                                                          .addBands(deadWoodLitterFactorBiome05)
                                                          .addBands(deadWoodLitterFactorBiome06)
                                                          .addBands(deadWoodLitterFactorBiome07)
                                                          .addBands(deadWoodLitterFactorBiome08)
                                                          .addBands(deadWoodLitterFactorBiome09)
                                                          .addBands(deadWoodLitterFactorBiome10)
                                                          .addBands(deadWoodLitterFactorBiome11)
                                                          .addBands(deadWoodLitterFactorBiome12)
                                                          .addBands(deadWoodLitterFactorBiome13)
                                                          .addBands(deadWoodLitterFactorBiome14)
                                                          .addBands(deadWoodLitterFactorBiome98)
                                                          .addBands(deadWoodLitterFactorBiome99)
                                                          .addBands(deadWoodLitterFactorBiome0);

// merge into one image
var mosaicMeanMap = allBiomeFactorsComposite.reduce(ee.Reducer.sum());
Map.addLayer(mosaicMeanMap);

var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);

Export.image.toAsset({
  image:mosaicMeanMap,
  description:'DeadWood_Litter_Factor_to_Asset',
  assetId:'users/leonidmoore/ForestBiomass/DeadWoodLitter/DeadWood_Litter_Ratio_Map',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform: [0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});


// Export.image.toCloudStorage({
//   image:mosaicMeanMap,
//   description:'DeadWood_Litter_Factor_Export_to_CloudStorage',
//   fileNamePrefix: 'DeadWoodLitterFactor/DeadWood_Litter_Biome_Level_Factor_Map',
//   region:unboundedGeo,
//   bucket:"crowtherlab_gcsb_lidong",
//   crs:'EPSG:4326',
//   crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
//   maxPixels:1e13,
//   fileFormat:'GeoTIFF'
// });

var deadWoodLitterLowerBiome01 = ee.Image(biomeLayer).eq(1).multiply(1.15).rename("Biome01");
var deadWoodLitterLowerBiome02 = ee.Image(biomeLayer).eq(2).multiply(1.15).rename("Biome02");
var deadWoodLitterLowerBiome03 = ee.Image(biomeLayer).eq(3).multiply(1.15).rename("Biome03");
var deadWoodLitterLowerBiome04 = ee.Image(biomeLayer).eq(4).multiply(1.3).rename("Biome04");
var deadWoodLitterLowerBiome05 = ee.Image(biomeLayer).eq(5).multiply(1.3).rename("Biome05");
var deadWoodLitterLowerBiome06 = ee.Image(biomeLayer).eq(6).multiply(1.68).rename("Biome06");
var deadWoodLitterLowerBiome07 = ee.Image(biomeLayer).eq(7).multiply(1.15).rename("Biome07");
var deadWoodLitterLowerBiome08 = ee.Image(biomeLayer).eq(8).multiply(1.3).rename("Biome08");
var deadWoodLitterLowerBiome09 = ee.Image(biomeLayer).eq(9).multiply(1.15).rename("Biome09");
var deadWoodLitterLowerBiome10 = ee.Image(biomeLayer).eq(10).multiply(1.3).rename("Biome10");
var deadWoodLitterLowerBiome11 = ee.Image(biomeLayer).eq(11).multiply(1.68).rename("Biome11");
var deadWoodLitterLowerBiome12 = ee.Image(biomeLayer).eq(12).multiply(1.02).rename("Biome12");
var deadWoodLitterLowerBiome13 = ee.Image(biomeLayer).eq(13).multiply(1.02).rename("Biome13");
var deadWoodLitterLowerBiome14 = ee.Image(biomeLayer).eq(14).multiply(1.15).rename("Biome14");
var deadWoodLitterLowerBiome98 = ee.Image(biomeLayer).eq(98).multiply(1).rename("Biome98");
var deadWoodLitterLowerBiome99 = ee.Image(biomeLayer).eq(99).multiply(1).rename("Biome99");
var deadWoodLitterLowerBiome0 = ee.Image(biomeLayer).eq(0).multiply(1).rename("Biome0"); // after we did unmask, the masked regions was turned into 0, and we replaced it with 1

// merge these biome level factor maps into one
var allBiomeLowerComposite = deadWoodLitterLowerBiome01.addBands(deadWoodLitterLowerBiome02)
                                                          .addBands(deadWoodLitterLowerBiome03)
                                                          .addBands(deadWoodLitterLowerBiome04)
                                                          .addBands(deadWoodLitterLowerBiome05)
                                                          .addBands(deadWoodLitterLowerBiome06)
                                                          .addBands(deadWoodLitterLowerBiome07)
                                                          .addBands(deadWoodLitterLowerBiome08)
                                                          .addBands(deadWoodLitterLowerBiome09)
                                                          .addBands(deadWoodLitterLowerBiome10)
                                                          .addBands(deadWoodLitterLowerBiome11)
                                                          .addBands(deadWoodLitterLowerBiome12)
                                                          .addBands(deadWoodLitterLowerBiome13)
                                                          .addBands(deadWoodLitterLowerBiome14)
                                                          .addBands(deadWoodLitterLowerBiome98)
                                                          .addBands(deadWoodLitterLowerBiome99)
                                                          .addBands(deadWoodLitterLowerBiome0);

// merge into one image
var mosaicLowerMap = allBiomeLowerComposite.reduce(ee.Reducer.sum());
Map.addLayer(mosaicLowerMap);

var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);

Export.image.toAsset({
  image:mosaicLowerMap,
  description:'DeadWood_Litter_Lower_to_Asset',
  assetId:'users/leonidmoore/ForestBiomass/DeadWoodLitter/DeadWood_Litter_Ratio_Lower_Map',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform: [0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});


var deadWoodLitterUpperBiome01 = ee.Image(biomeLayer).eq(1).multiply(1.3).rename("Biome01");
var deadWoodLitterUpperBiome02 = ee.Image(biomeLayer).eq(2).multiply(1.3).rename("Biome02");
var deadWoodLitterUpperBiome03 = ee.Image(biomeLayer).eq(3).multiply(1.3).rename("Biome03");
var deadWoodLitterUpperBiome04 = ee.Image(biomeLayer).eq(4).multiply(1.37).rename("Biome04");
var deadWoodLitterUpperBiome05 = ee.Image(biomeLayer).eq(5).multiply(1.37).rename("Biome05");
var deadWoodLitterUpperBiome06 = ee.Image(biomeLayer).eq(6).multiply(1.94).rename("Biome06");
var deadWoodLitterUpperBiome07 = ee.Image(biomeLayer).eq(7).multiply(1.3).rename("Biome07");
var deadWoodLitterUpperBiome08 = ee.Image(biomeLayer).eq(8).multiply(1.37).rename("Biome08");
var deadWoodLitterUpperBiome09 = ee.Image(biomeLayer).eq(9).multiply(1.3).rename("Biome09");
var deadWoodLitterUpperBiome10 = ee.Image(biomeLayer).eq(10).multiply(1.37).rename("Biome10");
var deadWoodLitterUpperBiome11 = ee.Image(biomeLayer).eq(11).multiply(1.94).rename("Biome11");
var deadWoodLitterUpperBiome12 = ee.Image(biomeLayer).eq(12).multiply(1.4).rename("Biome12");
var deadWoodLitterUpperBiome13 = ee.Image(biomeLayer).eq(13).multiply(1.4).rename("Biome13");
var deadWoodLitterUpperBiome14 = ee.Image(biomeLayer).eq(14).multiply(1.3).rename("Biome14");
var deadWoodLitterUpperBiome98 = ee.Image(biomeLayer).eq(98).multiply(1).rename("Biome98");
var deadWoodLitterUpperBiome99 = ee.Image(biomeLayer).eq(99).multiply(1).rename("Biome99");
var deadWoodLitterUpperBiome0 = ee.Image(biomeLayer).eq(0).multiply(1).rename("Biome0"); // after we did unmask, the masked regions was turned into 0, and we replaced it with 1

// merge these biome level factor maps into one
var allBiomeUpperComposite = deadWoodLitterUpperBiome01.addBands(deadWoodLitterUpperBiome02)
                                                          .addBands(deadWoodLitterUpperBiome03)
                                                          .addBands(deadWoodLitterUpperBiome04)
                                                          .addBands(deadWoodLitterUpperBiome05)
                                                          .addBands(deadWoodLitterUpperBiome06)
                                                          .addBands(deadWoodLitterUpperBiome07)
                                                          .addBands(deadWoodLitterUpperBiome08)
                                                          .addBands(deadWoodLitterUpperBiome09)
                                                          .addBands(deadWoodLitterUpperBiome10)
                                                          .addBands(deadWoodLitterUpperBiome11)
                                                          .addBands(deadWoodLitterUpperBiome12)
                                                          .addBands(deadWoodLitterUpperBiome13)
                                                          .addBands(deadWoodLitterUpperBiome14)
                                                          .addBands(deadWoodLitterUpperBiome98)
                                                          .addBands(deadWoodLitterUpperBiome99)
                                                          .addBands(deadWoodLitterUpperBiome0);

// merge into one image
var mosaicUpperMap = allBiomeUpperComposite.reduce(ee.Reducer.sum());
Map.addLayer(mosaicLowerMap);

var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);

Export.image.toAsset({
  image:mosaicUpperMap,
  description:'DeadWood_Litter_Upper_to_Asset',
  assetId:'users/leonidmoore/ForestBiomass/DeadWoodLitter/DeadWood_Litter_Ratio_Upper_Map',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform: [0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});