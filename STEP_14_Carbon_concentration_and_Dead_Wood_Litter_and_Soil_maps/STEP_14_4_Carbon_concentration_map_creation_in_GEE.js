// the data of the carbon concentration was stored in the csv in the folder :Data/WoodCarbonConcentration/WoodCarbonConcentrations.csv

// load the composite
var compositeImageNew = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec");
// extract the biome layer
var biomeLayer = compositeImageNew.select("WWF_Biome");

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// define the average of the  carbon concentration in each biome
var medianQuantileScalerBiome01 = biomeLayer.eq(1).multiply(0.456).toFloat();
var medianQuantileScalerBiome02 = biomeLayer.eq(2).multiply(0.457);
var medianQuantileScalerBiome03 = biomeLayer.eq(3).multiply(0.4725);
var medianQuantileScalerBiome04 = biomeLayer.eq(4).multiply(0.465);
var medianQuantileScalerBiome05 = biomeLayer.eq(5).multiply(0.501);
var medianQuantileScalerBiome06 = biomeLayer.eq(6).multiply(0.48);
var medianQuantileScalerBiome07 = biomeLayer.eq(7).multiply(0.456);
var medianQuantileScalerBiome08 = biomeLayer.eq(8).multiply(0.465);
var medianQuantileScalerBiome09 = biomeLayer.eq(9).multiply(0.456);
var medianQuantileScalerBiome10 = biomeLayer.eq(10).multiply(0.465);
var medianQuantileScalerBiome11 = biomeLayer.eq(11).multiply(0.468);
var medianQuantileScalerBiome12 = biomeLayer.eq(12).multiply(0.4775);
var medianQuantileScalerBiome13 = biomeLayer.eq(13).multiply(0.4775);
var medianQuantileScalerBiome14 = biomeLayer.eq(14).multiply(0.456)

// merge these biome level factor maps into one
print(2)
var allBiomeCompositeMedian = ee.ImageCollection.fromImages([medianQuantileScalerBiome01,
                                                             medianQuantileScalerBiome02,
                                                             medianQuantileScalerBiome03,
                                                             medianQuantileScalerBiome04,
                                                             medianQuantileScalerBiome05,
                                                             medianQuantileScalerBiome06,
                                                             medianQuantileScalerBiome07,
                                                             medianQuantileScalerBiome08,
                                                             medianQuantileScalerBiome09,
                                                             medianQuantileScalerBiome10,
                                                             medianQuantileScalerBiome11,
                                                             medianQuantileScalerBiome12,
                                                             medianQuantileScalerBiome13,
                                                             medianQuantileScalerBiome14]);
// merge into one image
var woodCarbonConentration = allBiomeCompositeMedian.reduce(ee.Reducer.sum()).rename('WoodCarbonConcentration');
Map.addLayer(woodCarbonConentration);



var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);


Export.image.toAsset({
  image:woodCarbonConentration,
  description:'Wood_Carbon_Concentration_Map_to_Asset',
  assetId:'users/leonidmoore/ForestBiomass/Biome_level_Wood_Carbon_Conentration_Map',
  region:unboundedGeo,
  crs:'EPSG:4326',
  crsTransform: [0.008333333333333333,0,-180,0,-0.008333333333333333,90],
  maxPixels:1e13
});