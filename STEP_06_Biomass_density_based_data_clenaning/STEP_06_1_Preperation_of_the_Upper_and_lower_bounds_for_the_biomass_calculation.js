!!!!!!!!!!!!!!!!!
// any issue about the acces of the data please contact Lidong. He will give you the access to those files on google earth engine.

// load the 30m resolution Hansen tree cover map
var forestCover30m = ee.Image("users/crowtherlab/Composite/HansenEtAl_TreeCover2010");
// define the boundary for the maps 
var unboundedGeo = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);

// print the image information
print(forestCover30m)
// visualize the map
Map.addLayer(forestCover30m)
// build a mask to mask out the pixels with tree cover smaller than 10 percent
var treeCovermask30m = forestCover30m.gte(10);
// apply the mask on the raw tree cover 30m 
var maskedForestCover30m = forestCover30m.mask(treeCovermask30m);
Map.addLayer(maskedForestCover30m)
//reduce resolution of the map from 30m to 1km by averaging the non null pixels
// define the projection for the 1km . load an existing map as standard projection
var compositeImage = ee.Image("projects/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec").select('HansenEtAl_TreeCover_Year2010');
var stdProj = compositeImage.projection();
// print(compositeImage.bandNames())
// define the cover with percentage larger than 10percent as true forest
var treeCoverMean1km = maskedForestCover30m
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 10000
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: stdProj
    });

// build a mask to mask out the pixels with tree cover larger than 0 percent
var treeCovermask30mAll = forestCover30m.gt(0);
// apply the mask on the raw tree cover 30m 
var maskedForestCover30mAll = forestCover30m.mask(treeCovermask30mAll);

var treeCoverMax1km = maskedForestCover30mAll
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.max(),
      maxPixels: 10000
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: stdProj
    });

// aggreate the cover map into 1km resolution with all the 30m pixels larger than 0 percent
var treeCoverStdMean1km = forestCover30m
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 10000
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: stdProj
    });
    
// export the map of the lower canopy cover scaler
Export.image.toAsset({
	image: treeCoverStdMean1km.divide(treeCoverMax1km),
	description: '20221004_Forest_Cover_scaler_Max_to_Assest',
	assetId: 'users/leonidmoore/ForestBiomass/ForestCover/Forest_cover_MaxScaler',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});

// calculate the scaler to Assests of the upper canopy cover scaler
var meanScalerImage = treeCoverStdMean1km.divide(treeCoverMean1km);
Export.image.toAsset({
	image: meanScalerImage,
	description: '20221004_Forest_Cover_scaler_Mean_to_Assest',
	assetId: 'users/leonidmoore/ForestBiomass/ForestCover/Forest_cover_MeanScaler',
	region: unboundedGeo,
  crs:'EPSG:4326',
  crsTransform:[0.008333333333333333,0,-180,0,-0.008333333333333333,90],
	maxPixels: 1e13
});