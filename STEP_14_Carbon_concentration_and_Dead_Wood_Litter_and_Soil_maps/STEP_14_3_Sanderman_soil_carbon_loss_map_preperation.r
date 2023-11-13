library(raster)
library(snowfall)
# paper is from Sanderman 2017. Soil carbon debt of 12,000 years of human land use
# https://www.pnas.org/content/114/36/9575

# STEP 1. generate the soil carbon stock difference map
#####################################################################
# load the present and no land using soil carbon stock density map for 0-200cm
presentSoilCarbon = raster("Data/SoilCarbonModel/SandermanSoilCarbonLossMap/SOCS_0_200cm_year_2010AD_10km.tif")
noLandUseSoilCarbon = raster("Data/SoilCarbonModel/SandermanSoilCarbonLossMap/SOCS_0_200cm_year_NoLU_10km.tif")

# get the difference by subtracting the present with the no land use carbon map
carbonStockLoss = noLandUseSoilCarbon - presentSoilCarbon
carbonStockLoss[carbonStockLoss<0] =0
# disaggregate from 1km to 10km
disaggregatedCarbonStockLoss = disaggregate(carbonStockLoss,fact=10)
# and write to local folder
writeRaster(disaggregatedCarbonStockLoss,"Data/SoilCarbonModel/SandermanSoilCarbonLossMap/SOCS_0_200cm_year_Diff_1km_Present_subtract_NoLU.tif",overwrite=T)
# afterwards, this map was uploaded to GEE for the number calculation 

disaggregatedCarbonStockPresent = disaggregate(presentSoilCarbon,fact=10)
writeRaster(disaggregatedCarbonStockPresent,"Data/SoilCarbonModel/SandermanSoilCarbonLossMap/SOCS_0_200cm_1km_Present.tif",overwrite=T)


# STEP 2. generate the soil carbon stock absolute error map based on Sandermann's code
#####################################################################
## Derive cumulative SOCS for 0-2 m ----
sum_SOCS = function(tifs, depthT = c(30,70,100), year, depth.sel=200){
  out.tif = paste0("Data/SoilCarbonModel/Maps/SOCS_0_", c(depthT[1], sum(depthT[1:2]), sum(depthT[1:3])), "cm_year_",year,"_10km_Obs_Error.tif")
  s = stack(tifs)
  s = as(s, "SpatialGridDataFrame")
  for(i in 1:ncol(s)){ s@data[,i] = ifelse(s@data[,i]<0, 0, s@data[,i]) }
  x = list(NULL)
  for(i in 1:(length(tifs)-1)){
    x[[i]] = rowMeans(s@data[,i:(i+1)], na.rm=TRUE)*depthT[i]/100 
  }
  for(k in 1:length(depthT)){
    if(depthT[k]==100){
      s$SOCS = rowSums(as.data.frame(x[1:3]), na.rm=TRUE)
    }
    if(depthT[k]==70){
      s$SOCS = rowSums(as.data.frame(x[1:2]), na.rm=TRUE)
    }
    if(depthT[k]==30){
      s$SOCS = rowSums(as.data.frame(x[1]), na.rm=TRUE)
    }
    ## tones / ha
    writeGDAL(s["SOCS"], out.tif[3], type="Int16", mvFlag=-32767, options="COMPRESS=DEFLATE")
  }
}
periods = c("2010AD")
year = "2010AD"
tif.lst <- lapply(periods, function(x){paste0("Data/SoilCarbonModel/Maps/abs_error_OCDENS_",c(0,30,100,200),"cm_",year,"_10km_ll.tif")})
# source("SoilCarbonModel/WHRC_functions.r") we didn't use the source function directly, since we did adaptation on the file name and path
sfInit(parallel=TRUE, cpus=length(periods))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("sum_SOCS", "tif.lst", "periods")
missing.lst <- sfLapply(1:length(periods), function(i){sum_SOCS(tif.lst[[i]], year=periods[i])})
sfStop()

## for NoLU and 2010:
sum_SOCS(tif.lst[[1]], year=periods[1], depth.sel=c(30,100,200))

# read the worte abs error map
presentSoilCarbon = raster("Data/SoilCarbonModel/SandermanSoilCarbonLossMap/SOCS_0_200cm_year_2010AD_10km.tif")
absErrorMap = raster("Data/SoilCarbonModel/Maps/SOCS_0_200cm_year_2010AD_10km_Obs_Error.tif")
noLandUseSoilCarbon = raster("Data/SoilCarbonModel/SandermanSoilCarbonLossMap/SOCS_0_200cm_year_NoLU_10km.tif")

relativeErrorMap = absErrorMap/presentSoilCarbon
# replace all the pixels with the relative error lager than 1 by 1
relativeErrorMap[relativeErrorMap>1] =1
# generate the lower and upper boundary of the present soil carbon, use the  noland use to subtract the present SOCS
lowerBoundMap = presentSoilCarbon-(presentSoilCarbon*relativeErrorMap)
upperBoundMap = presentSoilCarbon+(presentSoilCarbon*relativeErrorMap)
# generate  
lowerCarbonStockLoss = noLandUseSoilCarbon - upperBoundMap
upperCarbonStockLoss = noLandUseSoilCarbon - lowerBoundMap
# replace the negative values by zero
lowerCarbonStockLoss[lowerCarbonStockLoss<0] = 0
upperCarbonStockLoss[upperCarbonStockLoss<0] = 0
# disaggregate from 1km to 10km
disaggregatedLowerCarbonStockLoss = disaggregate(lowerCarbonStockLoss,fact=10)
disaggregatedUpperCarbonStockLoss = disaggregate(upperCarbonStockLoss,fact=10)
# write the two boundary layers into local folder 
writeRaster(disaggregatedLowerCarbonStockLoss,"Data/SoilCarbonModel/maps/SOCS_0_200cm_Diff_1km_Present_subtract_NoLU_Lower.tif",overwrite=T)
writeRaster(disaggregatedUpperCarbonStockLoss,"Data/SoilCarbonModel/maps/SOCS_0_200cm_Diff_1km_Present_subtract_NoLU_Upper.tif",overwrite=T)


