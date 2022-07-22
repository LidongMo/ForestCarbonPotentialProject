### ForestCarbonPotentialProject

This is the readMe text file for the forest biomass paper
Code will be arranged by steps with corresponding code.


# STEP 1: Binomial correction of the species from GFBI and the global wood density database
     In this step we need to mention that the BIOMSS R package was updated since we started our project in 2018. The changed the binomial correction approach which result in the wood density of 
# STEP 2: Biome level allometric equation for the non-tropical biomes
    In this step we downloaded the allometric data from the GlobAllomeTree database, then applied the data cleaning and pseudo allometric approach to generate the biome level allometric equations for the extra-tropical biomes.

STEP 3: Biomass calculation for each individual from GFBi database
    Considering the file size, data use policy and code running time, in this part of the code was constructed based on the sample data from GFBi database.

STEP 4: Biomass aggregation from individual to plot and pixel level
    Considering the file size, data use policy and code running time, in this part of the code was constructed based on the sample data from GFBi database.

STEP 5: Biomass data cleaning based on plot size and tree density
    Considering the file size, data use policy and code running time, in this part of the code was constructed based on the sample data from GFBi database.
    After we acquired the full shapefile "GFBI_Full_MAD25_AllYear_Mean_Aggre.shp", we uploaded it to Google earth engine for the following calculation.

STEP 6: Extract covariates for all the points from STEP 5
This was conducted through the Python API of Google earth engine.


STEP 7: GS1 Model grid search (potential and present biomass)
    For the GS1 model, which was conducted in the early stage of this project, we did a 15 model ensemble approach for each biome.
    In each step we also generate the shapefiles for each biome, which were stored in the folder 'ShapeFilesBiome'.
    After the grid search finished, the parameters were manully input into the Goolge earthn engine UI on webpage.

STEP 8: GS1 Model mapping for present
    This part was coded with the language JavaScript on the Goolge earthn engine UI on webpage.
 
STEP 9: GS1 Model mapping for potential
    This part was coded with the language JavaScript on the Goolge earthn engine UI on webpage.  

STEP 10: GS1 Model present and potential maps ensemble calculation
    This part was coded with the language JavaScript on the Goolge earthn engine UI on webpage.

STEP 11: GS2 Model grid search and mapping
    This was conducted by Python API of GEE in jupyter notebook.
    First we did grid subsampling of the GS data points inside the protective area, with 100 repeatition. All those subsamples were uploaded to GEE for futher grid search and mapping, and ensemble calculation
    Second, we conducted all these grid search, maping and ensemble calculation on Jupyter notebook through Python API of GEE

STEP 12: SD1 Model grid search and potential Mapping

STEP 13: SD2 Model grid search and potential Mapping

STEP 15: Maps exportation 

STEP 16: Statistics for all the biomass products

STEP 17: Figures making

