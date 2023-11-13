# Forest carbon potential project
## This is the code for the Mo et al., 2013 in Nature.

****************************************************************************************************
This is the README text file for the forest carbon potential paper.
The code is organized into steps, each accompanied by the corresponding code.
****************************************************************************************************


### STEP 1: Binomial Correction of Species from GFBI and the Global Wood Density Database
	In this step, it is important to note that the BIOMSS R package was updated since we started our project in 2018. They changed the binomial correction approach, which resulted in modifications to the wood density data.

### STEP 2: Biome-Level Allometric Equations for Non-Tropical Biomes
	In this step, we downloaded the allometric data from the GlobAllomeTree database, then applied data cleaning and a pseudo allometric approach to generate biome-level allometric equations for the extra-tropical biomes.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
From STEP 3 to STEP 5, only a subset of the data has been provided due to the policy of the GFBI database.
Starting from STEP 6, you can access the data through Google Earth Engine.

### STEP 3: Biomass Calculation for Each Individual from GFBI Database
    Considering the file size, data use policy, and code running time, this part of the code was constructed based on the sample data from the GFBI database.

### STEP 4: Biomass Aggregation from Individual to Plot and Pixel Level.
    Considering the file size, data use policy, and code running time, this part of the code was constructed based on the sample data from the GFBI database.

### STEP 5: Biomass Data Cleaning Based on Plot Size and Tree Density
    Considering the file size, data use policy, and code running time, this part of the code was constructed based on the sample data from the GFBI database.
    After acquiring the full shapefile "GFBI_Full_MAD25_AllYear_Mean_Aggre.shp," we uploaded it to Google Earth Engine for subsequent calculations.

### STEP 6: Extract the Upper and Lower Carbon Density for Each Pixel from STEP 5
    This process was conducted through the Python API of Google Earth Engine.
    STEP 6_1: Build the upper and lower bounds based on forest cover.
    STEP 6_2: Extract carbon density using upper and lower bound information.
    STEP 6_3: Merge the data into a single file.

### STEP 7: Data Cleaning and Spatial Subsampling of Ground-Sourced Data
    STEP 7_1: Clean the data based on carbon density, create Figure S13, and export the shapefile of the cleaned data to GEE (Google Earth Engine).
    STEP 7_2: Spatially subsample the ground-sourced data for two GS model types.

### STEP 8: GS Modeling Conducted via Python API of GEE in Jupyter Notebook
    Note: For any inquiries regarding data accessibility on Google Earth Engine, please contact us ASAP.
        STEP 8_1: Preparation of the training data, grid search, and model prediction for Type 1 model based on GS data.
        STEP 8_2: Preparation of the training data, grid search, and model prediction for Type 2 model based on GS data.

    Detailed Description of the Sub-Steps:
        Firstly, we conducted a grid search through the Python API of GEE (Google Earth Engine) on a Jupyter notebook. This was to find the parameter sets that yielded the best model performance for each subsample-based model.
        Afterwards, we used these parameter sets to conduct predictions, which generated 100 ensemble biomass maps.
        Finally, we calculated the ensemble mean, variation coefficients, and 95% confidence interval maps for the following analysis.

### STEP 9: Satellite-Derived (SD) Modeling Conducted via Python API of GEE in Jupyter Notebook
    Note: For any inquiries regarding data accessibility on Google Earth Engine, please contact us.
        STEP 9_1: We used the adjusted Harmonized biomass map, ESA_CCI biomass map, and Walker et al. biomass map. For convenience, we abbreviated them as HM, SD, and WK, respectively. Here, 'SD' stands for satellite-derived, which may cause some confusion. To generate points for subsequent model training, we randomly generated 1,000,000 points on the ESA map in areas with forest cover greater than 0. Spatial subsampling for the Type 1 model utilized random subsampling, whereas the Type 2 model employed spatial subsampling to account for the uneven distribution of data.

        STEP 9_2/4/6: Preparation of Training Data, Grid Search, Model Prediction for Type 1 Model Based on SD Data
        STEP 9_3/5/7: Preparation of Training Data, Grid Search, Model Prediction for Type 2 Model Based on SD Data

        Detailed Description of the Sub-Steps:
            Firstly, we conducted a grid search through the Python API of GEE (Google Earth Engine) on a Jupyter notebook to find the parameter sets with the best model performance for each subsample-based model.
            Afterwards, we used these parameter sets to conduct predictions, which generated 100 ensemble biomass maps.
            Finally, we calculated the ensemble mean, variation coefficients, and 95% confidence interval maps for the subsequent analysis.

### STEP 10: Potential Carbon Stock Statistics for TGB, PGB, and Soil Carbon
    TGB = aboveground + belowground (only living biomass considered)
    PGB = aboveground + belowground + deadwood/litter
    STEP 10_1-10: Statistics for the potential from the two model types using data sourced from GS upper, GS lower, HM, SD, and WK.
    STEP 10_11: Statistics for soil carbon potential.

### STEP 11: Potential Carbon Stock Statistics for AGB (Above Ground Biomass) Carbon
    STEP 11_1-10: Statistics for the potential from the two model types using data sourced from GS upper, GS lower, HM, SD, and WK.

### STEP 12: Potential Carbon Stock Statistics at Country/Region Level
    In this step, we calculated the relevant carbon potential values for each country using GGoogle Earth Engine.

### STEP 13: Evaluation of the Uncertainty for Each Model at Different Carbon Stock Components (e.g., AGB, BGB, Deadwood, and Litter).

### STEP 14: Creation of Maps for Carbon Concentration and Dead Wood/Litter
    STEP 14_1: Calculation of the dead wood and litter ratio in each biome.
    STEP 14_2: Creation of the map for dead wood and litter in Google Earth Engine.
    STEP 14_3: Preparation of the Sanderman soil carbon potential map from Sanderman et al.
    STEP 14_4: Creation of the carbon concentration map, used for converting biomass into carbon stock.

### STEP 15: Spatial Autocorrelation Analysis
    STEP 15_1/2: Moran's I test for each model.
    STEP 15_3-12: Leave-one-out analysis.

### STEP 16: Statistical Analysis at Biome Level, Land Use Types, etc.

### STEP 17: Code for Creating Figures and Tables for All Main and Supplementary Figures

### STEP 18: Potential in Global Forest Plantations.

### STEP 19: Analysis of the Representativeness of the Data Based on Convex Hull Analysis
    STEP 15_1: Analysis for all covariates.
    STEP 15_2: Analysis for covariates related to human activities.
