# GlenFinglas
Data cleaning, visualisation, and analysis for MSc Dissertation on Glen Finglas landscape-scale upland grazing experiment in the Southern Highlands. This repository contains no data, which can be obtained on request from the James Hutton Institute, [Copernicus Data Portal](https://ssoidp.copernicus.eu/idp/umsso20/login?faq), and [E-OBS website](https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php).

# Definitions

**block**: 6 large areas (A-F) across different habitats and 3 discontinuous locations  
**plot**: 4 (I-IV) grazing intensity experimental treatments within each **block**  
**points**: depending on survey with either 25 (triennial) or 75-83 (annual) per **plot**

# Contents

* **`AMEM` Folder**
	* `GF_climate_EOBS.R`: 
	* `GF_climate.R`: 
	* `NVCclimate.csv`:  
* **`Data Cleaning` Folder**
	* `GF_climate_EOBS.R`: 
	* `GF_climate.R`: 
	* `NVCclimate.csv`:  
* **`Models` Folder**
	* `GF_climate_EOBS.R`: 
	* `GF_climate.R`: 
	* **`Extra` Folder**
* **`QGIS` Folder**: QGIS projects documenting visualisation, workflows, and processes for data cleaning and processing. Given QGIS, they are non-replicable, but visualisation and plotting still works with the right data objects.
	* `droneGeoref_Updated.qgz`: georeferencing, merging, and visualising drone data
	* `plotCleaningMapping.qgz`: manually cleaning up spatial issues in raw shapefiles and generating Figure 2 from section 3.1 Study Area and Grazing
	* `RasterExplore.qgz`: exploring DEM data from EEA-10 and generating Figure 1 from section 3.1 Study Area and Grazing