# GlenFinglas
Data cleaning, visualisation, and analysis for MSc Dissertation on Glen Finglas landscape-scale upland grazing experiment in the Southern Highlands titled **Graze Anatomy: A Bayesian Spatiotemporal Analysis of Livestock Effects on Scottish Upland Vegetation**. This repository contains no data, which can be obtained on request from the James Hutton Institute, [Copernicus Data Portal](https://ssoidp.copernicus.eu/idp/umsso20/login?faq), and [E-OBS website](https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php). As such, file import and export locations need to be updated to match local file structures for replication.

# Definitions

**block**: 6 large areas (A-F) across different habitats and 3 discontinuous locations  
**plot**: 4 (I-IV) grazing intensity experimental treatments within each **block**  
**points**: depending on survey with either 25 (triennial) or 75-83 (annual) per **plot**

# Contents

* **`AMEM` Folder**: Calculate "approximate marginal effects at means" (AMEM) for logit-link models
	* `MarginalEffects_TurnCom.R`: Implement AMEM for beta distribution NVC community turnover models
	* `MarginalEffectsBeta_Alt.R`: *Poorly documented* alternative attempt at calculating "true" MEM for precipitation in the NVC sub-community turnover beta model. Largely based off of `MarginalEffectsBeta.R` with more involved predictions. Strange results mean that this is incomplete and needs refinement.
	* `MarginalEffectsBeta_Alt2.R`: *Poorly documented* alternative attempt at calculating "true" MEM for precipitation in the NVC sub-community turnover beta model. Based off of `MarginalEffectsBeta_Alt.R` with a less granular spatial mesh to speed up computations. Strange results mean that this is incomplete and needs refinement.
	* `MarginalEffectsBeta.R`: Implement AMEM applied for beta distribution NVC sub-community turnover models
* **`Data Cleaning` Folder**: All the work required to get raw data into usable forms for anlaysis
	* **`CEH Plant Database` Folder**: original CEH data for [vascular plant](https://catalogue.ceh.ac.uk/documents/9f097d82-7560-4ed2-af13-604a9110cf6d) and [bryophyte](https://www.brc.ac.uk/biblio/bryoatt-attributes-british-and-irish-mosses-liverworts-and-hornworts-spreadsheet) IDs and traits
	* `EOBS_DataImport.R`: Read in E-OBS data from online and save intermediate output within the study area
	* `EOBS_VarDerive.R`: Convert E-OBS raw data to our spatiotemporal scale to make usable and derive supplementary climate variables
	* `FullDatasetMerge.R`: Combine all cleaned data (NVC, Diversity, DEM, NDWI, E-OBS, and spatial objects) to create a final dataset. Derive NVC turnover that requires the unified dataset
	* `NVCSurveyClean.R`: Clean raw, annual NVC data for analysis
	* `PinFrameSurveyClean.R`: Clean species data and compute biodiversity measures
	* `SpatialProcess_DEM_NDWI.R`: Convert spatial raster data to point-based values for analysis
* **`Models` Folder**: All run models either reported or mentioned in the dissertation
	* **`Extra` Folder**: additional modelling mentioned in the write-up
		* `BryoDiv_AR_TriMesh_Inlabru.R`: Explore alternative spatial meshes
		* `BryoRichness_Inlabru.R`: Bryophyte richness model
		* `glm_glmm_gam.R`: Fit frequentist models to test suitability for height data
		* `Inlaspacetime.R`: Attempt non-separable spacetime models
		* `Turnover_Binom_Inlabru.R`: Binomial variant of NVC sub-community turnover
		* `Turnover_Com_Binom_Inlabru.R`: Binomial variant of NVC community turnover
		* `VascRichness_Inlabru.R`: Vascular plant richness model
	* `BryoDiv_Inlabru.R`: Final INLA models and checks for bryophyte diversity. This was our **main** script for refining approach and variable selection, so the most detail can be found here
	* `Height_Inlabru.R`: Final INLA models and checks for vegetation height
	* `TurnoverBeta_Inlabru.R`: Final INLA models and checks for NVC sub-community turnover
	* `TurnoverComBeta_Inlabru.R`: Final INLA models and checks for NVC community turnover
	* `VascDiv_Inlabru.R`: Final INLA models and checks for vascular plant diversity
* **`QGIS` Folder**: QGIS projects documenting visualisation, workflows, and processes for data cleaning and processing. Given QGIS, they are non-replicable, but visualisation and plotting still works with the right data objects.
	* `droneGeoref_Updated.qgz`: georeferencing, merging, and visualising drone data
	* `pointCleaning.qgz`: manually cleaning up spatial issues in raw shapefiles and generating Figure 2 from section 3.1 Study Area and Grazing
	* `RasterExplore.qgz`: exploring DEM data from EEA-10 and generating Figure 1 from section 3.1 Study Area and Grazing