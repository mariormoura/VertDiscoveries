# VertDiscoveries
## Raw data and code base for the research approach on future discoveries of terrestrial vertebrates.

##### Raw data (spreadsheets in .csv format):
- SpeciesLevelPredictors.csv (Raw data on the species-level attributes used to predict discovery probability).
- AssemblageDataSubset.zip (Raw data on vertebrate assemblages at 220, 440, and 880 km of spatial resolution).

##### Shapefiles:
- gridcells_110km.zip (Equal area grid cell shapefiles at 110 km of spatial resolution).
- gridcells_220km.zip (Equal area grid cell shapefiles at 220 km of spatial resolution).
- gridcells_440km.zip (Equal area grid cell shapefiles at 440 km of spatial resolution).
- gridcells_880km.zip (Equal area grid cell shapefiles at 880 km of spatial resolution).
- world_limit.zip (Shapefile containing world limits).
- hatched_lines.zip (Shapefile of hatched lines over the Earth).
- landcover.zip (Shapefile containing Earth landcover, modified from https://ecoregions2017.appspot.com/).
- countries.zip (Shapefile of global administrative boundaries, modified from http://gadm.org/)

##### R-code:
- Moura_RScript1_PredictingSpeciesLevelDiscProbs.R
- Moura_RScript2_PredictingTaxonLevelDiscoveries.R
- Moura_RScript3_PredictingAssemblageLevelDiscoveries.R
- Moura_RScript4_PredictingBioregionLevelDiscoveries.R
- Moura_RScript5_PredictingCountryLevelDiscoveries.R
- Moura_RScript6_RobustnessToSamplingCondition.R
- Moura_RScript7_TimePeriodOfSpeciesDiscovery.R
- Moura_RScript8_ValidationProcedure.R
- Moura_RScript9_PlottingMainTextFigures.R

##### References:
Moura, M. R., & Jetz, W. (2021). Shortfalls and opportunities in terrestrial vertebrate species discovery. Nature Ecology and Evolution, 5:631–639.
DOI https://doi.org/10.1038/s41559-021-01411-5
