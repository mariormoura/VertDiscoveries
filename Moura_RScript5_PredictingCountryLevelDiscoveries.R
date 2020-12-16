############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com


# SUPPORTING SCRIPT 5: PREDICTING COUNTRY-LEVEL DISCOVERIES
############################################################################################################################
# SUPPORTING SCRIPT 5: PREDICTING COUNTRY-LEVEL DISCOVERIES

# Steps in this script:
#  1. Prepare the files and directory structure needed for the analysis.
#  2. Get the per country percentage cover for each grid cell and spatial resolution.
#  3. Compute the discovery metrics for each country.
#  4. Get the standardized metrics of discovery potential at the country-level.

# First, clean workspace:
rm(list=ls()); gc()

# Install and load R packages needed to run the analysis:
needed_packages<-c("rgdal", "utils", "dplyr", "data.table", "filesstrings")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages, new.packages)

#####


# STEP 1 - PREPARE THE FILES AND DIRECTORY STRUCTURE NEEDED FOR THE ANALYSIS
##########################################################################################################################
# STEP 1 - PREPARE THE FILES AND DIRECTORY STRUCTURE NEEDED FOR THE ANALYSIS

# Clean the workspace and set the working directory:
rm(list=ls())
setwd("DefineYourDirectory")

# After downloading, please, define the directory where the zipped files are located:
zip_files_dir<-getwd() #  specify it mannually if different from the working directory

# Create the directory 'CountryLevelAnalysis':
dir.create("CountryLevelAnalysis", showWarnings = F)

# Create the directory 'Shapefiles' and unzip the shapefiles:
dir.create("Shapefiles", showWarnings = F)
utils::unzip(paste0(zip_files_dir, "/gridcells_110km.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/gridcells_220km.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/gridcells_440km.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/gridcells_880km.zip"), exdir="Shapefiles")

# Create the directory 'AssemblageLevelAnalysis':
dir.create("AssemblageLevelAnalysis", showWarnings = F)

# Unzip the assemblage-level estimates and move them to their respective directory:
utils::unzip(paste0(zip_files_dir, "/AssemblageLevelEstimates.zip"), exdir="AssemblageLevelAnalysis") 

# Unzip the country-level files and move them to their respective directory:
utils::unzip(paste0(zip_files_dir, "/CountryLevelEstimates.zip"), exdir="CountryLevelAnalysis") 

# The per country discovery metrics are based on the assemblage-level estimates.
# To compute the discovery metrics for each country, it is needed to known the country percentage cover for each grid cell.
# Country boundaries followed the Global Administrative Units database, version 3.6 (https://gadm.org/)

# Before running this script, there are files that need to be downloaded:
# 1. Shapefiles of equal-area grid cells at four different spatial resolution (110, 220, 440, and 880 km width).
# 2. Per country percentange cover in each 110 x 110 km grid cell.
# 3. Assemblage-level estimates of discovery probability.

#####

# STEP 2 - GET THE PER COUNTRY PERCENTAGE COVER FOR EACH GRID CELL AND EACH SPATIAL RESOLUTION
##########################################################################################################################
# STEP 2 - GET THE PER COUNTRY PERCENTAGE COVER FOR EACH GRID CELL AND EACH SPATIAL RESOLUTION

# Clean the workspace and set the working directory:
rm(list=ls())
setwd("DefineYourDirectory")

# Load the 360grid informing the proportion of land cover per country in each grid cell at 110km of spatial resolution:
country_prop<-fread("CountryLevelAnalysis/CountryLevelProportions_110km.csv", stringsAsFactors=TRUE, encoding="UTF-8")
# Informs the percentage of each grid cell that is covered by each country.

# Load the grid cells at 110 km of spatial resolution and get the identifiers at different spatial resolutions:
grid_cells<-readOGR(dsn='Shapefiles',layer='gridcells_110km')@data[,1:5]

# Merge the grid cell-level attributes to the country-level object:
country_prop$Cell_Id110<-as.character(country_prop$Cell_Id110)
grid_cells$Cell_Id110<-as.character(grid_cells$Cell_Id110)
grid_cells$Cell_Id220<-as.character(grid_cells$Cell_Id220)
grid_cells$Cell_Id440<-as.character(grid_cells$Cell_Id440)
grid_cells$Cell_Id880<-as.character(grid_cells$Cell_Id880)
country_prop_merged<-merge(country_prop, grid_cells, by='Cell_Id110', all.x=TRUE)
rm(grid_cells, country_prop)
#country_prop_merged[(country_prop_merged$Cell_Id220=="4560"),2]

# Recompute the country-level data based on the grid cell identifiers at different spatial resolutions (220, 440, 880km):
country_prop_220km<-country_prop_merged[, .(sum(PropCell, na.rm=T), .N), by = .(Cell_Id220, ISO3)]
names(country_prop_220km)<-c('Cell_Id', 'ISO3', 'PropCell', 'N_110km_cells')
country_prop_220km$PropCell<-country_prop_220km$PropCell/4

country_prop_440km<-country_prop_merged[, .(sum(PropCell, na.rm=T), .N), by = .(Cell_Id440, ISO3)] 
names(country_prop_440km)<-c('Cell_Id', 'ISO3', 'PropCell', 'N_110km_cells')
country_prop_440km$PropCell<-country_prop_440km$PropCell/16

country_prop_880km<-country_prop_merged[, .(sum(PropCell, na.rm=T), .N), by = .(Cell_Id880, ISO3)]
names(country_prop_880km)<-c('Cell_Id', 'ISO3', 'PropCell', 'N_110km_cells')
country_prop_880km$PropCell<-country_prop_880km$PropCell/64

# The total sum of UnknownSR across assemblages needs to equal the total sum of UnknownSR across countries.
# For that purpose, we need to register the proportion of land covered by each country within each grid cell. 
# The per country cover in coastal grid cells needs to be corrected by the grid cell land area instead of grid cell area.
# For inland cells, the grid cell land cover equals the grid cell area, therefore, their per country cover is already 'land-corrected'.

for(i in 1:3){ # i spatial resolution
  
  if(i==1){country_prop_cells<-country_prop_220km}
  if(i==2){country_prop_cells<-country_prop_440km}
  if(i==3){country_prop_cells<-country_prop_880km}
  
  # Create an empty dataframe to store the recomputed values:
  replaced_cells<-as.data.frame(matrix(nrow=1, ncol=4)); names(replaced_cells)<-names(country_prop_cells)
  country_prop_cells$Cell_Id<-as.factor(as.character(country_prop_cells$Cell_Id))
  
  # Recompute the per country percentange cover based on grid cell land area:
  for(j in 1:nlevels(country_prop_cells$Cell_Id)){ # j grid cells
    
    # Select all rows corresponding to a same grid cell:
    selected_cells<-country_prop_cells[country_prop_cells$Cell_Id==levels(country_prop_cells$Cell_Id)[j]]
    #selected_cells<-country_prop_cells[country_prop_cells$Cell_Id==6901]
    
    # Ignore computations for grid cells fully located in the ocean:
    selected_cells$PropCell<- (selected_cells$PropCell) / (sum(selected_cells$PropCell, na.rm=T)) # rescale
    replaced_cells<-rbind(replaced_cells, selected_cells)
    
  } # end of j for loop (across grid cells)
  
  replaced_cells<-replaced_cells[-1,]
  replaced_cells<-droplevels(replaced_cells)
  
  if(i==1){fwrite(replaced_cells, "CountryLevelAnalysis/CountryLevelProportions_220km.csv")}
  if(i==2){fwrite(replaced_cells, "CountryLevelAnalysis/CountryLevelProportions_440km.csv")}
  if(i==3){fwrite(replaced_cells, "CountryLevelAnalysis/CountryLevelProportions_880km.csv")}
  
} # end of i for loop (spatial resolutions)

#####

# STEP 3 - COMPUTE THE DISCOVERY METRICS FOR EACH COUNTRY
##########################################################################################################################
# STEP 3 - COMPUTE THE DISCOVERY METRICS FOR EACH COUNTRY

# Clean the workspace and set the working directory:
rm(list=ls())
setwd("DefineYourDirectory")

# Function to load R data using user-specified filename:
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Compute the standardization for each vertebrate class:
for(i in 1:3){ # i spatial resolutions
  
  # Import the data needed:
  if(i==1){ # 220 km of spatial resolution
    Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_220km_WithoutStandardization.RData")
    country_prop<-fread("CountryLevelAnalysis/CountryLevelProportions_220km.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  }
  
  if(i==2){ # 440 km of spatial resolution
    Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_440km_WithoutStandardization.RData")
    country_prop<-fread("CountryLevelAnalysis/CountryLevelProportions_440km.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  }
  
  if(i==3){ # 880 km of spatial resolution
    Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_880km_WithoutStandardization.RData")
    country_prop<-fread("CountryLevelAnalysis/CountryLevelProportions_880km.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  }
  
  # Filter the assemblage data to include the subsampling level with highest model performance (externally informed):
  vert_filtered<-Vert_data[Vert_data$Subsampling=="5",] # as detected in the model validation procedure
  vert_filtered<-droplevels(vert_filtered)
  
  # Merge the country proportion of each grid cell to the discovery potential results at 200 km:
  country_prop$Cell_Id<-as.character(country_prop$Cell_Id)
  vert_filtered$Cell_Id<-as.character(vert_filtered$Cell_Id)
  country_prop_merged<-merge(country_prop, vert_filtered, by='Cell_Id', all.x=T, allow.cartesian=TRUE)
  country_prop_merged<-country_prop_merged[!is.na(country_prop_merged$UnknownSR),] # remove cells without estimates of discovery potential
  
  # Compute the per country totals at 220km:
  per_country_results<-country_prop_merged[, .(.N,
                                        mean(PropUnknown*PropCell, na.rm=T), 
                                        mean(PropUnknown_L95*PropCell, na.rm=T), 
                                        mean(PropUnknown_U95*PropCell, na.rm=T),
                                        sum(UnknownSR*PropCell, na.rm=T),
                                        sum(UnknownSR_L95*PropCell, na.rm=T),
                                        sum(UnknownSR_U95*PropCell, na.rm=T)),
                                          by = .(Class, ISO3)]
  
  # Rename the columns:
  names(per_country_results)<-c("Class", "ISO3", "N_cells",
                                "PropUnknown", "PropUnknown_L95", "PropUnknown_U95", 
                              "UnknownSR", "UnknownSR_L95", "UnknownSR_U95")
  
  # Export the standardized values of 'PropKnown' and 'UnknownSR' per higher-level grouping:
  if(i==1){fwrite(per_country_results, "CountryLevelAnalysis/CountryLevelEstimates_220km_WithoutStandardization.csv")}
  if(i==2){fwrite(per_country_results, "CountryLevelAnalysis/CountryLevelEstimates_440km_WithoutStandardization.csv")}
  if(i==3){fwrite(per_country_results, "CountryLevelAnalysis/CountryLevelEstimates_880km_WithoutStandardization.csv")}
  
} # end of i for loop (spatial resolutions)

# Check if the total sums are equal for the assemblage- and country-level estimates.
# First, load the results for one spatial resolution (choose one, 220, 440, or 880 km)
Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_220km_WithoutStandardization.RData")
per_country_results<-fread("CountryLevelAnalysis/CountryLevelEstimates_220km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")
Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_440km_WithoutStandardization.RData")
per_country_results<-fread("CountryLevelAnalysis/CountryLevelEstimates_440km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")
Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_880km_WithoutStandardization.RData")
per_country_results<-fread("CountryLevelAnalysis/CountryLevelEstimates_880km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")

# Compute the total UnknownSR across assemblages for each vertebrate group:
Vert_data<-Vert_data[Vert_data$Subsampling=="5",] # as detected in the model validation procedure
Vert_data %>% dplyr::group_by(Class) %>% 
  dplyr::summarise(TotalUnknown=sum(UnknownSR, na.rm=T), 
                   TotalUnknown_L95=sum(UnknownSR_L95, na.rm=T),
                   TotalUnknown_U95=sum(UnknownSR_U95, na.rm=T))
per_country_results %>%  dplyr::group_by(Class) %>% 
  dplyr::summarise(TotalUnknown=sum(UnknownSR, na.rm=T), 
                   TotalUnknown_L95=sum(UnknownSR_L95, na.rm=T),
                   TotalUnknown_U95=sum(UnknownSR_U95, na.rm=T))

# The totals of UnknownSR should be equal, but minor deviations might occur.
# Such deviations resulted from island grid cells not accounted for among the country-level estimates.

#####

# STEP 4 - GET THE STANDARDIZED METRICS OF DISCOVERY POTENTIAL AT THE COUNTRY LEVEL
##########################################################################################################################
# STEP 4 - GET THE STANDARDIZED METRICS OF DISCOVERY POTENTIAL AT THE COUNTRY LEVEL

# Clean the workspace and set the working directory:
rm(list=ls())
setwd("DefineYourDirectory")

# The sensitivity analysis at the assemblage-level revealed that 'UnknownSR' and 'PropUnknown' are underestimated.
# Since the country-level estimates are derived from the assemblage-level ones, 
# we applied the same standardization procedure to the metrics of discovery potential.
## We divided UnknownSR by the total number of estimated discoveries (i.e., 
## sum of UnknownSR across taxa) to provide the estimated percent of total discoveries.
## PropUnknown was standardized to vary between 0 and 1, by first subtracting the minimum observed for each vertebrate class and then
## dividing by the respective range of PropUnknown. The value of 1 indicated the taxon with the highest proportion of unknown species 
## (whatever such number might be), and not necessarily a taxon with 100% of unknown species.

# Perform the standardization for each vertebrate class:
for(i in 1:3){ # i spatial resolutions
  
  # Create a empty data frame to receive the standardized values of "PropUnknown" and 'UnknownSR':
  data_names<-names(fread("CountryLevelAnalysis/CountryLevelEstimates_880km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8"))
  Std_vert_summary<-as.data.frame(matrix(nrow=1, ncol=length(data_names), NA))
  names(Std_vert_summary)<-data_names
  
  for (j in 1:5){ # j vertebrate groups
    
    # Load the assemblage data and select the metrics computed for the subsampling level of '5':
    if(i==1){Vert_data<-fread("CountryLevelAnalysis/CountryLevelEstimates_220km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")}
    if(i==2){Vert_data<-fread("CountryLevelAnalysis/CountryLevelEstimates_440km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")}
    if(i==3){Vert_data<-fread("CountryLevelAnalysis/CountryLevelEstimates_880km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")}
    
    # Filter the assemblage data to include one subsampling level and one vertebrate group:
    Vert_data$Class<-as.factor(Vert_data$Class)
    vert_filtered<-Vert_data %>% filter(Class==levels(Vert_data$Class)[j])
    vert_filtered<-droplevels(vert_filtered)
    
    # Make the 'UnknownSR' relative to the total:
    track_denominator<-sum(vert_filtered$UnknownSR, na.rm=T) # keep the same denominator
    vert_filtered$UnknownSR<-vert_filtered$UnknownSR/track_denominator
    vert_filtered$UnknownSR_L95<-vert_filtered$UnknownSR_L95/track_denominator # keep the same denominator
    vert_filtered$UnknownSR_U95<-vert_filtered$UnknownSR_U95/track_denominator # keep the same denominator
    
    # Make the 'PropUnknown' to vary between 0 and 1:
    track_denominator<-(max(vert_filtered$PropUnknown, na.rm=T) - min(vert_filtered$PropUnknown, na.rm=T)) # keep the same denominator
    track_min<-min(vert_filtered$PropUnknown, na.rm=T)
    vert_filtered$PropUnknown<-(vert_filtered$PropUnknown - track_min) / track_denominator
    vert_filtered$PropUnknown_L95<-(vert_filtered$PropUnknown_L95 - track_min) / track_denominator
    vert_filtered$PropUnknown_U95<-(vert_filtered$PropUnknown_U95 - track_min) / track_denominator
    
    # Bind
    Std_vert_summary<-rbind(Std_vert_summary, vert_filtered)
    
  } # end of j for loop (vertebrate groups)
  
  Std_vert_summary<-Std_vert_summary[-1,] # remove the first  (empty) row
  
  # The standardized PropUnknown is bounded by the maximum of 1 due to the nature of our standardization by the range:
  Std_vert_summary<-as.data.frame(Std_vert_summary)
  Std_vert_summary[Std_vert_summary$PropUnknown_U95>=1,6]<-1
  Std_vert_summary<-as.data.table(Std_vert_summary)
  
  # Export the standardized values of 'PropKnown' and 'UnknownSR' per higher-level grouping:
  if(i==1){fwrite(Std_vert_summary, "CountryLevelAnalysis/CountryLevelEstimates_220km.csv")}
  if(i==2){fwrite(Std_vert_summary, "CountryLevelAnalysis/CountryLevelEstimates_440km.csv")}
  if(i==3){fwrite(Std_vert_summary, "CountryLevelAnalysis/CountryLevelEstimates_880km.csv")}
  
  rm(Std_vert_summary)
  
} # end of i for loop (spatial resolutions)

#####