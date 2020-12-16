############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com


# SUPPORTING SCRIPT 3: PREDICTING ASSEMBLAGE-LEVEL DISCOVERIES
############################################################################################################################
# SUPPORTING SCRIPT 3: PREDICTING ASSEMBLAGE-LEVEL DISCOVERIES

# Steps in this script:
#  1. Prepare the files and directory structure needed for the analysis.
#  2. Apply the subsampling algorithm to account for range-size driven variation in representation.
#  3. Compute the discovery metrics for species assemblages at different spatial resolutions and subsampling levels.
#  4. Get the assemblage-level metrics of discovery potential for all vertebrates combined.
#  5. Get the standardized metrics of discovery potential at selected subsampling level.

# First, clean workspace:
rm(list=ls()); gc()

# Install and load R packages needed to run the analysis:
needed_packages<-c("rgdal", "utils", "data.table", "filesstrings")
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

# Before running this script, there are files that need to be downloaded:
# 1. Shapefiles of equal-area grid cells at four different spatial resolution (110, 220, 440, and 880 km width).
# 2. Spatial intersection of species geographic range and equal-area grid cell.
# 3. Species-level estimates of discovery probability.

# After downloading, please, define the directory where the zipped files are located:
zip_files_dir<-getwd() #  specify it mannually if different from the working directory

# Create the directory 'Shapefiles' and unzip the shapefiles:
dir.create("Shapefiles", showWarnings = F)
utils::unzip(paste0(zip_files_dir, "/gridcells_220km.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/gridcells_440km.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/gridcells_880km.zip"), exdir="Shapefiles")

# Create the directory 'AssemblageLevelAnalysis':
dir.create("AssemblageLevelAnalysis", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Amphibia", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Reptilia", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Mammalia", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Aves", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Amphibia/AssemblageData", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Reptilia/AssemblageData", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Mammalia/AssemblageData", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Aves/AssemblageData", showWarnings = F)

# Unzip the assemblage data and move them to their respective directory:
utils::unzip(paste0(zip_files_dir, "/AssemblageDataSubset.zip"), exdir="AssemblageLevelAnalysis")
utils::unzip(paste0(zip_files_dir, "/AssemblageLevelAnalysis/AssemblageDataAmph.zip"), 
             exdir=paste0(zip_files_dir, "/AssemblageLevelAnalysis/Amphibia/AssemblageData"))
utils::unzip(paste0(zip_files_dir, "/AssemblageLevelAnalysis/AssemblageDataRept.zip"), 
             exdir=paste0(zip_files_dir, "/AssemblageLevelAnalysis/Reptilia/AssemblageData"))
utils::unzip(paste0(zip_files_dir, "/AssemblageLevelAnalysis/AssemblageDataMamm.zip"), 
             exdir=paste0(zip_files_dir, "/AssemblageLevelAnalysis/Mammalia/AssemblageData"))
utils::unzip(paste0(zip_files_dir, "/AssemblageLevelAnalysis/AssemblageDataAves.zip"), 
             exdir=paste0(zip_files_dir, "/AssemblageLevelAnalysis/Aves/AssemblageData"))

#####

# STEP 2 - APPLY THE SUBSAMPLING ALGORITHM TO ACCOUNT FOR RANGE-SIZE DRIVEN VARIATION IN REPRESENTATION
##########################################################################################################################
# STEP 2 - APPLY THE SUBSAMPLING ALGORITHM TO ACCOUNT FOR RANGE-SIZE DRIVEN VARIATION IN REPRESENTATION
rm(list=ls())
setwd("DefineYourDirectory")

# Create list objects to indicate the level of pseudoreplication: 
spatial_scale<-c("220km", "440km", "880km")
sub_list<-list() # list of subsampling values (maximum number of subsampled occurrences per species)
sub_list[[1]]<-c(1, 5, 10, 50, 100, 200) # 220km
sub_list[[2]]<-c(1, 5, 10, 50, 100) # 440km
sub_list[[3]]<-c(1, 5, 10) # 880km

# Create directories to store the subsampled assemblages:
dir.create("AssemblageLevelAnalysis/Amphibia/AssemblageData/Subsampling", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Reptilia/AssemblageData/Subsampling", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Mammalia/AssemblageData/Subsampling", showWarnings = F)
dir.create("AssemblageLevelAnalysis/Aves/AssemblageData/Subsampling", showWarnings = F)

# List of directories and files to read:
file_list<-list()
file_list[[1]]<-"AssemblageLevelAnalysis/Amphibia/AssemblageData/AssemblageLevelData_Amph_"
file_list[[2]]<-"AssemblageLevelAnalysis/Reptilia/AssemblageData/AssemblageLevelData_Rept_"
file_list[[3]]<-"AssemblageLevelAnalysis/Mammalia/AssemblageData/AssemblageLevelData_Mamm_"
file_list[[4]]<-"AssemblageLevelAnalysis/Aves/AssemblageData/AssemblageLevelData_Aves_"
csv_scale<-c("220km.csv", "440km.csv", "880km.csv")

# List of directories and files to save:
export_list<-list()
export_list[[1]]<-"AssemblageLevelAnalysis/Amphibia/AssemblageData/Subsampling/"
export_list[[2]]<-"AssemblageLevelAnalysis/Reptilia/AssemblageData/Subsampling/"
export_list[[3]]<-"AssemblageLevelAnalysis/Mammalia/AssemblageData/Subsampling/"
export_list[[4]]<-"AssemblageLevelAnalysis/Aves/AssemblageData/Subsampling/"

# Apply the subsampling algorithm across all vertebrate groups and spatial scales:
for (i in 1:4){ # i vertebrate groups
  
  for (j in 1:3){ # j spatial resolutions
  
    # Read assemblage data
    assemblage_data<-read.csv(paste0(file_list[[i]], csv_scale[[j]]))
    
    # Load the species-level data:
    trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
    trait_data$Class<-factor(trait_data$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
    trait_data<-trait_data[trait_data$Class==levels(trait_data$Class)[i],] # filter one vertebrate group
    # trait_data<-trait_data[trait_data$YearOfDescription>=1759 & trait_data$YearOfDescription<=2014,]
    trait_data<-droplevels(trait_data) # drop unused levels
    trait_data<-trait_data[complete.cases(trait_data), ] # remove species without predictor data
  
    # Subsample only the species with trait values available:
    assemblage_data<-as.data.table(assemblage_data)
    assemblage_data<-assemblage_data[Binomial %in% as.character(trait_data$Binomial)]
    assemblage_data<-droplevels(assemblage_data)
    
    # Calculate the number of occurrences per species and add the information to the datatable:
    gr_cells_per_spp<-assemblage_data[, .(n_occ = length(unique(Cell_Id))), by = Binomial]
    assemblage_data = merge(assemblage_data, gr_cells_per_spp, by = 'Binomial', all.x = TRUE)
    
    # Apply the subsampling algorithm using different levels of pseudoreplication:
    for(k in 1:length(sub_list[[j]])){ # k = sub_list value giving the number of pseudoreplications
      ssam = list()
      assemblage_data_not_subsampled<-assemblage_data[n_occ < sub_list[[j]][k], ]
      assemblage_data_subsampled<-assemblage_data[n_occ >= sub_list[[j]][k], ]
      
      for (l in 1:100) { # j = number of iteractions
        sub = assemblage_data_subsampled[assemblage_data_subsampled[, sample(.I, sub_list[[j]][k]), by=Binomial][[2]],]
        sub[, Iter := l]
        assemblage_data_not_subsampled[, Iter := l]
        sub<-rbind(sub, assemblage_data_not_subsampled)
        ssam[[l]] = sub
        
        } # end of l for loop (across iteractions)
      
      ssam = rbindlist(ssam)
      save(ssam, file = paste0(export_list[[i]], # directory
                               paste0(format("subsampled_"), # filename
                               paste0(spatial_scale[j]), '_', # spatial scale
                               sub_list[[j]][k], # number of subsampled species occurrences 
                               ".Rdata"))) # file extension
      rm(ssam, assemblage_data_not_subsampled, assemblage_data_subsampled)
      
      } # end of k for loop (across number of pseudoreplication values)
  
    } # end of j for loop (spatial resolutions)
  
} # end of i for loop (vertebrate groups)

#####

# STEP 3 - COMPUTE THE DISCOVERY METRICS FOR SPECIES ASSEMBLAGES AT DIFFERENT SPATIAL RESOLUTIONS AND SUBSAMPLING LEVELS
##########################################################################################################################
# STEP 3 - COMPUTE THE DISCOVERY METRICS FOR SPECIES ASSEMBLAGES AT DIFFERENT SPATIAL RESOLUTIONS AND SUBSAMPLING LEVELS
rm(list=ls())
setwd("DefineYourDirectory")

# Load the species-level estimates needed to perform the taxon-level analysis:
trait_data<-fread("SpeciesLevelData.csv", stringsAsFactors=TRUE, encoding="UTF-8")
trait_data<-trait_data[, .(Binomial, DiscProb, DiscProb_L95, DiscProb_U95)]
trait_data<-trait_data[complete.cases(trait_data), ] # remove species without discovery probabilily (e.g. described in 1758)
trait_data$Class<-factor(trait_data$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# We have 4 columns in the remaining dataset, each of which is explained below:
# Binomial: Sspecies scientific name
# DiscProb: discovery probability in the year 2015
# DiscProb_L95: lower bound of the discovery probability in the year 2015
# DiscProb_U95: upper bound of the discovery probability in the year 2015

# The average discovery probability of a set of species can be read as the proportion of known species (Propknown).
# The complement of that number is the proportion of unknown species (PropUnknown).
# Based on the known species richness (KnownSR) and PropKnown, one can obtain the number of unknown species (UnknownSR).
# Please, follow the rule of three below:
### PropKnown -> KnownSR
### 100%_SR   ->   ?
### 100%_SR = (KnownSR * 100)/PropKnown

# From the relationship above, we can extract the UnknownSR:
### UnknownSR = 100%_SR - KnownSR

# The 'PropUnknown' and 'UnknownSR' are the two discovery metrics were will compute below:

# To reduce influence of outliers in averaging discovery probability values, we will use geometric mean instead of arithmetic mean.
# Create a function to compute the geometric mean for a set of species:
geometric.mean <- function(x,na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm)) }

# Create list objects to indicate the level of pseudoreplication at each spatial resolution: 
spatial_scale<-c("220km", "440km", "880km")
sub_list<-list() # list of subsampling values (maximum number of subsampled occurrences per species)
sub_list[[1]]<-c(1, 5, 10, 50, 100, 200) # 220km
sub_list[[2]]<-c(1, 5, 10, 50, 100) # 440km
sub_list[[3]]<-c(1, 5, 10) # 880km

# List of directories and files to read the observed assemblage data:
file_list<-list()
file_list[[1]]<-"AssemblageLevelAnalysis/Amphibia/AssemblageData/AssemblageLevelData_Amph_"
file_list[[2]]<-"AssemblageLevelAnalysis/Reptilia/AssemblageData/AssemblageLevelData_Rept_"
file_list[[3]]<-"AssemblageLevelAnalysis/Mammalia/AssemblageData/AssemblageLevelData_Mamm_"
file_list[[4]]<-"AssemblageLevelAnalysis/Aves/AssemblageData/AssemblageLevelData_Aves_"
csv_scale<-c("220km.csv", "440km.csv", "880km.csv")

# List of directories and files containing the subsampled assemblage data:
subsampled_list<-list()
subsampled_list[[1]]<-"AssemblageLevelAnalysis/Amphibia/AssemblageData/Subsampling/subsampled_"
subsampled_list[[2]]<-"AssemblageLevelAnalysis/Reptilia/AssemblageData/Subsampling/subsampled_"
subsampled_list[[3]]<-"AssemblageLevelAnalysis/Mammalia/AssemblageData/Subsampling/subsampled_"
subsampled_list[[4]]<-"AssemblageLevelAnalysis/Aves/AssemblageData/Subsampling/subsampled_"

# Compute the metrics of discovery potential at the assemblage level:
for (i in 1:4){ # i vertebrate groups
  
  # Create empty lists to store the assemblage-level estimates of discovery potential:
  Results_220km<-list()
  Results_440km<-list()
  Results_880km<-list()
  
  for (j in 1:3){ # j spatial resolutions
    
  # Read assemblage data
  assemblage_data<-read.csv(paste0(file_list[[i]], csv_scale[[j]]))
  assemblage_data<-as.data.table(assemblage_data)
  
  # Get the average discovery traits for observed species assemblages (without subsampling):
  assemblage_data$Cell_Id<-as.character(assemblage_data$Cell_Id)
  assemblage_data<-assemblage_data[Binomial %in% as.character(trait_data$Binomial)] # remove species without 'discovery traits'
  assemblage_data<-droplevels(assemblage_data)
  
  # Add the species-level data to the assemblage data:
  assemblage_with_traits<-merge(assemblage_data, trait_data, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
  
  # Compute the average discovery probability of species in each assemblage:
  assemblage_data<-assemblage_with_traits[, .(.N, # Species richness per assemblage
                                       (1 - geometric.mean(DiscProb, na.rm=T)), # PropUnknown per assemblage
                                       (1 - geometric.mean(DiscProb_U95, na.rm=T)), # PropUnknown_L95, note that the upper and lower bounds are now inverted
                                       (1 - geometric.mean(DiscProb_L95, na.rm=T))), # PropUnknown_U95, note that the upper and lower bounds are now inverted
                                   by = .(Cell_Id)] 
  names(assemblage_data)<-c('Cell_Id', 'Richness', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95')
  rm(assemblage_with_traits)
  
  # Identify those grid cells with at least 5 species (observed):
  cells_to_keep<-assemblage_data[which(assemblage_data$Richness>=5),1]
  
  # Remove those assemblages with < 5 species:
  assemblage_data<-assemblage_data[Cell_Id %in% as.character(cells_to_keep$Cell_Id)]
    
  # Compute the estimated UnknownSR for each assemblage:
  assemblage_data$UnknownSR<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown)) - assemblage_data$Richness
  assemblage_data$UnknownSR_L95<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown_L95)) - assemblage_data$Richness
  assemblage_data$UnknownSR_U95<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown_U95)) - assemblage_data$Richness
  
  # Register the subsampling level:
  assemblage_data$Subsampling<-"Obs"
  
  # Store the results:
  if(j==1){Results_220km[[(length(sub_list[[j]])+1)]]<-assemblage_data}
  if(j==2){Results_440km[[(length(sub_list[[j]])+1)]]<-assemblage_data}
  if(j==3){Results_880km[[(length(sub_list[[j]])+1)]]<-assemblage_data}
  rm(assemblage_data)
  
  # Get the estimates of discovery potential across all levels of the subsampled assemblages: 
  for(k in 1:length(sub_list[[j]])){ # k = sub_list value giving the number of pseudoreplications
      
    # Specify the filename of the subsampled assemblage data and load it:
    myfile<-file.path(paste0(subsampled_list[[i]], spatial_scale[[j]], "_", sub_list[[j]][k], ".Rdata"))
    load(myfile) # load the subsampled assemblage data
    
    # Add the species-level data to the assemblage data:
    assemblage_with_traits<-merge(ssam, trait_data, by='Binomial', all.x=TRUE, allow.cartesian=TRUE) # join trait data and subsampled assemblages
    
    # Get the average subsampled richness (without removing species - see next step):
    subsampled_richness<-assemblage_with_traits[, .(.N), by = .(Iter, Cell_Id)]
    names(subsampled_richness)<-c('Iter', 'Cell_Id', 'Richness')
    avg_subsampled_richness<-subsampled_richness[, .(mean(Richness, na.rm=T)), by = .(Cell_Id)]
    names(avg_subsampled_richness)<-c('Cell_Id', 'Richness')
    
    # Remove species without discovery probability from the subsampled assemblages: 
    assemblage_with_traits<-assemblage_with_traits[!is.na(assemblage_with_traits$DiscProb), ]
    
    # Compute the average discovery probability of species across all rows (.N) within each level of Iter and Cell_Id:
    subsampled_assemblages<-assemblage_with_traits[, .((1 - geometric.mean(DiscProb, na.rm=T)), # PropUnknown per assemblage
                                                      (1 - geometric.mean(DiscProb_U95, na.rm=T)), # PropUnknown_L95, note that the upper and lower bounds are now inverted
                                                      (1 - geometric.mean(DiscProb_L95, na.rm=T))), # PropUnknown_U95, note that the upper and lower bounds are now inverted
                                                  by = .(Iter, Cell_Id)]
    names(subsampled_assemblages)<-c('Iter', 'Cell_Id', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95')
    rm(assemblage_with_traits)
    
    # Get the average assemblage-level metrics across all iteractions:
    avg_subsampled_assemblages<-subsampled_assemblages[, .(mean(PropUnknown, na.rm=T), 
                                                     mean(PropUnknown_L95, na.rm=T), 
                                                     mean(PropUnknown_U95, na.rm=T)), by = .(Cell_Id)]
    names(avg_subsampled_assemblages)<-c('Cell_Id', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95')
    rm(subsampled_assemblages)
    
    # Merge the subsampled metrics of assemblages richness and discovery potential:
    avg_subsampled_assemblages<-merge(avg_subsampled_richness, avg_subsampled_assemblages, by='Cell_Id', all.y=TRUE, allow.cartesian=TRUE) # join trait data and subsampled assemblages
    
    # Remove those subsampled assemblages with les than 5 species (exact same assemblages when using the observed assemblage data):
    avg_subsampled_assemblages<-avg_subsampled_assemblages[Cell_Id %in% as.character(cells_to_keep$Cell_Id)]
    
    # Compute the estimated UnknownSR for each subsampled assemblage using rule of three:
    avg_subsampled_assemblages$UnknownSR<-(avg_subsampled_assemblages$Richness/(1-avg_subsampled_assemblages$PropUnknown)) - avg_subsampled_assemblages$Richness
    avg_subsampled_assemblages$UnknownSR_L95<-(avg_subsampled_assemblages$Richness/(1-avg_subsampled_assemblages$PropUnknown_L95)) - avg_subsampled_assemblages$Richness
    avg_subsampled_assemblages$UnknownSR_U95<-(avg_subsampled_assemblages$Richness/(1-avg_subsampled_assemblages$PropUnknown_U95)) - avg_subsampled_assemblages$Richness
    
    # Register the subsampling level:
    avg_subsampled_assemblages$Subsampling<-sub_list[[j]][k]
    
    # Store the results according to the subsampling level:
    if(j==1){Results_220km[[k]]<-avg_subsampled_assemblages}
    if(j==2){Results_440km[[k]]<-avg_subsampled_assemblages}
    if(j==3){Results_880km[[k]]<-avg_subsampled_assemblages}
    rm(avg_subsampled_assemblages)
    
    } # end of k loop (subsampling levels)
  
  rm(cells_to_keep)
  
  # Export the results for amphibians
  if(i==1 & j==1){save(Results_220km, file = paste0('AssemblageLevelAnalysis/Amphibia/AssemblageLevelEstimates_220km_WithoutStandardization.RData'))}
  if(i==1 & j==2){save(Results_440km, file = paste0('AssemblageLevelAnalysis/Amphibia/AssemblageLevelEstimates_440km_WithoutStandardization.RData'))}
  if(i==1 & j==3){save(Results_880km, file = paste0('AssemblageLevelAnalysis/Amphibia/AssemblageLevelEstimates_880km_WithoutStandardization.RData'))}
  
  # Export the results for reptiles
  if(i==2 & j==1){save(Results_220km, file = paste0('AssemblageLevelAnalysis/Reptilia/AssemblageLevelEstimates_220km_WithoutStandardization.RData'))}
  if(i==2 & j==2){save(Results_440km, file = paste0('AssemblageLevelAnalysis/Reptilia/AssemblageLevelEstimates_440km_WithoutStandardization.RData'))}
  if(i==2 & j==3){save(Results_880km, file = paste0('AssemblageLevelAnalysis/Reptilia/AssemblageLevelEstimates_880km_WithoutStandardization.RData'))}
  
  # Export the results for mammals
  if(i==3 & j==1){save(Results_220km, file = paste0('AssemblageLevelAnalysis/Mammalia/AssemblageLevelEstimates_220km_WithoutStandardization.RData'))}
  if(i==3 & j==2){save(Results_440km, file = paste0('AssemblageLevelAnalysis/Mammalia/AssemblageLevelEstimates_440km_WithoutStandardization.RData'))}
  if(i==3 & j==3){save(Results_880km, file = paste0('AssemblageLevelAnalysis/Mammalia/AssemblageLevelEstimates_880km_WithoutStandardization.RData'))}
  
  # Export the results for birds
  if(i==4 & j==1){save(Results_220km, file = paste0('AssemblageLevelAnalysis/Aves/AssemblageLevelEstimates_220km_WithoutStandardization.RData'))}
  if(i==4 & j==2){save(Results_440km, file = paste0('AssemblageLevelAnalysis/Aves/AssemblageLevelEstimates_440km_WithoutStandardization.RData'))}
  if(i==4 & j==3){save(Results_880km, file = paste0('AssemblageLevelAnalysis/Aves/AssemblageLevelEstimates_880km_WithoutStandardization.RData'))}
  
  } # end of j loop (spatial resolutions)
  
} # end of i for loop (vertebrate groups)

#####

# STEP 4 - GET THE ASSEMBLAGE-LEVEL METRICS OF DISCOVERY POTENTIAL FOR ALL VERTEBRATES COMBINED
##########################################################################################################################
# STEP 4 - GET THE ASSEMBLAGE-LEVEL METRICS OF DISCOVERY POTENTIAL FOR ALL VERTEBRATES COMBINED
rm(list=ls())
setwd("DefineYourDirectory")

# List of directories and files to read the assemblage data:
RData_scale<-c("220km_WithoutStandardization.RData", "440km_WithoutStandardization.RData", "880km_WithoutStandardization.RData")

# Function to load R data using user-specified filename:
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

for(i in 1:3){ # i spatial resolutions

  # Specify the filename of the assemblage data and load it:
  Amph_data<-loadRData(file.path(paste0("AssemblageLevelAnalysis/Amphibia/AssemblageLevelEstimates_", RData_scale[[i]])))
  Rept_data<-loadRData(file.path(paste0("AssemblageLevelAnalysis/Reptilia/AssemblageLevelEstimates_", RData_scale[[i]])))
  Mamm_data<-loadRData(file.path(paste0("AssemblageLevelAnalysis/Mammalia/AssemblageLevelEstimates_", RData_scale[[i]])))
  Aves_data<-loadRData(file.path(paste0("AssemblageLevelAnalysis/Aves/AssemblageLevelEstimates_", RData_scale[[i]])))
  
  # Bind the lists containing the assemblage data into one large data.table:
  Amph_data<-rbindlist(Amph_data); Amph_data$Class<-"Amphibia"
  Rept_data<-rbindlist(Rept_data); Rept_data$Class<-"Reptilia"
  Mamm_data<-rbindlist(Mamm_data); Mamm_data$Class<-"Mammalia"
  Aves_data<-rbindlist(Aves_data); Aves_data$Class<-"Aves"
  
  # Bind the assemblage-level results of the vertebrate groups:
  Vert_data<-rbind(Amph_data, Rept_data, Mamm_data, Aves_data)
  
  # Get the assemblage-level metrics for all terrestrial vertebrate combined: 
  Vert_data<-Vert_data[, .(sum(Richness, na.rm=T),
                            mean(PropUnknown, na.rm=T),
                            mean(PropUnknown_L95, na.rm=T),
                            mean(PropUnknown_U95, na.rm=T),
                            sum(UnknownSR, na.rm=T),
                            sum(UnknownSR_L95, na.rm=T),
                            sum(UnknownSR_U95, na.rm=T)),
                       by = .(Subsampling, Cell_Id)]
  names(Vert_data)<-c('Subsampling', 'Cell_Id', 'Richness', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95', "UnknownSR", "UnknownSR_L95", "UnknownSR_U95")
  
  # Reorder columns of the Vert_data object and bind all results into a single data.table:
  Vert_data<-Vert_data[,c(2:9,1)]
  Vert_data$Class<-"Vertebrates"
  All_data<-rbind(Amph_data, Rept_data, Mamm_data, Aves_data, Vert_data)
  rm(Amph_data, Rept_data, Mamm_data, Aves_data)
  
  if(i==1){save(All_data, file="AssemblageLevelAnalysis/AssemblageLevelEstimates_220km_WithoutStandardization.RData")}
  if(i==2){save(All_data, file="AssemblageLevelAnalysis/AssemblageLevelEstimates_440km_WithoutStandardization.RData")}
  if(i==3){save(All_data, file="AssemblageLevelAnalysis/AssemblageLevelEstimates_880km_WithoutStandardization.RData")}
  
} # end of i for loops (across spatial resolutions)

#####

# STEP 5 - GET THE STANDARDIZED METRICS OF DISCOVERY POTENTIAL AT SELECTED SUBSAMPLING LEVEL
##########################################################################################################################
# STEP 5 - GET THE STANDARDIZED METRICS OF DISCOVERY POTENTIAL AT SELECTED SUBSAMPLING LEVEL
rm(list=ls())
setwd("DefineYourDirectory")

# The sensitivity analysis at the assemblage-level pointed out that the Subsampling level of '5' showed the highest model performance.
# The sensitivity analysis also revealed that 'UnknownSR' and 'PropUnknown' are underestimated.
## To provide more informative metrics of discovery potential, we divided UnknownSR by the total number of estimated discoveries (i.e., 
## sum of UnknownSR across taxa) to provide the estimated percent of total discoveries.
## PropUnknown was standardized to vary between 0 and 1, by first subtracting the minimum observed for each vertebrate class and then
## dividing by the respective range of PropUnknown. The value of 1 indicated the taxon with the highest proportion of unknown species 
## (whatever such number might be), and not necessarily a taxon with 100% of unknown species.
  
# Function to load R data using user-specified filename:
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Perform the standardization for each vertebrate class:
for(i in 1:3){ # i spatial resolutions
  
  # Create a empty data frame to receive the standardized values of "PropUnknown" and 'UnknownSR':
  data_names<-names(loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_220km_WithoutStandardization.RData"))
  Std_vert_summary<-as.data.frame(matrix(nrow=1, ncol=length(data_names), NA))
  names(Std_vert_summary)<-data_names
  
  for (j in 1:5){ # j vertebrate groups
    
    # Load the assemblage data and select the metrics computed for the subsampling level of '5':
    if(i==1){Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_220km_WithoutStandardization.RData")}
    if(i==2){Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_440km_WithoutStandardization.RData")}
    if(i==3){Vert_data<-loadRData(file="AssemblageLevelAnalysis/AssemblageLevelEstimates_880km_WithoutStandardization.RData")}
    
    # Filter the assemblage data to include one subsampling level and one vertebrate group:
    Vert_data<-Vert_data[Vert_data$Subsampling=="5",]
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
  Std_vert_summary[Std_vert_summary$PropUnknown_U95>=1,5]<-1
  Std_vert_summary<-as.data.table(Std_vert_summary)
  
  # Export the standardized values of 'PropKnown' and 'UnknownSR' per higher-level grouping:
  if(i==1){fwrite(Std_vert_summary, "AssemblageLevelAnalysis/AssemblageLevelEstimates_220km.csv")}
  if(i==2){fwrite(Std_vert_summary, "AssemblageLevelAnalysis/AssemblageLevelEstimates_440km.csv")}
  if(i==3){fwrite(Std_vert_summary, "AssemblageLevelAnalysis/AssemblageLevelEstimates_880km.csv")}
  
  rm(Std_vert_summary)
} # end of i for loop (spatial resolutions)

#####