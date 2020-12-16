############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com

# SUPPORTING SCRIPT 1: PREDICTING SPECIES-LEVEL DISCOVERY PROBABILITY
############################################################################################################################
# SUPPORTING SCRIPT 1: PREDICTING SPECIES-LEVEL DISCOVERY PROBABILITY

# Steps in this script:
#  1. Load, understand and prepare the dataset for analysis.
#  2. Check different family error distribution and create model formulas with all possible predictor combinations.
#  3. Run the model averaging procedure using all possible accelerated failure time models.
#  4. Get the average weighted coefficients using the outputs from the model averaging procedure.
#  5. Get species-level predictions of discovery probability and estimated year of discovery for all AFT models.
#  6. Get the average weighted discovery probability and estimated year of discovery.

# First, clean workspace:
rm(list=ls()); gc()
setwd("DefineYourDirectory")

# Install and load R packages needed to run the analysis:
needed_packages<-c("survival", "flexsurv", "MuMIn", "usdm", "plyr", "data.table", "foreach", "doParallel")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages,new.packages)

#####

# STEP 1 - LOAD, UNDERSTAND, AND PREPARE THE DATASET FOR ANALYSIS
##########################################################################################################################
# STEP 1 - LOAD, UNDERSTAND, AND PREPARE THE DATASET FOR ANALYSIS

# Clean the workspace and set the working directory:
rm(list=ls())
setwd("DefineYourDirectory")

# Load the response and predictor variables:
trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")

# We have 16 columns in the dataset, each of which is explained below:
# Binomial: species scientific name
# Class: taxonomic class
# Order: taxonomic order
# Family: taxonomic family
# YearOfDescription: year in which the species was formally described for the first time.
# BodySize: maximum body size (body length, for amphibians; body mass for reptiles, mammals, and birds)
# RangeSize: geographic range size as the number of occurrences in a 110 x 110 km equal area grid cell scheme
# RangeRarity: average within-range endemism richness known at the year in which the species was discovered
# AnnuPrecip: average within-range annual precipitation
# AnnuMeanTemp: average within-range annual mean temperature
# TempSeasonality: average within-range temperature seasonality (standard deviation ×100)
# PrecipSeasonality: average within-range precipitation seasonality (coefficient of variation)
# Elevation: average within-range elevation
# HumanDensity: average within-range human population density at the year in which the species was discovered
# ActivityFamily: average within-range taxonomic activity per family at the year in which the species was discovered
# ActivityBioregion: average within-range taxonomic activity per bioregion at the year in which the species was discovered

# Filter the species data to comprise those described between 1759 and 2014:
trait_data<-trait_data[trait_data$YearOfDescription>=1759 & trait_data$YearOfDescription<=2014,]

# Get the response variable (time-to-event) standardized between 0 and 1:
trait_data$Time<-(trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
trait_data$Censor<-1 # Censor variable (it informs if the event happened)

# Filter the dataset to represent only one vertebrate group:
trait_data<-trait_data[trait_data$Class=="Aves",] # choose one ("Amphibia", "Reptilia", "Mammalia", "Aves")
trait_data<-droplevels(trait_data) # drop unused levels

# Remove species without data available:
trait_data<-trait_data[complete.cases(trait_data), ]

# Log10 transform and standardize predictors (mean 0 and SD = 1):
trackthis<-ncol(trait_data)
trait_data$LogBodySize<-scale(log10(trait_data$BodySize), center=T, scale=T)
trait_data$LogRangeSize<-scale(log10(trait_data$RangeSize), center=T, scale=T)
trait_data$LogRarity<-scale(log10(trait_data$RangeRarity+1), center=T, scale=T)
trait_data$LogAPP<-scale(log10(trait_data$AnnuPrecip+1), center=T, scale=T)
trait_data$LogAMT<-scale(log10((trait_data$AnnuMeanTemp+2730)/10), center=T, scale=T) # convert from Kelvin to Celsius degree
trait_data$LogTS<-scale(log10(trait_data$TempSeasonality+1), center=T, scale=T)
trait_data$LogPS<-scale(log10(trait_data$PrecipSeasonality+1), center=T, scale=T)
trait_data$LogElevM<-scale(log10(trait_data$Elevation+1), center=T, scale=T)
trait_data$LogPopD<-scale(log10(trait_data$HumanDensity+1), center=T, scale=T)
trait_data$LogTaxPerClade<-scale(log10(trait_data$ActivityFamily+1), center=T, scale=T)
trait_data$LogTaxPerBiome<-scale(log10(trait_data$ActivityBioregion+1), center=T, scale=T)

# Check multicollinearity among the predictor variables:
usdm::vif(trait_data[,trackthis:ncol(trait_data)])

#####

# STEP 2 - CHECK DIFFERENT FAMILY ERROR DISTRIBUTION AND CREATE MODEL FORMULAS WITH ALL POSSIBLE PREDICTOR COMBINATIONS
##########################################################################################################################
# STEP 2 - CHECK DIFFERENT FAMILY ERROR DISTRIBUTION AND CREATE MODEL FORMULAS WITH ALL POSSIBLE PREDICTOR COMBINATIONS
rm(list=setdiff(ls(), c("trait_data")))

# Identify the best error distribution to be used for the Accelerated Failure Time (AFT) Model:
Model_fit<-as.data.frame(matrix(ncol=3, nrow=6))
names(Model_fit)<-(c("Distrib", "npars", "BIC"))
Model_fit[,1]<-c("exp", "weibull", "lnorm", "llogis", "gamma", "gompertz")
Model_fit[1,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="exp", data=trait_data))$npars # get the numbers of parameters
Model_fit[1,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="exp", data=trait_data)) # the the BIC of the model
Model_fit[2,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="weibull", data=trait_data))$npars
Model_fit[2,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="weibull", data=trait_data))
Model_fit[3,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="lnorm", data=trait_data))$npars
Model_fit[3,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="lnorm", data=trait_data))
Model_fit[4,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="llogis", data=trait_data))$npars
Model_fit[4,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="llogis", data=trait_data))
Model_fit[5,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="gamma", data=trait_data))$npars
Model_fit[5,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="gamma", data=trait_data))
Model_fit[6,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="gompertz", data=trait_data))$npars
Model_fit[6,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="gompertz", data=trait_data))

# Verify which error distribution best fits the data:
Model_fit$deltaBIC<-Model_fit[,3]-min(Model_fit$BIC, na.rm=T)
Model_fit$wBIC<-Weights(Model_fit$BIC)
Model_fit<-Model_fit[order(Model_fit$BIC, decreasing=F),]
Model_fit
# Conclusion: when using the full dataset, the Gompertz family was identifed for amphibians, reptiles, and mammals, and the Weibull for birds.
# Sensitivity analysis on the influence of time period of species discovery revealed that the Gompertz family is also predominantly selected
# for birds (if earlier described species are discarded). To allow cross group comparisons, we kept Gompertz family for all groups.
# Results were qualitatively similar when running AFT models for birds with the Weibull family.

# Create model formulas with all possible combinations of predictor variables:
my_predictors<-c("LogBodySize", "LogRangeSize", "LogRarity", "LogAPP", "LogAMT", "LogTS", "LogPS", "LogElevM", "LogPopD", "LogTaxPerClade", "LogTaxPerBiome")
model_combinations<-list()
for(i in 1:length(my_predictors)){ # create all predictor combinations
  model_combinations[[i]]<-combn(my_predictors, i, simplify=FALSE)
} 

# Merge all lists of predictor combinations:
model_combinations<-c(model_combinations[[1]], model_combinations[[2]], model_combinations[[3]], model_combinations[[4]],
                      model_combinations[[5]], model_combinations[[6]], model_combinations[[7]], model_combinations[[8]],
                      model_combinations[[9]], model_combinations[[10]], model_combinations[[11]])

# Create and store model formulas for each predictor combination:
formlist<-list() 
for(i in 1:length(model_combinations)){ # create model formulas with the combinations of explanatory variables
  formlist[[i]]<-as.formula(paste("Surv(Time, Censor)~", paste((model_combinations)[[i]], collapse= "+")))
}
formlist[c(1,100,1000,2047)] # see what it looks like

#####

# STEP 3 - RUN THE MODEL AVERAGING PROCEDURE USING ALL POSSIBLE THE ACCELERATED FAILURE TIME (AFT) MODELS
##########################################################################################################################
# STEP 3 - RUN THE MODEL AVERAGING PROCEDURE USING ALL POSSIBLE THE ACCELERATED FAILURE TIME (AFT) MODELS
rm(list=setdiff(ls(), c("trait_data", "formlist", "Model_fit")))

# Create a dataframe to store the outputs of the survival model averaging:
list_of_families<-list("exponential", "weibull", "lognormal", "llogis", "gamma", "gompertz")

# Create a dataframe to store the coefficients, standard errors, and BIC of each AFT model:
set_of_models<-as.data.frame(matrix(nrow=length(formlist), ncol=28))  
colnames(set_of_models)<-(c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",
                            "LogTS", "LogPS", "LogElevM", "LogPopD", "LogTaxPerBiome", "LogRarity",
                            "LogBodySize_SE", "LogTaxPerClade_SE", "LogRangeSize_SE", "LogAMT_SE", "LogAPP_SE",
                            "LogTS_SE", "LogPS_SE", "LogElevM_SE", "LogPopD_SE", "LogTaxPerBiome_SE", "LogRarity_SE",
                            "Parameter1", "Parameter1_SE", "Parameter2", "Parameter2_SE", "FamilyError", "BIC"))

# Same as above, but to store upper and lower limits of the coefficients:
set_of_models_L95<-as.data.frame(matrix(nrow=length(formlist), ncol=15))  
colnames(set_of_models_L95)<-c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",  "LogTS", "LogPS", "LogElevM",
                               "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2", "FamilyError", "BIC")
set_of_models_U95<-as.data.frame(matrix(nrow=length(formlist), ncol=15))  
colnames(set_of_models_U95)<-c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",  "LogTS", "LogPS", "LogElevM",
                               "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2", "FamilyError", "BIC")

# Remember the best family error distribution or mannualy set one of your choice:
bestfamily<-Model_fit$Distrib[1]
# Or define it mannualy here ("exponential", "weibull", "lognormal", "llogis", "gamma", "gompertz")
# bestfamily<-"gompertz"

# Extract the BIC, and standardized coefficents for each AFT model:
for(i in 1:length(formlist)){
  
  fit<-flexsurvreg(formlist[[i]], dist=bestfamily, data=trait_data) # build the AFT model
  set_of_models[i, (ncol(set_of_models)-1)]<-bestfamily # store the family error distribution (the same for all models)
  set_of_models_L95[i, (ncol(set_of_models_L95)-1)]<-bestfamily # store the family error distribution (the same for all models)
  set_of_models_U95[i, (ncol(set_of_models_U95)-1)]<-bestfamily # store the family error distribution (the same for all models)
  set_of_models[i,ncol(set_of_models)]<-BIC(fit) # store the BIC for the model i
  set_of_models_L95[i,ncol(set_of_models_L95)]<-BIC(fit) # store the BIC for the model i
  set_of_models_U95[i,ncol(set_of_models_U95)]<-BIC(fit) # store the BIC for the model i
  model_coefs<-fit$res # get the standardized coefficients
  
  # Store the coefficient of predictors included in the AFT model i:
  if (length(model_coefs[which(row.names(model_coefs)=="LogBodySize"),1])==1){
    set_of_models[i,1]<-model_coefs[(row.names(model_coefs)=="LogBodySize"),1]
    set_of_models_L95[i,1]<-model_coefs[(row.names(model_coefs)=="LogBodySize"),2]
    set_of_models_U95[i,1]<-model_coefs[(row.names(model_coefs)=="LogBodySize"),3]
    }
  if (length(model_coefs[which(row.names(model_coefs)=="LogTaxPerClade"),1])==1){
    set_of_models[i,2]<-model_coefs[(row.names(model_coefs)=="LogTaxPerClade"),1]
    set_of_models_L95[i,2]<-model_coefs[(row.names(model_coefs)=="LogTaxPerClade"),2]
    set_of_models_U95[i,2]<-model_coefs[(row.names(model_coefs)=="LogTaxPerClade"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogRangeSize"),1])==1){
    set_of_models[i,3]<-model_coefs[(row.names(model_coefs)=="LogRangeSize"),1]
    set_of_models_L95[i,3]<-model_coefs[(row.names(model_coefs)=="LogRangeSize"),2]
    set_of_models_U95[i,3]<-model_coefs[(row.names(model_coefs)=="LogRangeSize"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogAMT"),1])==1){
    set_of_models[i,4]<-model_coefs[(row.names(model_coefs)=="LogAMT"),1]
    set_of_models_L95[i,4]<-model_coefs[(row.names(model_coefs)=="LogAMT"),2]
    set_of_models_U95[i,4]<-model_coefs[(row.names(model_coefs)=="LogAMT"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogAPP"),1])==1){
    set_of_models[i,5]<-model_coefs[(row.names(model_coefs)=="LogAPP"),1]
    set_of_models_L95[i,5]<-model_coefs[(row.names(model_coefs)=="LogAPP"),2]
    set_of_models_U95[i,5]<-model_coefs[(row.names(model_coefs)=="LogAPP"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogTS"),1])==1){
    set_of_models[i,6]<-model_coefs[(row.names(model_coefs)=="LogTS"),1]
    set_of_models_L95[i,6]<-model_coefs[(row.names(model_coefs)=="LogTS"),2]
    set_of_models_U95[i,6]<-model_coefs[(row.names(model_coefs)=="LogTS"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogPS"),1])==1){
    set_of_models[i,7]<-model_coefs[(row.names(model_coefs)=="LogPS"),1]
    set_of_models_L95[i,7]<-model_coefs[(row.names(model_coefs)=="LogPS"),2]
    set_of_models_U95[i,7]<-model_coefs[(row.names(model_coefs)=="LogPS"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogElevM"),1])==1){
    set_of_models[i,8]<-model_coefs[(row.names(model_coefs)=="LogElevM"),1]
    set_of_models_L95[i,8]<-model_coefs[(row.names(model_coefs)=="LogElevM"),2]
    set_of_models_U95[i,8]<-model_coefs[(row.names(model_coefs)=="LogElevM"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogPopD"),1])==1){
    set_of_models[i,9]<-model_coefs[(row.names(model_coefs)=="LogPopD"),1]
    set_of_models_L95[i,9]<-model_coefs[(row.names(model_coefs)=="LogPopD"),2]
    set_of_models_U95[i,9]<-model_coefs[(row.names(model_coefs)=="LogPopD"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogTaxPerBiome"),1])==1){
    set_of_models[i,10]<-model_coefs[(row.names(model_coefs)=="LogTaxPerBiome"),1]
    set_of_models_L95[i,10]<-model_coefs[(row.names(model_coefs)=="LogTaxPerBiome"),2]
    set_of_models_U95[i,10]<-model_coefs[(row.names(model_coefs)=="LogTaxPerBiome"),3]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogRarity"),1])==1){
    set_of_models[i,11]<-model_coefs[(row.names(model_coefs)=="LogRarity"),1]
    set_of_models_L95[i,11]<-model_coefs[(row.names(model_coefs)=="LogRarity"),2]
    set_of_models_U95[i,11]<-model_coefs[(row.names(model_coefs)=="LogRarity"),3]
  }
  
  # Store the std. error of predictors included in the AFT model i:
  if (length(model_coefs[which(row.names(model_coefs)=="LogBodySize"),4])==1){
    set_of_models[i,12]<-model_coefs[(row.names(model_coefs)=="LogBodySize"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogTaxPerClade"),4])==1){
    set_of_models[i,13]<-model_coefs[(row.names(model_coefs)=="LogTaxPerClade"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogRangeSize"),4])==1){
    set_of_models[i,14]<-model_coefs[(row.names(model_coefs)=="LogRangeSize"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogAMT"),4])==1){
    set_of_models[i,15]<-model_coefs[(row.names(model_coefs)=="LogAMT"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogAPP"),4])==1){
    set_of_models[i,16]<-model_coefs[(row.names(model_coefs)=="LogAPP"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogTS"),4])==1){
    set_of_models[i,17]<-model_coefs[(row.names(model_coefs)=="LogTS"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogPS"),4])==1){
    set_of_models[i,18]<-model_coefs[(row.names(model_coefs)=="LogPS"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogElevM"),4])==1){
    set_of_models[i,19]<-model_coefs[(row.names(model_coefs)=="LogElevM"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogPopD"),4])==1){
    set_of_models[i,20]<-model_coefs[(row.names(model_coefs)=="LogPopD"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogTaxPerBiome"),4])==1){
    set_of_models[i,21]<-model_coefs[(row.names(model_coefs)=="LogTaxPerBiome"),4]
  }
  if (length(model_coefs[which(row.names(model_coefs)=="LogRarity"),4])==1){
    set_of_models[i,22]<-model_coefs[(row.names(model_coefs)=="LogRarity"),4]
  }
  
  # Get the coefficient and std. error for the scale, shape, and/or rate parameters, when present:
  if (fit$dlist$name=="exp"){
    set_of_models[i,23]<-model_coefs[1,1] # get the coef of parameter 1
    set_of_models[i,24]<-model_coefs[1,4] # get the std. error of parameter 1
    set_of_models_L95[i,12]<-model_coefs[1,2] # get the lower coef of parameter 1
    set_of_models_U95[i,12]<-model_coefs[1,3] # get the upper coef of parameter 1
    }
  if (fit$dlist$name!="exp"){
    set_of_models[i,23]<-model_coefs[1,1] # get the coef of parameter 1
    set_of_models[i,24]<-model_coefs[1,4] # get the std. error of parameter 1
    set_of_models[i,25]<-model_coefs[2,1] # get the coef of parameter 2
    set_of_models[i,26]<-model_coefs[2,1] # get the std. error of parameter 2
    
    set_of_models_L95[i,12]<-model_coefs[1,2] # get the lower coef of parameter 1
    set_of_models_L95[i,13]<-model_coefs[2,2] # get the lower coef of parameter 2
    set_of_models_U95[i,12]<-model_coefs[1,3] # get the upper coef of parameter 1
    set_of_models_U95[i,13]<-model_coefs[2,3] # get the upper coef of parameter 2
    
  }
  
  print(i) # to view the progress of the analysis
  
} # end of i for loop across formlist elements

# Compute the wBIC:
set_of_models$wBIC<-Weights(set_of_models$BIC) # compute the wBIC
set_of_models_L95$wBIC<-Weights(set_of_models_L95$BIC) # compute the wBIC
set_of_models_U95$wBIC<-Weights(set_of_models_U95$BIC) # compute the wBIC

# Register the taxonomic class.
set_of_models$Class<-unique(trait_data$Class)
set_of_models_L95$Class<-unique(trait_data$Class)
set_of_models_U95$Class<-unique(trait_data$Class)

# Export the results:
dir.create("ModelAveraging", showWarnings=F)

if(unique(set_of_models$Class)=="Amphibia"){
    fwrite(set_of_models, "ModelAveraging/set_of_models_Amphibia.csv")
    fwrite(set_of_models_L95, "ModelAveraging/set_of_models_L95_Amphibia.csv")
    fwrite(set_of_models_U95, "ModelAveraging/set_of_models_U95_Amphibia.csv")
    }

if(unique(set_of_models$Class)=="Reptilia"){
    fwrite(set_of_models, "ModelAveraging/set_of_models_Reptilia.csv")
    fwrite(set_of_models_L95, "ModelAveraging/set_of_models_L95_Reptilia.csv")
    fwrite(set_of_models_U95, "ModelAveraging/set_of_models_U95_Reptilia.csv")
  }
  
if(unique(set_of_models$Class)=="Mammalia"){
    fwrite(set_of_models, "ModelAveraging/set_of_models_Mammalia.csv")
    fwrite(set_of_models_L95, "ModelAveraging/set_of_models_L95_Mammalia.csv")
    fwrite(set_of_models_U95, "ModelAveraging/set_of_models_U95_Mammalia.csv")
  }
  
if(unique(set_of_models$Class)=="Aves"){
    fwrite(set_of_models, "ModelAveraging/set_of_models_Aves.csv")
    fwrite(set_of_models_L95, "ModelAveraging/set_of_models_L95_Aves.csv")
    fwrite(set_of_models_U95, "ModelAveraging/set_of_models_U95_Aves.csv")
  }

#####

# STEP 4. GET THE AVERAGE WEIGHTED COEFFICIENTS USING THE OUTPUTS FROM THE MODEL AVERAGING PROCEDURE
##########################################################################################################################
# STEP 4. GET THE AVERAGE WEIGHTED COEFFICIENTS USING THE OUTPUTS FROM THE MODEL AVERAGING PROCEDURE
rm(list=setdiff(ls(), c("trait_data")))

# Re-read files (if re-running particular code sections):
my_csv_file<-paste0("ModelAveraging/set_of_models_", unique(trait_data$Class), ".csv")
set_of_models<-read.csv(my_csv_file)

# Get the difference in BIC between and best model and all the others:
#set_of_models<-na.omit(set_of_models)
set_of_models<-set_of_models[order(set_of_models$wBIC, decreasing=T),] # reorder according to wBIC
set_of_models$delta_BIC<-set_of_models$BIC - min(set_of_models$BIC, na.rm=T) # compute the delta BIC
set_of_models<-as.data.frame(set_of_models) # convert from data.table to data.frame

# Get the averaged model coefficients based on wBIC:
avg_model<-as.data.frame(matrix(nrow=5, ncol=13))
colnames(avg_model)<-(c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP", "LogTS", "LogPS",
                        "LogElevM", "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2"))
row.names(avg_model)<-c("Avg.coef","Std.error", "Lower_IC", "Upper_IC", "Label")

for(i in 1:11){ # averaged weighted std. coefficients for predictors
  candidate_models<-set_of_models[which(set_of_models[,i]!=0),]  # select all models for a given predictor
  candidate_models$wBIC<-Weights(candidate_models$BIC) # rescale the wBIC
  avg_model[1,i]<-sum(candidate_models[,i]*candidate_models$wBIC) # standardized coefficients of the weighted model
} # averaged coefficients

# Compute the unconditional variance estimator for the averaged model coefficients:
for(i in 12:22){
  candidate_models<-set_of_models[which(set_of_models[,i]!=0),] # select all models for a given predictor
  weighted_var<-as.data.frame(matrix(nrow=nrow(candidate_models), ncol=1)) 
  
  for(j in 1:nrow(candidate_models)){ # j models, varying from 1 to 2^(n_covariates - 1) 
    delta_coef<-(candidate_models[j,(i-11)] - avg_model[1,(i-11)])^2 # error due to model selection uncertainty
    var_error<-(candidate_models[j,i])^2 # error in parameter estimation of the model 'j'
    weighted_var[j,1]<-(sqrt(var_error + delta_coef))*candidate_models$wBIC[j]
  }
  
  avg_model[2,(i-11)]<-sum(weighted_var[,1]) # unconditional standard error
} 

# Compute the averaged weighted coefficients for the shape/rate/scale parameters (when present):
avg_model[1,12]<-sum(set_of_models[,23]*set_of_models$wBIC)
if(unique(set_of_models$FamilyError)!="exp"){avg_model[1,13]<-sum(set_of_models[,25]*set_of_models$wBIC)}

# Compute the unconditional variance estimator for the rate/shape/scale parameters (when present):
weighted_var<-as.data.frame(matrix(nrow=nrow(set_of_models), ncol=1)) 
for(j in 1:nrow(set_of_models)){ # j models, varying from 1 to 2^(n_covariates - 1) 
  delta_coef<-(set_of_models[j,23] - avg_model[1,12])^2 # error due to model selection uncertainty
  var_error<-(set_of_models[j,24])^2 # error in parameter estimation of the model 'j'
  weighted_var[j,1]<-(sqrt(var_error + delta_coef))*set_of_models$wBIC[j]
}
avg_model[2,12]<-sum(weighted_var[,1]) # unconditional standard error

# Same as above, but for the second parameter (not present for exponential family error distribution):
if(unique(set_of_models$FamilyError)!="exp"){
  weighted_var<-as.data.frame(matrix(nrow=nrow(set_of_models), ncol=1)) 
  for(j in 1:nrow(set_of_models)){ # j models, varying from 1 to 2^(n_covariates - 1) 
    delta_coef<-(set_of_models[j,25] - avg_model[1,13])^2 # error due to model selection uncertainty
    var_error<-(set_of_models[j,26])^2 # error in parameter estimation of the model 'j'
    weighted_var[j,1]<-(sqrt(var_error + delta_coef))*set_of_models$wBIC[j]
    }
  avg_model[2,13]<-sum(weighted_var[,1]) # unconditional standard error
  }

for(i in 1:13){ # unconditional confidence intervals
  avg_model[3,i]<-avg_model[1,i]-(1.96*avg_model[2,i]) # lower confidence interval
  avg_model[4,i]<-avg_model[1,i]+(1.96*avg_model[2,i]) # upper confidence interval
} # unconditional confidence intervals

# Prepare the data.frame for posterior visualization and saved it:
avg_model<-as.data.frame(t(avg_model))
avg_model$Label<-as.factor(row.names(avg_model))

# Register the taxonomic class:
avg_model$Class<-unique(set_of_models$Class)

# Export the average weighted coefficients:
if(unique(avg_model$Class)=="Amphibia"){write.csv(avg_model,"ModelAveraging/AvgModel_Coefs_Amphibia.csv", row.names=F, fileEncoding="UTF-8")}
if(unique(avg_model$Class)=="Reptilia"){write.csv(avg_model,"ModelAveraging/AvgModel_Coefs_Reptilia.csv", row.names=F, fileEncoding="UTF-8")}
if(unique(avg_model$Class)=="Mammalia"){write.csv(avg_model,"ModelAveraging/AvgModel_Coefs_Mammalia.csv", row.names=F, fileEncoding="UTF-8")}
if(unique(avg_model$Class)=="Aves"){write.csv(avg_model,"ModelAveraging/AvgModel_Coefs_Aves.csv", row.names=F, fileEncoding="UTF-8")}

#####

# STEP 5 - GET SPECIES-LEVEL PREDICTIONS OF DISCOVERY PROBABILITY AND ESTIMATED YEAR OF DISCOVERY FOR ALL AFT MODELS
##########################################################################################################################
# STEP 5 - GET SPECIES-LEVEL PREDICTIONS OF DISCOVERY PROBABILITY AND ESTIMATED YEAR OF DISCOVERY FOR ALL AFT MODELS
rm(list=setdiff(ls(), c("trait_data")))
# Re-run STEP1 if the obect 'trait_data' is not in your workspace 

# Load the table containg all model coefs:
my_csv_file<-paste0("ModelAveraging/set_of_models_", unique(trait_data$Class), ".csv")
my_csv_file_L95<-paste0("ModelAveraging/set_of_models_L95_", unique(trait_data$Class), ".csv")
my_csv_file_U95<-paste0("ModelAveraging/set_of_models_U95_", unique(trait_data$Class), ".csv")
set_of_models<-read.csv(my_csv_file)
set_of_models_L95<-read.csv(my_csv_file_L95)
set_of_models_U95<-read.csv(my_csv_file_U95)

# Replace NA values by 0 in the 'set_of_models' objects:
set_of_models[is.na(set_of_models)] <- 0
set_of_models_L95[is.na(set_of_models_L95)] <- 0
set_of_models_U95[is.na(set_of_models_U95)] <- 0

# Define an object with percentile subdivisions (to be used to latter):
pct<-seq(0, 1, by=0.004) # set the by argument to 0.01 to reduce computational cost
# The 260 years of discovery exploration will not be accurately reflected if pct has much less than 260 units

# Create an empty dataframe to store the estimated discovery probability (95% CI) in the year 2015 (all species, all AFT models):
DiscProb_dt<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data)))
DiscProb_dt_L95<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data)))
DiscProb_dt_U95<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data)))
DiscProb_dt[is.na(DiscProb_dt)] = 0
DiscProb_dt_L95[is.na(DiscProb_dt_L95)] = 0
DiscProb_dt_U95[is.na(DiscProb_dt_U95)] = 0

# Same as above, but to store the predicted year of discovery (= the year in which discovery probability is 0.50):
YearPred_dt<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data))) 
YearPred_dt[is.na(YearPred_dt)] = 0

# Get the discovery probability (95% CI) in 2015 and estimated year of discovery for all species and AFT models:
cl<-makePSOCKcluster(detectCores()-1)
registerDoParallel(cl)
getDoParWorkers()
a<-Sys.time()
for(i in 1:nrow(set_of_models)){
  
  if(unique(set_of_models$FamilyError)=="exp"){
    TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  qexp(pct, rate=exp(as.numeric(set_of_models[i,23] +
                                                                  trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                  trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                  trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                  trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                  trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                  trait_data$LogTS[k]*set_of_models[i,6] +
                                                                  trait_data$LogPS[k]*set_of_models[i,7] +
                                                                  trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                  trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                  trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                  trait_data$LogRarity[k]*set_of_models[i,11])))
                                }
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pexp(pct, rate=exp(as.numeric(set_of_models[i,23] +
                                                                trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                trait_data$LogTS[k]*set_of_models[i,6] +
                                                                trait_data$LogPS[k]*set_of_models[i,7] +
                                                                trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                trait_data$LogRarity[k]*set_of_models[i,11])))
                              }
    ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pexp(pct, rate=exp(as.numeric(set_of_models_L95[i,12] +
                                                                    trait_data$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                    trait_data$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                    trait_data$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                    trait_data$LogAMT[k]*set_of_models_L95[i,4] +
                                                                    trait_data$LogAPP[k]*set_of_models_L95[i,5] +
                                                                    trait_data$LogTS[k]*set_of_models_L95[i,6] +
                                                                    trait_data$LogPS[k]*set_of_models_L95[i,7] +
                                                                    trait_data$LogElevM[k]*set_of_models_L95[i,8] +
                                                                    trait_data$LogPopD[k]*set_of_models_L95[i,9] +
                                                                    trait_data$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                    trait_data$LogRarity[k]*set_of_models_L95[i,11])))
                                  }
    ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pexp(pct, rate=exp(as.numeric(set_of_models_U95[i,12] +
                                                                    trait_data$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                    trait_data$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                    trait_data$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                    trait_data$LogAMT[k]*set_of_models_U95[i,4] +
                                                                    trait_data$LogAPP[k]*set_of_models_U95[i,5] +
                                                                    trait_data$LogTS[k]*set_of_models_U95[i,6] +
                                                                    trait_data$LogPS[k]*set_of_models_U95[i,7] +
                                                                    trait_data$LogElevM[k]*set_of_models_U95[i,8] +
                                                                    trait_data$LogPopD[k]*set_of_models_U95[i,9] +
                                                                    trait_data$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                    trait_data$LogRarity[k]*set_of_models_U95[i,11])))
                                  }
      }
  
  if(unique(set_of_models$FamilyError)=="weibull"){
    TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  qweibull(pct, shape=set_of_models[i,23], 
                                           scale=as.numeric(set_of_models[i,25]*
                                                              exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                    trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                    trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                    trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                    trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                    trait_data$LogTS[k]*set_of_models[i,6] +
                                                                    trait_data$LogPS[k]*set_of_models[i,7] +
                                                                    trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                    trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                    trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                    trait_data$LogRarity[k]*set_of_models[i,11])))
                                }
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pweibull(pct, shape=set_of_models[i,23], 
                                         scale=as.numeric(set_of_models[i,25]*
                                                            exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                  trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                  trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                  trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                  trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                  trait_data$LogTS[k]*set_of_models[i,6] +
                                                                  trait_data$LogPS[k]*set_of_models[i,7] +
                                                                  trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                  trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                  trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                  trait_data$LogRarity[k]*set_of_models[i,11])))
                              }
    ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pweibull(pct, shape=set_of_models_L95[i,12], 
                                             scale=as.numeric(set_of_models_L95[i,13]*
                                                                exp(trait_data$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                      trait_data$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                      trait_data$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                      trait_data$LogAMT[k]*set_of_models_L95[i,4] +
                                                                      trait_data$LogAPP[k]*set_of_models_L95[i,5] +
                                                                      trait_data$LogTS[k]*set_of_models_L95[i,6] +
                                                                      trait_data$LogPS[k]*set_of_models_L95[i,7] +
                                                                      trait_data$LogElevM[k]*set_of_models_L95[i,8] +
                                                                      trait_data$LogPopD[k]*set_of_models_L95[i,9] +
                                                                      trait_data$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                      trait_data$LogRarity[k]*set_of_models_L95[i,11])))
                                  }
    ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pweibull(pct, shape=set_of_models_U95[i,12], 
                                             scale=as.numeric(set_of_models_U95[i,13]*
                                                                exp(trait_data$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                      trait_data$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                      trait_data$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                      trait_data$LogAMT[k]*set_of_models_U95[i,4] +
                                                                      trait_data$LogAPP[k]*set_of_models_U95[i,5] +
                                                                      trait_data$LogTS[k]*set_of_models_U95[i,6] +
                                                                      trait_data$LogPS[k]*set_of_models_U95[i,7] +
                                                                      trait_data$LogElevM[k]*set_of_models_U95[i,8] +
                                                                      trait_data$LogPopD[k]*set_of_models_U95[i,9] +
                                                                      trait_data$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                      trait_data$LogRarity[k]*set_of_models_U95[i,11])))
    
                                  }
  }
  
  if(unique(set_of_models$FamilyError)=="lnorm"){
    TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  qlnorm(pct, meanlog=as.numeric(set_of_models[i,23]*
                                                                   exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                         trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                         trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                         trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                         trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                         trait_data$LogTS[k]*set_of_models[i,6] +
                                                                         trait_data$LogPS[k]*set_of_models[i,7] +
                                                                         trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                         trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                         trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                         trait_data$LogRarity[k]*set_of_models[i,11],
                                                                       sdlog=set_of_models[i,25])))
                                }
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                plnorm(pct, meanlog=as.numeric(set_of_models[i,23]*
                                                                 exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                       trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                       trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                       trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                       trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                       trait_data$LogTS[k]*set_of_models[i,6] +
                                                                       trait_data$LogPS[k]*set_of_models[i,7] +
                                                                       trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                       trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                       trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                       trait_data$LogRarity[k]*set_of_models[i,11],
                                                                     sdlog=set_of_models[i,25])))
                              }
    ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    plnorm(pct, meanlog=as.numeric(set_of_models_L95[i,12]*
                                                                     exp(trait_data$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                           trait_data$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                           trait_data$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                           trait_data$LogAMT[k]*set_of_models_L95[i,4] +
                                                                           trait_data$LogAPP[k]*set_of_models_L95[i,5] +
                                                                           trait_data$LogTS[k]*set_of_models_L95[i,6] +
                                                                           trait_data$LogPS[k]*set_of_models_L95[i,7] +
                                                                           trait_data$LogElevM[k]*set_of_models_L95[i,8] +
                                                                           trait_data$LogPopD[k]*set_of_models_L95[i,9] +
                                                                           trait_data$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                           trait_data$LogRarity[k]*set_of_models_L95[i,11],
                                                                         sdlog=set_of_models_L95[i,13])))
                                  }
    ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    plnorm(pct, meanlog=as.numeric(set_of_models_U95[i,12]*
                                                                     exp(trait_data$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                           trait_data$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                           trait_data$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                           trait_data$LogAMT[k]*set_of_models_U95[i,4] +
                                                                           trait_data$LogAPP[k]*set_of_models_U95[i,5] +
                                                                           trait_data$LogTS[k]*set_of_models_U95[i,6] +
                                                                           trait_data$LogPS[k]*set_of_models_U95[i,7] +
                                                                           trait_data$LogElevM[k]*set_of_models_U95[i,8] +
                                                                           trait_data$LogPopD[k]*set_of_models_U95[i,9] +
                                                                           trait_data$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                           trait_data$LogRarity[k]*set_of_models_U95[i,11],
                                                                         sdlog=set_of_models_U95[i,13])))
                                  }
   
  }
  
  if(unique(set_of_models$FamilyError)=="llogis"){
    TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  qllogis(pct, shape=set_of_models[i,23], 
                                          scale=as.numeric(set_of_models[i,25]*
                                                             exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                   trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                   trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                   trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                   trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                   trait_data$LogTS[k]*set_of_models[i,6] +
                                                                   trait_data$LogPS[k]*set_of_models[i,7] +
                                                                   trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                   trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                   trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                   trait_data$LogRarity[k]*set_of_models[i,11])))
                                }
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pllogis(pct, shape=set_of_models[i,23], 
                                        scale=as.numeric(set_of_models[i,25]*
                                                           exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                 trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                 trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                 trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                 trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                 trait_data$LogTS[k]*set_of_models[i,6] +
                                                                 trait_data$LogPS[k]*set_of_models[i,7] +
                                                                 trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                 trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                 trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                 trait_data$LogRarity[k]*set_of_models[i,11])))
                              }
    ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pllogis(pct, shape=set_of_models_L95[i,12], 
                                            scale=as.numeric(set_of_models_L95[i,13]*
                                                               exp(trait_data$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                     trait_data$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                     trait_data$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                     trait_data$LogAMT[k]*set_of_models_L95[i,4] +
                                                                     trait_data$LogAPP[k]*set_of_models_L95[i,5] +
                                                                     trait_data$LogTS[k]*set_of_models_L95[i,6] +
                                                                     trait_data$LogPS[k]*set_of_models_L95[i,7] +
                                                                     trait_data$LogElevM[k]*set_of_models_L95[i,8] +
                                                                     trait_data$LogPopD[k]*set_of_models_L95[i,9] +
                                                                     trait_data$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                     trait_data$LogRarity[k]*set_of_models_L95[i,11])))
                                  }
    ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pllogis(pct, shape=set_of_models_U95[i,12], 
                                            scale=as.numeric(set_of_models_U95[i,13]*
                                                               exp(trait_data$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                     trait_data$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                     trait_data$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                     trait_data$LogAMT[k]*set_of_models_U95[i,4] +
                                                                     trait_data$LogAPP[k]*set_of_models_U95[i,5] +
                                                                     trait_data$LogTS[k]*set_of_models_U95[i,6] +
                                                                     trait_data$LogPS[k]*set_of_models_U95[i,7] +
                                                                     trait_data$LogElevM[k]*set_of_models_U95[i,8] +
                                                                     trait_data$LogPopD[k]*set_of_models_U95[i,9] +
                                                                     trait_data$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                     trait_data$LogRarity[k]*set_of_models_U95[i,11])))
                                  }
  }
  
  if(unique(set_of_models$FamilyError)=="gamma"){
    TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  qgamma(pct, shape=set_of_models[i,23], 
                                         rate=as.numeric(set_of_models[i,25]*
                                                           exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                 trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                 trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                 trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                 trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                 trait_data$LogTS[k]*set_of_models[i,6] +
                                                                 trait_data$LogPS[k]*set_of_models[i,7] +
                                                                 trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                 trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                 trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                 trait_data$LogRarity[k]*set_of_models[i,11])))
                                }
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pgamma(pct, shape=set_of_models[i,23], 
                                       rate=as.numeric(set_of_models[i,25]*
                                                         exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                               trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                               trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                               trait_data$LogAMT[k]*set_of_models[i,4] +
                                                               trait_data$LogAPP[k]*set_of_models[i,5] +
                                                               trait_data$LogTS[k]*set_of_models[i,6] +
                                                               trait_data$LogPS[k]*set_of_models[i,7] +
                                                               trait_data$LogElevM[k]*set_of_models[i,8] +
                                                               trait_data$LogPopD[k]*set_of_models[i,9] +
                                                               trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                               trait_data$LogRarity[k]*set_of_models[i,11])))
                              }
    ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pgamma(pct, shape=set_of_models_L95[i,12], 
                                           rate=as.numeric(set_of_models_L95[i,13]*
                                                             exp(trait_data$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                   trait_data$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                   trait_data$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                   trait_data$LogAMT[k]*set_of_models_L95[i,4] +
                                                                   trait_data$LogAPP[k]*set_of_models_L95[i,5] +
                                                                   trait_data$LogTS[k]*set_of_models_L95[i,6] +
                                                                   trait_data$LogPS[k]*set_of_models_L95[i,7] +
                                                                   trait_data$LogElevM[k]*set_of_models_L95[i,8] +
                                                                   trait_data$LogPopD[k]*set_of_models_L95[i,9] +
                                                                   trait_data$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                   trait_data$LogRarity[k]*set_of_models_L95[i,11])))
                                  }
    ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data),
                                  .combine = 'cbind',
                                  .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                    pgamma(pct, shape=set_of_models_U95[i,12], 
                                           rate=as.numeric(set_of_models_U95[i,13]*
                                                             exp(trait_data$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                   trait_data$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                   trait_data$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                   trait_data$LogAMT[k]*set_of_models_U95[i,4] +
                                                                   trait_data$LogAPP[k]*set_of_models_U95[i,5] +
                                                                   trait_data$LogTS[k]*set_of_models_U95[i,6] +
                                                                   trait_data$LogPS[k]*set_of_models_U95[i,7] +
                                                                   trait_data$LogElevM[k]*set_of_models_U95[i,8] +
                                                                   trait_data$LogPopD[k]*set_of_models_U95[i,9] +
                                                                   trait_data$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                   trait_data$LogRarity[k]*set_of_models_U95[i,11])))
                                  }
    }
  
  if(unique(set_of_models$FamilyError)=="gompertz"){
    TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data),
            .combine = 'cbind',
            .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
              qgompertz(pct, shape=set_of_models[i,23],
                        rate=as.numeric(set_of_models[i,25]*
                                        exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                        trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                        trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                        trait_data$LogAMT[k]*set_of_models[i,4] +
                                        trait_data$LogAPP[k]*set_of_models[i,5] +
                                        trait_data$LogTS[k]*set_of_models[i,6] +
                                        trait_data$LogPS[k]*set_of_models[i,7] +
                                        trait_data$LogElevM[k]*set_of_models[i,8] +
                                        trait_data$LogPopD[k]*set_of_models[i,9] +
                                        trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                        trait_data$LogRarity[k]*set_of_models[i,11])))
            }
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  pgompertz(pct, shape=set_of_models[i,23],
                                            rate=as.numeric(set_of_models[i,25]*
                                                              exp(trait_data$LogBodySize[k]*set_of_models[i,1] +
                                                                    trait_data$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                    trait_data$LogRangeSize[k]*set_of_models[i,3] +
                                                                    trait_data$LogAMT[k]*set_of_models[i,4] +
                                                                    trait_data$LogAPP[k]*set_of_models[i,5] +
                                                                    trait_data$LogTS[k]*set_of_models[i,6] +
                                                                    trait_data$LogPS[k]*set_of_models[i,7] +
                                                                    trait_data$LogElevM[k]*set_of_models[i,8] +
                                                                    trait_data$LogPopD[k]*set_of_models[i,9] +
                                                                    trait_data$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                    trait_data$LogRarity[k]*set_of_models[i,11])))
                                }
    ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pgompertz(pct, shape=set_of_models_L95[i,12],
                                          rate=as.numeric(set_of_models_L95[i,13]*
                                                            exp(trait_data$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                  trait_data$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                  trait_data$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                  trait_data$LogAMT[k]*set_of_models_L95[i,4] +
                                                                  trait_data$LogAPP[k]*set_of_models_L95[i,5] +
                                                                  trait_data$LogTS[k]*set_of_models_L95[i,6] +
                                                                  trait_data$LogPS[k]*set_of_models_L95[i,7] +
                                                                  trait_data$LogElevM[k]*set_of_models_L95[i,8] +
                                                                  trait_data$LogPopD[k]*set_of_models_L95[i,9] +
                                                                  trait_data$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                  trait_data$LogRarity[k]*set_of_models_L95[i,11])))
                              }
    ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pgompertz(pct, shape=set_of_models_U95[i,12],
                                          rate=as.numeric(set_of_models_U95[i,13]*
                                                            exp(trait_data$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                  trait_data$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                  trait_data$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                  trait_data$LogAMT[k]*set_of_models_U95[i,4] +
                                                                  trait_data$LogAPP[k]*set_of_models_U95[i,5] +
                                                                  trait_data$LogTS[k]*set_of_models_U95[i,6] +
                                                                  trait_data$LogPS[k]*set_of_models_U95[i,7] +
                                                                  trait_data$LogElevM[k]*set_of_models_U95[i,8] +
                                                                  trait_data$LogPopD[k]*set_of_models_U95[i,9] +
                                                                  trait_data$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                  trait_data$LogRarity[k]*set_of_models_U95[i,11])))
                              }
    }
  
  # colnames = species & rownames = discovery probability values
  TimeToDisc_per_spp<-as.data.frame(TimeToDisc_per_spp)
  names(TimeToDisc_per_spp)<-trait_data$Binomial; rownames(TimeToDisc_per_spp)<-pct
  
  # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
  ProbDisc_per_spp<-as.data.frame(ProbDisc_per_spp)
  ProbDisc_per_spp_L95<-as.data.frame(ProbDisc_per_spp_L95)
  ProbDisc_per_spp_U95<-as.data.frame(ProbDisc_per_spp_U95)
  names(ProbDisc_per_spp)<-trait_data$Binomial; rownames(ProbDisc_per_spp)<-pct 
  names(ProbDisc_per_spp_L95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
  names(ProbDisc_per_spp_U95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
  
  # Extract the discovery probability in 2015 for all species: 
  present_year<-(2015-1758)/(260) # get the year 2015 in a scale of 0 to 1 (0 = 1758, 1 = 2018)
  row_selected<-which(row.names(ProbDisc_per_spp)>=present_year)[1] # select the row of the ProbDisc_per_spp object the represents the year 2015
  ProbDisc_2015<-t(ProbDisc_per_spp[row_selected,]) # extract the discovery probability in 2015 according to the model i
  ProbDisc_2015_L95<-t(ProbDisc_per_spp_L95[row_selected,]) # extract the lower bound of the discovery probability in 2015 according to the model i
  ProbDisc_2015_U95<-t(ProbDisc_per_spp_U95[row_selected,]) # extract the upper bound of the discovery probability in 2015 according to the model i
  
  #and the predicted year for the discovery probability of 0.50:
  TimeToDisc_per_spp<-(TimeToDisc_per_spp*(2018-1758))+1758 # convert row values to year
  threshold_selected<-which(row.names(TimeToDisc_per_spp)>=0.50) # select the rows with discovery probability >=0.50
  thershold_selected<-threshold_selected[1] # select the first row informing the discovery probability >=0.50
  YearPred<-t(TimeToDisc_per_spp[threshold_selected, ]) # year predicted for a given species description (99.5% of probability)
  
  # Update the 'DiscProb_dt' and 'YearPred_dt' objects with the discovery probability and estimated year of discovery for each species (k) according to the model (i):
  for (k in 1:nrow(trait_data)){set(DiscProb_dt, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015[k])}
  for (k in 1:nrow(trait_data)){set(DiscProb_dt_L95, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015_L95[k])}
  for (k in 1:nrow(trait_data)){set(DiscProb_dt_U95, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015_U95[k])}
  for (k in 1:nrow(trait_data)){set(YearPred_dt, i=as.integer(i), j=as.integer(k), value=YearPred[k])}
  
} # end of i for loop (across all AFT models)
stopCluster(cl)

# How long it took:
b<-Sys.time()
b-a # about 10 hours using in a notebook with 2.11 Ghz, 24GB RAM, 8 cores.

# Register the taxonomic class and export the average weighted coefficients:
DiscProb_dt$Class<-unique(set_of_models$Class)
DiscProb_dt_L95$Class<-unique(set_of_models$Class)
DiscProb_dt_U95$Class<-unique(set_of_models$Class)
YearPred_dt$Class<-unique(set_of_models$Class)

# Export the outputs
if(unique(DiscProb_dt$Class)=="Amphibia"){
  fwrite(DiscProb_dt, "ModelAveraging/DiscProb_dt_Amphibia.csv")
  fwrite(DiscProb_dt_L95, "ModelAveraging/DiscProb_dt_L95_Amphibia.csv")
  fwrite(DiscProb_dt_U95, "ModelAveraging/DiscProb_dt_U95_Amphibia.csv")
  fwrite(YearPred_dt, "ModelAveraging/YearPred_dt_Amphibia.csv")
  }

if(unique(DiscProb_dt$Class)=="Reptilia"){
  fwrite(DiscProb_dt, "ModelAveraging/DiscProb_dt_Reptilia.csv")
  fwrite(DiscProb_dt_L95, "ModelAveraging/DiscProb_dt_L95_Reptilia.csv")
  fwrite(DiscProb_dt_U95, "ModelAveraging/DiscProb_dt_U95_Reptilia.csv")
  fwrite(YearPred_dt, "ModelAveraging/YearPred_dt_Reptilia.csv")
  }

if(unique(DiscProb_dt$Class)=="Mammalia"){
  fwrite(DiscProb_dt, "ModelAveraging/DiscProb_dt_Mammalia.csv")
  fwrite(DiscProb_dt_L95, "ModelAveraging/DiscProb_dt_L95_Mammalia.csv")
  fwrite(DiscProb_dt_U95, "ModelAveraging/DiscProb_dt_U95_Mammalia.csv")
  fwrite(YearPred_dt, "ModelAveraging/YearPred_dt_Mammalia.csv")
  }

if(unique(DiscProb_dt$Class)=="Aves"){
  fwrite(DiscProb_dt, "ModelAveraging/DiscProb_dt_Aves.csv")
  fwrite(DiscProb_dt_L95, "ModelAveraging/DiscProb_dt_L95_Aves.csv")
  fwrite(DiscProb_dt_U95, "ModelAveraging/DiscProb_dt_U95_Aves.csv")
  fwrite(YearPred_dt, "ModelAveraging/YearPred_dt_Aves.csv")
  }

#####

# STEP 6 - GET THE AVERAGE WEIGHTED DISCOVERY PROBABILITY AND ESTIMATED YEAR OF DISCOVERY
##########################################################################################################################
# STEP 6 - GET THE AVERAGE WEIGHTED DISCOVERY PROBABILITY AND ESTIMATED YEAR OF DISCOVERY
rm(list=setdiff(ls(), c("trait_data")))
# Re-run STEP1 if the obect 'trait_data' is not in your workspace 

# Read the data.table containg the BIC of each AFT model:
my_csv_file<-paste0("ModelAveraging/set_of_models_", unique(trait_data$Class), ".csv")
set_of_models<-read.csv(my_csv_file)

# Read the files containing the discovery probability, and estimated year of description for each AFT model:
DiscProb_file<-paste0("ModelAveraging/DiscProb_dt_", unique(trait_data$Class), ".csv")
DiscProb_L95_file<-paste0("ModelAveraging/DiscProb_dt_L95_", unique(trait_data$Class), ".csv")
DiscProb_U95_file<-paste0("ModelAveraging/DiscProb_dt_U95_", unique(trait_data$Class), ".csv")
YearPred_file<-paste0("ModelAveraging/YearPred_dt_", unique(trait_data$Class), ".csv")
DiscProb_dt<-fread(DiscProb_file, stringsAsFactors=TRUE, encoding="UTF-8")
DiscProb_dt_L95<-fread(DiscProb_L95_file, stringsAsFactors=TRUE, encoding="UTF-8")
DiscProb_dt_U95<-fread(DiscProb_U95_file, stringsAsFactors=TRUE, encoding="UTF-8")
YearPred_dt<-fread(YearPred_file, stringsAsFactors=TRUE, encoding="UTF-8")

# Remove 'Class' column before subsequent computations:
DiscProb_dt<-DiscProb_dt[ , Class:=NULL]
DiscProb_dt_L95<-DiscProb_dt_L95[ , Class:=NULL]
DiscProb_dt_U95<-DiscProb_dt_U95[ , Class:=NULL]
YearPred_dt<-YearPred_dt[ , Class:=NULL]

# Get the weigthed values of estimated year of discovery and discovery probability in 2015:
spp_columns<-names(DiscProb_dt)
DiscProb_dt[, (spp_columns) := lapply(.SD, function(x) 
  x * set_of_models$wBIC*100), .SDcols = spp_columns]
AvgW_DiscProb_dt<-colSums(DiscProb_dt)/100

spp_columns<-names(DiscProb_dt_L95)
DiscProb_dt_L95[, (spp_columns) := lapply(.SD, function(x) 
  x * set_of_models$wBIC*100), .SDcols = spp_columns]
AvgW_DiscProb_dt_L95<-colSums(DiscProb_dt_L95)/100

spp_columns<-names(DiscProb_dt_U95)
DiscProb_dt_U95[, (spp_columns) := lapply(.SD, function(x) 
  x * set_of_models$wBIC*100), .SDcols = spp_columns]
AvgW_DiscProb_dt_U95<-colSums(DiscProb_dt_U95)/100

spp_columns<-names(YearPred_dt)
YearPred_dt[, (spp_columns) := lapply(.SD, function(x) 
  x * set_of_models$wBIC*100), .SDcols = spp_columns]
AvgW_YearPred_dt<-colSums(YearPred_dt)/100

# Match the discovery traits with the observed trait data (species are in the same order as in trait_data):
trait_data$DiscProb<-AvgW_DiscProb_dt
trait_data$DiscProb_L95<-AvgW_DiscProb_dt_L95
trait_data$DiscProb_U95<-AvgW_DiscProb_dt_U95
trait_data$YearPred<-AvgW_YearPred_dt

# Filter the columns of interest:
trait_data<-trait_data[, .(Binomial, YearOfDescription, Time, Class, Order, Family, Genus, DiscProb, DiscProb_L95, DiscProb_U95)] # Authority column is missing

# Export the average weighted discovery metrics for the vertebrate group under consideration:
if(unique(trait_data$Class)=="Amphibia"){write.csv(trait_data,"ModelAveraging/SpeciesLevelPredictions_Amphibia.csv", row.names=F, fileEncoding="UTF-8")}
if(unique(trait_data$Class)=="Reptilia"){write.csv(trait_data,"ModelAveraging/SpeciesLevelPredictions_Reptilia.csv", row.names=F, fileEncoding="UTF-8")}
if(unique(trait_data$Class)=="Mammalia"){write.csv(trait_data,"ModelAveraging/SpeciesLevelPredictions_Mammalia.csv", row.names=F, fileEncoding="UTF-8")}
if(unique(trait_data$Class)=="Aves"){write.csv(trait_data,"ModelAveraging/SpeciesLevelPredictions_Aves.csv", row.names=F, fileEncoding="UTF-8")}

#####

# STEP 7 - GET ALL SPECIES-LEVEL PREDICTIONS IN A SINGLE FILE
##########################################################################################################################
# STEP 7 - GET ALL SPECIES-LEVEL PREDICTIONS IN A SINGLE FILE

# If all above code was ran for all four vertebrates groups.
# Read all species-level predictions and merge into a single file:
myfiles<-paste("ModelAveraging/SpeciesLevelPredictions_*.csv", sep="") 
SpeciesLevelPredictions<-lapply(Sys.glob(myfiles), fread, h=T)
SpeciesLevelPredictions<-rbindlist(SpeciesLevelPredictions)
fwrite(SpeciesLevelPredictions, "ModelAveraging/SpeciesLevelPredictions.csv")
