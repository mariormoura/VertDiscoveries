############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com


# SUPPORTING SCRIPT 6: ROBUSTNESS TO SAMPLING CONDITION
###########################################################################################################################
# SUPPORTING SCRIPT 6: ROBUSTNESS TO SAMPLING CONDITION

# Steps in this script:
# Species-level sensitivity analysis
#  1. Identify the best error distribution across increasing levels of right-truncation.
#  2. Run the model averaging procedure using species described across increasing levels of right-truncation.
#  3. Get the average weighted coefficients across increasing levels of right-truncation.
#  4. Get species-level predictions for all aft models and across increasing levels of right-truncation.
#  5. Get species-level average weighted discovery metrics across increasing levels of right-truncation.

# First, clean workspace:
rm(list=ls()); gc()

# Install and load R packages needed to run the analysis:
needed_packages<-c("survival", "flexsurv", "MuMIn", "usdm", "plyr", "data.table", "foreach", "doParallel", "dplyr")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages,new.packages)

#####

# STEP 1 - IDENTIFY THE BEST ERROR DISTRIBUTION ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
##########################################################################################################################
# STEP 1 - IDENTIFY THE BEST ERROR DISTRIBUTION ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION

# Clean the workspace and set the working directory:
rm(list=ls())
setwd("DefineYourDirectory")

# Create a vector with the percentanges of known-species described:
percentage_described<-seq(from=0.5, to=1, by=0.05)

# Create a dataframe to store the best error distributions identified for each vertebrate group and level of right-truncation:
Outputs<-as.data.frame(matrix(nrow=1, ncol=3))  
colnames(Outputs)<-(c("Class","PercSppDiscarded", "ErrorDistribution_1st"))

# Run the Accelerated failure time (AFT) models:
for (v in 1:4){ # v-number of vertebrate groups

  # Load the response and predictor variables:
  trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  
  # Filter the species data to comprise those described between 1759 and 2014:
  trait_data<-trait_data[trait_data$YearOfDescription>=1759 & trait_data$YearOfDescription<=2014,]
  
  # Get the response variable (time-to-event) standardized between 0 and 1:
  trait_data$Time<-(trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
  trait_data$Censor<-1 # Censor variable (it informs if the event happened)
  
  # Reorder levels of taxonomic class (just for organization purposes):
  trait_data$Class<-factor(trait_data$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
  
  # Filter the dataset to represent only one vertebrate group:
  trait_data<-trait_data[trait_data$Class==levels(trait_data$Class)[v],]
  trait_data<-droplevels(trait_data) # drop unused levels
  trait_data<-trait_data[complete.cases(trait_data), ] # remove species without data available
  
  # Create a dataframe to store the best error distribution for AFT models:
  PercDiscarded_output<-as.data.frame(matrix(nrow=length(percentage_described), ncol=3))  
  colnames(PercDiscarded_output)<-(c("Class","PercSppDiscarded", "ErrorDistribution_1st"))
  
  # For each level of right-truncation, run AFT null model (TimToEvente ~ 1) using different error distribution and extrac the BIC:
  for(j in 1:length(percentage_described)){
    
    # Filter the dataset according to the right-truncation level:
    year_threshold<-quantile(trait_data$YearOfDescription, probs=percentage_described[j], na.rm=T)
    trait_data_subset<-trait_data[trait_data$YearOfDescription<=year_threshold,]
    trait_data_subset<-droplevels(trait_data_subset)
    
    # Create an empty data.frame to store the model performance outputs for each error distribution:
    Model_fit<-as.data.frame(matrix(ncol=3, nrow=6))
    names(Model_fit)<-(c("Distrib", "npars", "BIC"))
    Model_fit[,1]<-c("exp", "weibull", "lnorm", "llogis", "gamma", "gompertz")
    
    # Run the null model (~1) for each error distribution and extract the BIC:
    Model_fit[1,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="exp", data=trait_data_subset))$npars # get the numbers of parameters
    Model_fit[1,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="exp", data=trait_data_subset)) # the the BIC of the model
    Model_fit[2,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="weibull", data=trait_data_subset))$npars
    Model_fit[2,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="weibull", data=trait_data_subset))
    Model_fit[3,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="lnorm", data=trait_data_subset))$npars
    Model_fit[3,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="lnorm", data=trait_data_subset))
    Model_fit[4,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="llogis", data=trait_data_subset))$npars
    Model_fit[4,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="llogis", data=trait_data_subset))
    Model_fit[5,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="gamma", data=trait_data_subset))$npars
    Model_fit[5,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="gamma", data=trait_data_subset))
    Model_fit[6,2]<-(flexsurvreg(Surv(Time, Censor)~1, dist="gompertz", data=trait_data_subset))$npars
    Model_fit[6,3]<-BIC(flexsurvreg(Surv(Time, Censor)~1, dist="gompertz", data=trait_data_subset))
    
    # Verify which error distribution best fits the data:
    Model_fit$deltaBIC<-Model_fit[,3]-min(Model_fit$BIC, na.rm=T)
    Model_fit$wBIC<-Weights(Model_fit$BIC)
    Model_fit<-Model_fit[order(Model_fit$BIC, decreasing=F),]
    
    # Store the results
    PercDiscarded_output[j,1]<-paste0(unique(trait_data_subset$Class)) # taxonomic class
    PercDiscarded_output[j,2]<-(1-percentage_described[j]) # percentage of recent-described spp. discarded
    PercDiscarded_output[j,3]<-Model_fit[1,1] # the best error distribution selected
    print(j)
    
  } # end of j for loop across levels of right-truncation
  
  Outputs<-rbind(Outputs, PercDiscarded_output)
  
} # end of i for loop across vertebrate groups

# Remove the first empty row:
FinalOutputs<-Outputs[-1,]
FinalOutputs$ErrorDistribution_1st<-as.factor(FinalOutputs$ErrorDistribution_1st)
FinalOutputs$Class<-as.factor(FinalOutputs$Class)

# Check for the most best selected error distribution for the full dataset:
FinalOutputs %>%
  dplyr::filter(PercSppDiscarded==0) %>%
  group_by(Class)

# Check for the most commonly selected error distribution across different levels of right-truncation:
FinalOutputs %>%
  dplyr::group_by(Class, ErrorDistribution_1st) %>% 
  dplyr::summarise(ErrorDistrib_1st=length(ErrorDistribution_1st))

# Conclusion: when using the full dataset, the Gompertz family was identifed for amphibians, reptiles, and mammals, and the Weibull for birds.
# For bird data, Gompertz family also emerged as the best error distributuon if right-truncation is added to the dataset.
# To allow cross group comparisons, we kept Gompertz family for all groups.

#####

# STEP 2 - RUN THE MODEL AVERAGING PROCEDURE USING INCREASING LEVELS OF RIGHT-TRUNCATION
##########################################################################################################################
# STEP 2 - RUN THE MODEL AVERAGING PROCEDURE USING INCREASING LEVELS OF RIGHT-TRUNCATION
rm(list=ls())

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
rm(model_combinations, my_predictors)

# Subset some models to run locally (e.g. the univariate models and the full model):
#formlist<-formlist[c(1:11,2047)] # please, omit this line to run the model averaging for all predictor combinations
formlist<-formlist[c(1,2046,2047)] # please, omit this line to run the model averaging for all predictor combinations

# Create a dataframe to store the coefficients, standard errors, and BIC of each AFT model:
list_set_of_models<-list()
list_set_of_models_L95<-list()
list_set_of_models_U95<-list()

# Create a vector with the percentanges of known-species described:
percentage_described<-seq(from=0.5, to=1, by=0.05)

# Extract the BIC, and standardized coefficents for each AFT model:
for(v in 1:4){ # v-number of vertebrate groups
  
  # Load the response and predictor variables:
  trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  trait_data<-trait_data[trait_data$YearOfDescription>=1759 & trait_data$YearOfDescription<=2014,]
  
  # Get the response variable (time-to-event) standardized between 0 and 1:
  trait_data$Time<-(trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
  trait_data$Censor<-1 # Censor variable (it informs if the event happened)
  
  # Reorder levels of taxonomic class (just for organization purposes):
  trait_data$Class<-factor(trait_data$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
  
  # Filter the dataset to represent only one vertebrate group:
  trait_data<-trait_data[trait_data$Class==levels(trait_data$Class)[v],]
  trait_data<-droplevels(trait_data) # drop unused levels
  trait_data<-trait_data[complete.cases(trait_data), ] # remove species without data available
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
  
  # Run the AFT models using species dataset with increasing levels of right-truncation:
  for(j in 1:length(percentage_described)){
    
    # Filter the dataset according to the right-truncation level:
    year_threshold<-quantile(trait_data$YearOfDescription, probs=percentage_described[j], na.rm=T)
    trait_data_subset<-trait_data[trait_data$YearOfDescription<=year_threshold,]
    trait_data_subset<-droplevels(trait_data_subset)
    
    # Create dataframes to store the coefficients, standard errors, and BIC of each AFT model:
    set_of_models<-as.data.frame(matrix(nrow=length(formlist), ncol=30))  
    set_of_models_L95<-as.data.frame(matrix(nrow=length(formlist), ncol=17))  
    set_of_models_U95<-as.data.frame(matrix(nrow=length(formlist), ncol=17))  
    colnames(set_of_models)<-(c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",
                                "LogTS", "LogPS", "LogElevM", "LogPopD", "LogTaxPerBiome", "LogRarity",
                                "LogBodySize_SE", "LogTaxPerClade_SE", "LogRangeSize_SE", "LogAMT_SE", "LogAPP_SE",
                                "LogTS_SE", "LogPS_SE", "LogElevM_SE", "LogPopD_SE", "LogTaxPerBiome_SE", "LogRarity_SE",
                                "Parameter1", "Parameter1_SE", "Parameter2", "Parameter2_SE", "FamilyError", "BIC", "PercSppDiscarded", "Class"))
    colnames(set_of_models_L95)<-c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",  "LogTS", "LogPS", "LogElevM",
                                   "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2", "FamilyError", "BIC", "PercSppDiscarded", "Class")
    colnames(set_of_models_U95)<-c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",  "LogTS", "LogPS", "LogElevM",
                                   "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2", "FamilyError", "BIC", "PercSppDiscarded", "Class")
    
    # Store id on taxa and level of right-truncation under test:
    set_of_models[,29]<-(1-percentage_described[j])
    set_of_models_L95[,16]<-(1-percentage_described[j])
    set_of_models_U95[,16]<-(1-percentage_described[j])
    set_of_models[,30]<-paste0(levels(trait_data$Class))
    set_of_models_L95[,17]<-paste0(levels(trait_data$Class))
    set_of_models_U95[,17]<-paste0(levels(trait_data$Class))
    
    # Define the best family error distribution as identifed in the STEP 1 of this script:
    bestfamily<-"gompertz"
    
    # Run AFT models for all possible predictor combinations:
    for(i in 1:length(formlist)){
      
      fit<-flexsurvreg(formlist[[i]], dist=bestfamily, data=trait_data_subset) # build the AFT model
      set_of_models[i,27]<-fit$dlist$name # store the family error distribution (the same for all models)
      set_of_models_L95[i,14]<-fit$dlist$name # store the family error distribution (the same for all models)
      set_of_models_U95[i,14]<-fit$dlist$name # store the family error distribution (the same for all models)
      set_of_models[i,28]<-BIC(fit) # store the BIC for the model i
      set_of_models_L95[i,15]<-BIC(fit) # store the BIC for the model i
      set_of_models_U95[i,15]<-BIC(fit) # store the BIC for the model i
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
    
    if(v==1){
      list_set_of_models[[j]]<-set_of_models
      list_set_of_models_L95[[j]]<-set_of_models_L95
      list_set_of_models_U95[[j]]<-set_of_models_U95
    } # store results for amphibians
    if(v==2){
      list_set_of_models[[(j+23)]]<-set_of_models
      list_set_of_models_L95[[(j+23)]]<-set_of_models_L95
      list_set_of_models_U95[[(j+23)]]<-set_of_models_U95
    } # store results for reptiles
    if(v==3){
      list_set_of_models[[(j+46)]]<-set_of_models
      list_set_of_models_L95[[(j+46)]]<-set_of_models_L95
      list_set_of_models_U95[[(j+46)]]<-set_of_models_U95
    } # store results for mammals
    if(v==4){
      list_set_of_models[[(j+69)]]<-set_of_models
      list_set_of_models_L95[[(j+69)]]<-set_of_models_L95
      list_set_of_models_U95[[(j+69)]]<-set_of_models_U95
    } # store results for birds
    
    rm(set_of_models, set_of_models_L95, set_of_models_U95, fit, model_coefs)
    
  } # end of j for loop (across levels of right-truncation)

} # end of v for loop (across vertebrate groups)

# Unlist results into a single dataframe and reorder levels of taxonomic class:
list_set_of_models<-data.table::rbindlist(list_set_of_models)
list_set_of_models_L95<-data.table::rbindlist(list_set_of_models_L95)
list_set_of_models_U95<-data.table::rbindlist(list_set_of_models_U95)
list_set_of_models$Class<-as.factor(list_set_of_models$Class)
list_set_of_models_L95$Class<-as.factor(list_set_of_models_L95$Class)
list_set_of_models_U95$Class<-as.factor(list_set_of_models_U95$Class)
list_set_of_models$Class<-factor(list_set_of_models$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
list_set_of_models_L95$Class<-factor(list_set_of_models_L95$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
list_set_of_models_U95$Class<-factor(list_set_of_models_U95$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
rm(trait_data, trait_data_subset, formlist, bestfamily) # remove unnecessary objects

# Compute the wBIC for each level of right-truncation and vertebrate group:
FinalOutput<-list()
FinalOutput_L95<-list()
FinalOutput_U95<-list()
for(v in 1:4){
  
  # Subset one vertebrate group:
  output_data1<-as.data.frame(list_set_of_models[list_set_of_models$Class==levels(list_set_of_models$Class)[v],])
  output_data1_L95<-as.data.frame(list_set_of_models_L95[list_set_of_models_L95$Class==levels(list_set_of_models_L95$Class)[v],])
  output_data1_U95<-as.data.frame(list_set_of_models_U95[list_set_of_models_U95$Class==levels(list_set_of_models_U95$Class)[v],])
  
  list1<-list()
  list2<-list()
  list3<-list()
  
  for(j in 1:length(percentage_described)){
    
    # Subset one level of right-truncation:
    output_data2<-output_data1[output_data1$PercSppDiscarded==(1-percentage_described[j]),]
    output_data2_L95<-output_data1_L95[output_data1_L95$PercSppDiscarded==(1-percentage_described[j]),]
    output_data2_U95<-output_data1_U95[output_data1_U95$PercSppDiscarded==(1-percentage_described[j]),]
    
    # Compute the wBIC:
    output_data2$wBIC<-Weights(output_data2$BIC) # compute the wBIC
    output_data2_L95$wBIC<-Weights(output_data2_L95$BIC) # compute the wBIC
    output_data2_U95$wBIC<-Weights(output_data2_U95$BIC) # compute the wBIC
    
    # Store the results
    list1[[j]]<-output_data2
    list2[[j]]<-output_data2_L95
    list3[[j]]<-output_data2_U95
    
  } # end of j for loop across levels of right-truncation
  
  rm(output_data1, output_data1_L95, output_data1_U95, output_data2, output_data2_L95, output_data2_U95)
  
  # Store the results
  FinalOutput[[v]]<-data.table::rbindlist(list1)
  FinalOutput_L95[[v]]<-data.table::rbindlist(list2)
  FinalOutput_U95[[v]]<-data.table::rbindlist(list3)
  rm(list1, list2, list3)
  
} # end of v for loop across vertebrate groups
rm(list_set_of_models, list_set_of_models_L95, list_set_of_models_U95) # remove unnecessary objects

# Unlist the results:
FinalOutput<-data.table::rbindlist(FinalOutput)
FinalOutput_L95<-data.table::rbindlist(FinalOutput_L95)
FinalOutput_U95<-data.table::rbindlist(FinalOutput_U95)

# Export the results:
dir.create("SensitivityAnalysis", showWarnings=F)
dir.create("SensitivityAnalysis/RobustnessToSamplingCondition", showWarnings=F)
fwrite(FinalOutput, "SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels.csv")
fwrite(FinalOutput_L95, "SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels_L95.csv")
fwrite(FinalOutput_U95, "SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels_U95.csv")

#####

# STEP 3. GET THE AVERAGE WEIGHTED COEFFICIENTS ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
##########################################################################################################################
# STEP 3. GET THE AVERAGE WEIGHTED COEFFICIENTS ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
rm(list=ls())

# Re-read files (if re-running particular code sections):
FinalOutput<-read.csv("SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels.csv")
FinalOutput$Class<-as.factor(FinalOutput$Class)
FinalOutput$Class<-factor(FinalOutput$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Create a vector with the percentanges of known-species described:
percentage_described<-seq(from=0.5, to=1, by=0.05)

# Get the average weighted coefficients:
AvgWeightedOutput<-list()
for (v in 1:4){
  
  # Filter one vertebrate group:
  FinalOutput_vert<-FinalOutput[FinalOutput$Class==levels(FinalOutput$Class)[v],]
  list1<-list()
  
  for (t in 1:length(percentage_described)){
    
    # Filter one level of right-truncation:
    set_of_models<-FinalOutput_vert[FinalOutput_vert$PercSppDiscarded==paste0(1-percentage_described[t]),]
    
    # Create empty dataframes to store the values:
    avg_model<-as.data.frame(matrix(nrow=5, ncol=11))
    colnames(avg_model)<-(c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP", 
                            "LogTS", "LogPS", "LogElevM", "LogPopD", "LogTaxPerBiome", "LogRarity"))
    row.names(avg_model)<-c("Avg.coef","Std.error", "Lower_IC", "Upper_IC", "Label")
    
    # Get the averaged model coefficients based on wBIC:
    for(i in 1:11){ # averaged coefficients
      candidate_models<-set_of_models[which(set_of_models[,i]!=0),]  # select all models for a given predictor
      candidate_models$wBIC<-Weights(candidate_models$BIC) # rescale the wAICc
      avg_model[1,i]<-sum(candidate_models[,i]*candidate_models$wBIC) # standardized coefficients of the weighted model
    } # averaged coefficients
    
    for(i in 12:22){ # unconditional variance estimator for the averaged model coefficients
      candidate_models<-set_of_models[which(set_of_models[,i]!=0),] # select all models for a given predictor
      weighted_var<-as.data.frame(matrix(nrow=nrow(candidate_models), ncol=1)) 
      
      for(j in 1:nrow(candidate_models)){ # j models, varying from 1 to 2^(n_covariates - 1) 
        delta_coef<-(candidate_models[j,(i-11)] - avg_model[1,(i-11)])^2 # error due to model selection uncertainty
        var_error<-(candidate_models[j,i])^2 # error in parameter estimation of the model 'j'
        weighted_var[j,1]<-(sqrt(var_error + delta_coef))*candidate_models$wBIC[j]
      }
      
      avg_model[2,(i-11)]<-sum(weighted_var[,1]) # unconditional standard error
    } # unconditional variance estimator for the averaged model coefficients
    
    for(i in 1:11){ # unconditional confidence intervals
      avg_model[3,i]<-avg_model[1,i]-(1.96*avg_model[2,i]) # lower confidence interval
      avg_model[4,i]<-avg_model[1,i]+(1.96*avg_model[2,i]) # upper confidence interval
    } # unconditional confidence intervals
    
    # Prepare the data.frame for posterior visualization and saved it:
    avg_model<-as.data.frame(t(avg_model))
    avg_model$Label<-as.factor(row.names(avg_model))
    avg_model$PercSppDiscarded<-paste0(1-percentage_described[t])
    avg_model$Class<-levels(FinalOutput$Class)[v]
   
    # Store the avg. coefs. values:
    list1[[t]]<-avg_model
  }
  
  # Store the coefs values for each vertebrate group in a list:
  AvgWeightedOutput[[v]]<-data.table::rbindlist(list1)
  
}
AvgWeightedOutput<-data.table::rbindlist(AvgWeightedOutput)

# Export the average weighted coefficients:
write.csv(AvgWeightedOutput,"SensitivityAnalysis/RobustnessToSamplingCondition/AverageWeightedCoefsAcrossTruncationLevels.csv", row.names=F, fileEncoding="UTF-8")

#####

# STEP 4 - GET SPECIES-LEVEL PREDICTIONS FOR ALL AFT MODELS AND ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
##########################################################################################################################
# STEP 4 - GET SPECIES-LEVEL PREDICTIONS FOR ALL AFT MODELS AND ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
rm(list=ls())

# Load the table containg all model coefs:
ModelParameters<-read.csv("SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels.csv")
ModelParameters_L95<-read.csv("SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels_L95.csv")
ModelParameters_U95<-read.csv("SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels_U95.csv")

# Replace NA values by 0 in the 'ModelParameters' objects:
ModelParameters[is.na(ModelParameters)] <- 0
ModelParameters_L95[is.na(ModelParameters_L95)] <- 0
ModelParameters_U95[is.na(ModelParameters_U95)] <- 0

# Define an object with percentile subdivisions (to be used to latter):
pct<-seq(0, 1, by=0.004) # set the by argument to 0.01 to reduce computational cost
# The 260 years of discovery exploration will not be accurately reflected if pct has much less than 260 units

# Create a vector with the percentanges of known-species described:
percentage_described<-seq(from=0.5, to=1, by=0.05)

# Get the species-level discovery metrics for each AFT model and level of right-truncation:
a<-Sys.time()
for(v in 1:4){ # v-number of vertebrate groups
  
  # Load the response and predictor variables:
  raw_trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  raw_trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>=1759 & raw_trait_data$YearOfDescription<=2014,]
  
  # Get the response variable (time-to-event) standardized between 0 and 1:
  raw_trait_data$Time<-(raw_trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
  raw_trait_data$Censor<-1 # Censor variable (it informs if the event happened)
  
  # Reorder levels of taxonomic class (just for organization purposes):
  raw_trait_data$Class<-factor(raw_trait_data$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
  
  # Filter the dataset to represent only one vertebrate group:
  raw_trait_data<-raw_trait_data[raw_trait_data$Class==levels(raw_trait_data$Class)[v],]
  raw_trait_data<-droplevels(raw_trait_data) # drop unused levels
  raw_trait_data<-raw_trait_data[complete.cases(raw_trait_data), ] # remove species without data available
  raw_trait_data$LogBodySize<-scale(log10(raw_trait_data$BodySize), center=T, scale=T)
  raw_trait_data$LogRangeSize<-scale(log10(raw_trait_data$RangeSize), center=T, scale=T)
  raw_trait_data$LogRarity<-scale(log10(raw_trait_data$RangeRarity+1), center=T, scale=T)
  raw_trait_data$LogAPP<-scale(log10(raw_trait_data$AnnuPrecip+1), center=T, scale=T)
  raw_trait_data$LogAMT<-scale(log10((raw_trait_data$AnnuMeanTemp+2730)/10), center=T, scale=T) # convert from Kelvin to Celsius degree
  raw_trait_data$LogTS<-scale(log10(raw_trait_data$TempSeasonality+1), center=T, scale=T)
  raw_trait_data$LogPS<-scale(log10(raw_trait_data$PrecipSeasonality+1), center=T, scale=T)
  raw_trait_data$LogElevM<-scale(log10(raw_trait_data$Elevation+1), center=T, scale=T)
  raw_trait_data$LogPopD<-scale(log10(raw_trait_data$HumanDensity+1), center=T, scale=T)
  raw_trait_data$LogTaxPerClade<-scale(log10(raw_trait_data$ActivityFamily+1), center=T, scale=T)
  raw_trait_data$LogTaxPerBiome<-scale(log10(raw_trait_data$ActivityBioregion+1), center=T, scale=T)
  
  # Create empty lists to store some outputs computed across increasing levels of right-truncation (see below):
  list1<-list() # this will receive the object DiscProb_dt
  list2<-list() # this will receive the object DiscProb_dt_L95
  list3<-list() # this will receive the object DiscProb_dt_U95
  list4<-list() # this will receive the object YearPred_dt
  
  for(j in 1:length(percentage_described)){
    
    # Filter the dataset according to the right-truncation level:
    year_threshold<-quantile(raw_trait_data$YearOfDescription, probs=percentage_described[j], na.rm=T)
    trait_data<-raw_trait_data[raw_trait_data$YearOfDescription<=year_threshold,]
    trait_data<-droplevels(trait_data)
    
    # Filter the ModelParameters object to include one vertebrate group and one level of right-truncation:
    set_of_models<-ModelParameters[ModelParameters$Class==levels(raw_trait_data$Class) & ModelParameters$PercSppDiscarded==paste0(1-percentage_described[j]),]
    set_of_models_L95<-ModelParameters_L95[ModelParameters_L95$Class==levels(raw_trait_data$Class) & ModelParameters$PercSppDiscarded==paste0(1-percentage_described[j]),]
    set_of_models_U95<-ModelParameters_U95[ModelParameters_U95$Class==levels(raw_trait_data$Class) & ModelParameters$PercSppDiscarded==paste0(1-percentage_described[j]),]
    
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
        
        # colnames = species & rownames = discovery probability values
        names(TimeToDisc_per_spp)<-trait_data$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
        
        # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
        names(ProbDisc_per_spp)<-trait_data$Binomial; rownames(ProbDisc_per_spp)<-pct 
        names(ProbDisc_per_spp_L95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
        names(ProbDisc_per_spp_U95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
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
        
        # colnames = species & rownames = discovery probability values
        names(TimeToDisc_per_spp)<-trait_data$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
        
        # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
        names(ProbDisc_per_spp)<-trait_data$Binomial; rownames(ProbDisc_per_spp)<-pct 
        names(ProbDisc_per_spp_L95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
        names(ProbDisc_per_spp_U95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
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
        
        # colnames = species & rownames = discovery probability values
        names(TimeToDisc_per_spp)<-trait_data$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
        
        # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
        names(ProbDisc_per_spp)<-trait_data$Binomial; rownames(ProbDisc_per_spp)<-pct 
        names(ProbDisc_per_spp_L95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
        names(ProbDisc_per_spp_U95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
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
        
        # colnames = species & rownames = discovery probability values
        names(TimeToDisc_per_spp)<-trait_data$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
        
        # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
        names(ProbDisc_per_spp)<-trait_data$Binomial; rownames(ProbDisc_per_spp)<-pct 
        names(ProbDisc_per_spp_L95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
        names(ProbDisc_per_spp_U95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
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
      names(TimeToDisc_per_spp)<-trait_data$Binomial; rownames(TimeToDisc_per_spp)<-pct # 
      
      # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
      names(ProbDisc_per_spp)<-trait_data$Binomial; rownames(ProbDisc_per_spp)<-pct 
      names(ProbDisc_per_spp_L95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
      names(ProbDisc_per_spp_U95)<-trait_data$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
      
      # Extract the discovery probability in 2015 for all species: 
      present_year<-(2015-1758)/(260) # get the year 2015 in a scale of 0 to 1 (0 = 1758, 1 = 2018)
      row_selected<-which(row.names(ProbDisc_per_spp)>=present_year)[1] # select the row of the ProbDisc_per_spp object the represents the year 2015
      ProbDisc_2015<-t(ProbDisc_per_spp[row_selected,]) # extract the discovery probability in 2015 according to the model i
      ProbDisc_2015_L95<-t(ProbDisc_per_spp_L95[row_selected,]) # extract the lower bound of the discovery probability in 2015 according to the model i
      ProbDisc_2015_U95<-t(ProbDisc_per_spp_U95[row_selected,]) # extract the upper bound of the discovery probability in 2015 according to the model i
      
      # Extract the predicted year for the discovery probability of 0.50:
      TimeToDisc_per_spp<-(TimeToDisc_per_spp*(2018-1758))+1758 # convert row values to year
      threshold_selected<-which(row.names(TimeToDisc_per_spp)>=0.50) # select the rows with discovery probability >=0.50
      thershold_selected<-threshold_selected[1] # select the first row informing the discovery probability >=0.50
      YearPred<-t(TimeToDisc_per_spp[threshold_selected, ]) # year predicted for a given species description (99.5% of probability)
      
      # Update the 'DiscProb_dt' and 'YearPred_dt' objects with the discovery probability and estimated year of discovery for each species (k) 
      # according to the model i (= rows of the data.table) and column j (= columns of the data.table = species):
      for (k in 1:nrow(trait_data)){data.table::set(DiscProb_dt, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015[k])}
      for (k in 1:nrow(trait_data)){data.table::set(YearPred_dt, i=as.integer(i), j=as.integer(k), value=YearPred[k])}
      for (k in 1:nrow(trait_data)){data.table::set(DiscProb_dt_L95, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015_L95[k])}
      for (k in 1:nrow(trait_data)){data.table::set(DiscProb_dt_U95, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015_U95[k])}
      
    } # end of i for loop (across all AFT models)
    stopCluster(cl)
    
    # Register the level of right-truncation (% of spp. discarded):
    DiscProb_dt$PercSppDiscarded<-paste0(1-percentage_described[j])
    DiscProb_dt_L95$PercSppDiscarded<-paste0(1-percentage_described[j])
    DiscProb_dt_U95$PercSppDiscarded<-paste0(1-percentage_described[j])
    YearPred_dt$PercSppDiscarded<-paste0(1-percentage_described[j])
    
    # Store the results in the empty lists prepared a priori:
    list1[[j]]<-DiscProb_dt
    list2[[j]]<-DiscProb_dt_L95
    list3[[j]]<-DiscProb_dt_U95
    list4[[j]]<-YearPred_dt
    
    }

  # Export the outputs
  if(v==1){
    save(list1, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_Amphibia.RData")
    save(list2, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_L95_Amphibia.RData")
    save(list3, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_U95_Amphibia.RData")
    save(list4, file="SensitivityAnalysis/RobustnessToSamplingCondition/YearPred_dt_Amphibia.RData")
  } # store estimates for amphibians
  if(v==2){
    save(list1, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_Reptilia.RData")
    save(list2, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_L95_Reptilia.RData")
    save(list3, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_U95_Reptilia.RData")
    save(list4, file="SensitivityAnalysis/RobustnessToSamplingCondition/YearPred_dt_Reptilia.RData")
  } # store estimates for reptiles
  if(v==3){
    save(list1, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_Mammalia.RData")
    save(list2, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_L95_Mammalia.RData")
    save(list3, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_U95_Mammalia.RData")
    save(list4, file="SensitivityAnalysis/RobustnessToSamplingCondition/YearPred_dt_Mammalia.RData")
  } # store estimates for mammals
  if(v==4){
    save(list1, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_Aves.RData")
    save(list2, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_L95_Aves.RData")
    save(list3, file="SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_U95_Aves.RData")
    save(list4, file="SensitivityAnalysis/RobustnessToSamplingCondition/YearPred_dt_Aves.RData")
  } # store estimates for birds
  
  } # end of v for loop across vertebrate groups
b<-Sys.time() # how long it took
b-a # about 1 hours for a subset of only 3 models running in a notebook with 2.11 Ghz, 24GB RAM, 8 cores.

#####

# STEP 5 - GET SPECIES-LEVEL AVERAGE WEIGHTED DISCOVERY METRICS ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
##########################################################################################################################
# STEP 5 - GET SPECIES-LEVEL AVERAGE WEIGHTED DISCOVERY METRICS ACROSS INCREASING LEVELS OF RIGHT-TRUNCATION
rm(list=ls())

# Function to load RData using user-specified object names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create a vector with the percentanges of known-species described:
percentage_described<-seq(from=0.5, to=1, by=0.05)
VertGroup<-c("Amphibia", "Reptilia", "Mammalia", "Aves")

# Create empty lists to store the output:
AmphList<-list() # to store output on amphibians
ReptList<-list() # to store output on reptiles
MammList<-list() # to store output on mammals
AvesList<-list() # to store output on birds

# Compute the average weighted discovery metrics for each level of right-truncation and vertebrate group:
for(v in 1:4){ # v-number of vertebrate groups
  
  # Load the response and predictor variables:
  raw_trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  raw_trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>=1759 & raw_trait_data$YearOfDescription<=2014,]
  
  # Get the response variable (time-to-event) standardized between 0 and 1:
  raw_trait_data$Time<-(raw_trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
  raw_trait_data$Censor<-1 # Censor variable (it informs if the event happened)
  
  # Filter the dataset to represent only one vertebrate group:
  raw_trait_data<-raw_trait_data[raw_trait_data$Class==VertGroup[v],]
  raw_trait_data<-droplevels(raw_trait_data) # drop unused levels
  raw_trait_data<-raw_trait_data[complete.cases(raw_trait_data), ] # remove species without data available
  raw_trait_data<-raw_trait_data[, .(Binomial, Class, YearOfDescription)]
  
  for(j in 1:length(percentage_described)){
    
    # Filter the dataset according to the right-truncation level:
    year_threshold<-quantile(raw_trait_data$YearOfDescription, probs=percentage_described[j], na.rm=T)
    trait_data<-raw_trait_data[raw_trait_data$YearOfDescription<=year_threshold,]
    trait_data<-droplevels(trait_data)
    
    # Read the data.table containg the BIC of each AFT model and filter the vertebrate group v and level of right-truncation t:
    set_of_models<-as.data.frame(fread("SensitivityAnalysis/RobustnessToSamplingCondition/ModelParametersAcrossTruncationLevels.csv", stringsAsFactors=TRUE, encoding="UTF-8"))
    set_of_models<-set_of_models[set_of_models$Class==VertGroup[v] & set_of_models$PercSppDiscarded==paste0(1-percentage_described[j]),]
    
    # Load the RData object for the vertebrate group v and level of right-truncation t:
    DiscProb_dt<-loadRData(file.path(paste0("SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_", VertGroup[v], ".RData")))[[j]]
    DiscProb_dt_L95<-loadRData(file.path(paste0("SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_L95_", VertGroup[v], ".RData")))[[j]]
    DiscProb_dt_U95<-loadRData(file.path(paste0("SensitivityAnalysis/RobustnessToSamplingCondition/DiscProb_dt_U95_", VertGroup[v], ".RData")))[[j]]
    YearPred_dt<-loadRData(file.path(paste0("SensitivityAnalysis/RobustnessToSamplingCondition/YearPred_dt_", VertGroup[v], ".RData")))[[j]]
    
    # Keep only the columns representing species in the dataset:
    cols<-names(DiscProb_dt); cols<-cols[-length(cols)] # the last column informs TimePeriod
    DiscProb_dt<-DiscProb_dt[, ..cols]
    DiscProb_dt_L95<-DiscProb_dt_L95[, ..cols]
    DiscProb_dt_U95<-DiscProb_dt_U95[, ..cols]
    YearPred_dt<-YearPred_dt[, ..cols]
    
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
    trait_data$PercSppDiscarded<-paste0(1-percentage_described[j])
    
    # Store the outputs for each vertebrate group:
    if(v==1){AmphList[[j]]<-trait_data}
    if(v==2){ReptList[[j]]<-trait_data}
    if(v==3){MammList[[j]]<-trait_data}
    if(v==4){AvesList[[j]]<-trait_data}
    
    rm(AvgW_DiscProb_dt, AvgW_DiscProb_dt_L95, AvgW_DiscProb_dt_U95, AvgW_YearPred_dt, trait_data,
       DiscProb_dt, DiscProb_dt_L95, DiscProb_dt_U95, YearPred_dt, set_of_models, cols, spp_columns)
    
  } # end of t for loop across levels of right-truncation
  
} # end of v for loop across vertebrate groups

# Export the results:
save(AmphList, file="SensitivityAnalysis/RobustnessToSamplingCondition/SpeciesLevelEstimatesAcrossTruncationLevels_Amphibia.RData")
save(ReptList, file="SensitivityAnalysis/RobustnessToSamplingCondition/SpeciesLevelEstimatesAcrossTruncationLevels_Reptilia.RData")
save(MammList, file="SensitivityAnalysis/RobustnessToSamplingCondition/SpeciesLevelEstimatesAcrossTruncationLevels_Mammalia.RData")
save(AvesList, file="SensitivityAnalysis/RobustnessToSamplingCondition/SpeciesLevelEstimatesAcrossTruncationLevels_Aves.RData")

#####