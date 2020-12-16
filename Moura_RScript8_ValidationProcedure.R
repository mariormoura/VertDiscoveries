############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com


# SUPPORTING SCRIPT 8: VALIDATION PROCEDURE (SUMMARY)
############################################################################################################################
# SUPPORTING SCRIPT 8: VALIDATION PROCEDURE (SUMMARY)

# Steps in this script:
# Species-level sensitivity analysis
#  1. Run AFT models across time periods, cross-validation partition sizes, and k-fold subdatasets.
#  2. Get the avg. weighted coefficients across time periods, cross-validation partition sizes, and k-fold subdatasets.
#  4. Get species-level predictions across AFT models, time periods, cross-validation partition sizes, and k-fold subdatasets.
#  5. Get species-level average weighted discovery metrics across different time periods.
#  6. Taxon-level analysis, compute the per-taxon discovery metrics across different time periods.
#  7. Assemblage-level analysis, compute the per-assemblage discovery metrics across different time periods.

# First, clean workspace:
rm(list=ls()); gc()

# Install and load R packages needed to run the analysis:
needed_packages<-c("survival", "flexsurv", "MuMIn", "usdm", "plyr", "data.table", "foreach", "doParallel", "dplyr")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages,new.packages)

#####

# STEP 1 - RUN AFT MODELS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
##########################################################################################################################
# STEP 1 - RUN AFT MODELS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
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

# Subset some models to run locally (e.g. the univariate models or the full model):
formlist<-formlist[c(4,2047)] # please, omit this line to run the model averaging for all predictor combinations

# Create an operator to get the opposite of %in%:
'%ni%' <- Negate('%in%')

# Create objects to guide for loops:
TimePeriod<-c(1758, seq(from=1760, to=1970, by=10)) # the starting date of each time period
VertGroups<-c("Amphibia", "Reptilia", "Mammalia", "Aves")
PercSppModelled<-c("25T75V", "50T50V", "75T25V", "90T10V") # T = training-fold / V = validation-fold
dir.create("ModelValidation", showWarnings=F)

for(v in 1:4){ # v-number of vertebrate groups
  
  # Perform computations for each cross-validation partition:
  for(c in 1:2){ # c-cross validation partition sizes
    
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
    
    # Create a dataframe to store the outputs for each AFT model, TimePeriod, KFold_partition, VertebrateGroup:
    FullList_set_of_models<-list()
    FullList_set_of_models_L95<-list()
    FullList_set_of_models_U95<-list()
    
    # Run the AFT models using species described within different time periods:
    for(t in 1:length(TimePeriod)){
      
      # Create a dataframe to store the coefficients, standard errors, and BIC of each AFT model:
      list_set_of_models<-list()
      list_set_of_models_L95<-list()
      list_set_of_models_U95<-list()
      
      # Filter species for a given time period:
      trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>TimePeriod[t],]
      
      # Randomly split species data in training- and validation-data subsets:
      set.seed(12)
      if(c==1){ # 25% of species in the training-fold and 75% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_fold<-list()
        species_fold[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/4)],1]
        species_fold[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/4)+1):(2*ceiling(length(random_rows)/4))],1]
        species_fold[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/4)+1):(3*ceiling(length(random_rows)/4))],1]
        species_fold[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/4)+1):length(random_rows)],1]
        
      } # 25% of species in the training-fold and 75% in the validation-fold
      if(c==2){ # 50% of species in the training-fold and 50% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_fold<-list()
        species_fold[[1]]<-trait_data[random_rows[1:floor(length(random_rows)/2)],1]
        species_fold[[2]]<-trait_data[random_rows[(floor(length(random_rows)/2)+1):length(random_rows)],1]
        
      } # 50% of species in the training-fold and 50% in the validation-fold
      if(c==3){ # 75% of species in the training-fold and 25% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_to_discard<-list()
        species_fold<-list()
        species_to_discard[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/4)],1]
        species_to_discard[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/4)+1):(2*ceiling(length(random_rows)/4))],1]
        species_to_discard[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/4)+1):(3*ceiling(length(random_rows)/4))],1]
        species_to_discard[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/4)+1):length(random_rows)],1]
        species_fold[[1]]<-trait_data[Binomial %ni% as.character(species_to_discard[[1]]$Binomial)]
        species_fold[[2]]<-trait_data[Binomial %ni% as.character(species_to_discard[[2]]$Binomial)]
        species_fold[[3]]<-trait_data[Binomial %ni% as.character(species_to_discard[[3]]$Binomial)]
        species_fold[[4]]<-trait_data[Binomial %ni% as.character(species_to_discard[[4]]$Binomial)]
        rm(species_to_discard)
        
      } # 75% of species in the training-fold and 25% in the validation-fold
      if(c==4){ # 90% of species in the training-fold and 10% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_to_discard<-list()
        species_fold<-list()
        species_to_discard[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/10)],1]
        species_to_discard[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/10)+1):(2*ceiling(length(random_rows)/10))],1]
        species_to_discard[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/10)+1):(3*ceiling(length(random_rows)/10))],1]
        species_to_discard[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/10)+1):(4*ceiling(length(random_rows)/10))],1]
        species_to_discard[[5]]<-trait_data[random_rows[(4*ceiling(length(random_rows)/10)+1):(5*ceiling(length(random_rows)/10))],1]
        species_to_discard[[6]]<-trait_data[random_rows[(5*ceiling(length(random_rows)/10)+1):(6*ceiling(length(random_rows)/10))],1]
        species_to_discard[[7]]<-trait_data[random_rows[(6*ceiling(length(random_rows)/10)+1):(7*ceiling(length(random_rows)/10))],1]
        species_to_discard[[8]]<-trait_data[random_rows[(7*ceiling(length(random_rows)/10)+1):(8*ceiling(length(random_rows)/10))],1]
        species_to_discard[[9]]<-trait_data[random_rows[(8*ceiling(length(random_rows)/10)+1):(9*ceiling(length(random_rows)/10))],1]
        species_to_discard[[10]]<-trait_data[random_rows[(9*ceiling(length(random_rows)/10)+1):length(random_rows)],1]
        species_fold[[1]]<-trait_data[Binomial %ni% as.character(species_to_discard[[1]]$Binomial)]
        species_fold[[2]]<-trait_data[Binomial %ni% as.character(species_to_discard[[2]]$Binomial)]
        species_fold[[3]]<-trait_data[Binomial %ni% as.character(species_to_discard[[3]]$Binomial)]
        species_fold[[4]]<-trait_data[Binomial %ni% as.character(species_to_discard[[4]]$Binomial)]
        species_fold[[5]]<-trait_data[Binomial %ni% as.character(species_to_discard[[5]]$Binomial)]
        species_fold[[6]]<-trait_data[Binomial %ni% as.character(species_to_discard[[6]]$Binomial)]
        species_fold[[7]]<-trait_data[Binomial %ni% as.character(species_to_discard[[7]]$Binomial)]
        species_fold[[8]]<-trait_data[Binomial %ni% as.character(species_to_discard[[8]]$Binomial)]
        species_fold[[9]]<-trait_data[Binomial %ni% as.character(species_to_discard[[9]]$Binomial)]
        species_fold[[10]]<-trait_data[Binomial %ni% as.character(species_to_discard[[10]]$Binomial)]
        rm(species_to_discard)
        
      } # 90% of species in the training-fold and 10% in the validation-fold
      for(d in 1:length(species_fold)){species_fold[[d]]<-droplevels(species_fold[[d]])} # discard non-used levels:
    
      # For loop across the k-fold partitions:
      for(b in 1:length(species_fold)){
        
        # Filter the species described within the time period of investigation:
        trait_data_subset<-trait_data[Binomial %in% as.character(species_fold[[b]]$Binomial)]
        
        # Create dataframes to store the coefficients, standard errors, and BIC of each AFT model:
        set_of_models<-as.data.frame(matrix(nrow=length(formlist), ncol=28))  
        set_of_models_L95<-as.data.frame(matrix(nrow=length(formlist), ncol=15))  
        set_of_models_U95<-as.data.frame(matrix(nrow=length(formlist), ncol=15))  
        colnames(set_of_models)<-(c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",
                                    "LogTS", "LogPS", "LogElevM", "LogPopD", "LogTaxPerBiome", "LogRarity",
                                    "LogBodySize_SE", "LogTaxPerClade_SE", "LogRangeSize_SE", "LogAMT_SE", "LogAPP_SE",
                                    "LogTS_SE", "LogPS_SE", "LogElevM_SE", "LogPopD_SE", "LogTaxPerBiome_SE", "LogRarity_SE",
                                    "Parameter1", "Parameter1_SE", "Parameter2", "Parameter2_SE", "FamilyError", "BIC"))
        colnames(set_of_models_L95)<-c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",  "LogTS", "LogPS", "LogElevM",
                                       "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2", "FamilyError", "BIC")
        colnames(set_of_models_U95)<-c("LogBodySize", "LogTaxPerClade", "LogRangeSize", "LogAMT", "LogAPP",  "LogTS", "LogPS", "LogElevM",
                                       "LogPopD", "LogTaxPerBiome", "LogRarity", "Parameter1", "Parameter2", "FamilyError", "BIC")
        
        # Store id on taxa and time period under test:
        set_of_models$TimePeriod<-TimePeriod[t]
        set_of_models_L95$TimePeriod<-TimePeriod[t]
        set_of_models_U95$TimePeriod<-TimePeriod[t]
        set_of_models$CV_PartitionSize<-PercSppModelled[c]
        set_of_models_L95$CV_PartitionSize<-PercSppModelled[c]
        set_of_models_U95$CV_PartitionSize<-PercSppModelled[c]
        set_of_models$KFold<-b
        set_of_models_L95$KFold<-b
        set_of_models_U95$KFold<-b
        set_of_models$Class<-paste0(levels(trait_data_subset$Class))
        set_of_models_L95$Class<-paste0(levels(trait_data_subset$Class))
        set_of_models_U95$Class<-paste0(levels(trait_data_subset$Class))
        
        # Define the best family error distribution as identifed in other steps outside this script:
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
        
        set_of_models$wBIC<-Weights(set_of_models$BIC) # compute the wBIC
        set_of_models_L95$wBIC<-Weights(set_of_models_L95$BIC) # compute the wBIC
        set_of_models_U95$wBIC<-Weights(set_of_models_U95$BIC) # compute the wBIC
        list_set_of_models[[b]]<-set_of_models
        list_set_of_models_L95[[b]]<-set_of_models_L95
        list_set_of_models_U95[[b]]<-set_of_models_U95
        rm(set_of_models, set_of_models_L95, set_of_models_U95, fit, model_coefs, trait_data_subset)
        
      } # end of b for loop (across training-folds)
      
      # Unbind the outputs for all kfold partitions in a single data.table:
      FullList_set_of_models[[t]]<-data.table::rbindlist(list_set_of_models)
      FullList_set_of_models_L95[[t]]<-data.table::rbindlist(list_set_of_models_L95)
      FullList_set_of_models_U95[[t]]<-data.table::rbindlist(list_set_of_models_U95)
      rm(trait_data, list_set_of_models, list_set_of_models_L95, list_set_of_models_U95)
      
    } # end of j for loop (across time periods)
    
    # Unbind the outputs for all kfold partitions and time periods in a single data.table:
    
    VertList<-data.table::rbindlist(FullList_set_of_models)
    VertList_L95<-data.table::rbindlist(FullList_set_of_models_L95)
    VertList_U95<-data.table::rbindlist(FullList_set_of_models_U95)
    save(VertList, file=paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(VertList_L95, file=paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_L95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(VertList_U95, file=paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_U95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    
    } # end of c for loop (across sizes of k-fold partitions)
    
    rm(raw_trait_data, FullList_set_of_models, FullList_set_of_models_L95, FullList_set_of_models_U95)
  
  } # end of v for loop (across vertebrate groups)

#####

# STEP 2 - GET THE AVG. WEIGHTED COEFFICIENTS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
##########################################################################################################################
# STEP 2 - GET THE AVG. WEIGHTED COEFFICIENTS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
rm(list=ls())

# Function to load RData using user-specified object names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create objects to guide for loops:
TimePeriod<-c(1758, seq(from=1760, to=1970, by=10)) # the starting date of each time period
VertGroups<-c("Amphibia", "Reptilia", "Mammalia", "Aves")
PercSppModelled<-c("25T75V", "50T50V", "75T25V", "90T10V") # T = training-fold / V = validation-fold
Kfold_partitions<-c(4, 2, 4, 10) # number of k-fold partitions represented in each levels of 'PercSppModelled'

# Perform computations for each vertebrate group:
for(v in 1:4){ # v-number of vertebrate groups
  list3<-list()
  
  # Perform computations for each cross-validation partition:
  for(c in 1:2){ # c-cross validation partition sizes
    list2<-list()
    
    # Perform computations for each time period:
    for(t in 1:length(TimePeriod)){
      list1<-list()
      
      # Perform computations for each k-fold partitions:
      for(b in 1:Kfold_partitions[c]){
        
        # Load the RData object for the vertebrate group v and time period t:
        FinalOutput<-loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_", VertGroups[v], "_", PercSppModelled[c], ".RData")))
        FinalOutput<-as.data.frame(FinalOutput)
        
        # Subset one time period and one k-fold partition:
        set_of_models<-FinalOutput[FinalOutput$TimePeriod==TimePeriod[t] & FinalOutput$KFold==b,]
        
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
        avg_model$TimePeriod<-TimePeriod[t]
        avg_model$CV_PartitionSize<-PercSppModelled[c]
        avg_model$KFold<-b
        avg_model$Class<-VertGroups[v]
        
        # Store the avg. coefs. values:
        list1[[b]]<-avg_model
        
      } # enf of b for loop (k-fold partitions)
      list2[[t]]<-data.table::rbindlist(list1)
          
    } # end of t for loop (time periods)
    
    list3[[c]]<-data.table::rbindlist(list2)
    
  } # end of c for loop (cross-validation partition sizes)
  
  # Unbind the final list and export:
  AvgWeightedOutput<-data.table::rbindlist(list3)
  write.csv(AvgWeightedOutput, paste0("ModelValidation/AverageWeightedCoefsAcrossTimeAcrossKfolds_", 
                                      VertGroups[v], "_", PercSppModelled[c], ".RData"), fileEncoding="UTF-8")
  
  } # end of v for loop (vertebrate groups)

#####

# STEP 3 - GET SPECIES-LEVEL PREDICTIONS ACROSS AFT MODELS, TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
##########################################################################################################################
# STEP 3 - GET SPECIES-LEVEL PREDICTIONS ACROSS AFT MODELS, TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
rm(list=ls())

# Function to load RData using user-specified object names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create objects to guide for loops:
TimePeriod<-c(1758, seq(from=1760, to=1970, by=10)) # the starting date of each time period
VertGroups<-c("Amphibia", "Reptilia", "Mammalia", "Aves")
PercSppModelled<-c("25T75V", "50T50V", "75T25V", "90T10V") # T = training-fold / V = validation-fold
Kfold_partitions<-c(4, 2, 4, 10) # number of k-fold partitions represented in each levels of 'PercSppModelled'

# Define an object with percentile subdivisions (to be used to latter):
pct<-seq(0, 1, by=0.004) # set the by argument to 0.01 to reduce computational cost
# The 260 years of discovery exploration will not be accurately reflected if pct has much less than 260 units

# Perform computations for each vertebrate group:
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
  raw_trait_data<-raw_trait_data[raw_trait_data$Class==VertGroups[v],]
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
  
  # Perform computations for each cross-validation partition:
  for(c in 1:2){ # c-cross validation partition sizes
    
    # Create empty lists to store some outputs computed across different time periods (see below):
    list1_2<-list() # this will receive the object DiscProb_dt for each time period
    list2_2<-list() # this will receive the object DiscProb_dt_L95 for each time period
    list3_2<-list() # this will receive the object DiscProb_dt_U95 for each time period
    list4_2<-list() # this will receive the object YearPred_dt for each time period
    
    # Perform computations for each time period:
    for(t in 1:length(TimePeriod)){
      
      # Filter species for a given time period:
      trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>TimePeriod[t],]
      
      # Randomly split species data in training- and validation-data subsets:
      set.seed(12)
      if(c==1){ # 25% of species in the training-fold and 75% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_fold<-list()
        species_fold[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/4)],1]
        species_fold[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/4)+1):(2*ceiling(length(random_rows)/4))],1]
        species_fold[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/4)+1):(3*ceiling(length(random_rows)/4))],1]
        species_fold[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/4)+1):length(random_rows)],1]
        
      } # 25% of species in the training-fold and 75% in the validation-fold
      if(c==2){ # 50% of species in the training-fold and 50% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_fold<-list()
        species_fold[[1]]<-trait_data[random_rows[1:floor(length(random_rows)/2)],1]
        species_fold[[2]]<-trait_data[random_rows[(floor(length(random_rows)/2)+1):length(random_rows)],1]
        
      } # 50% of species in the training-fold and 50% in the validation-fold
      if(c==3){ # 75% of species in the training-fold and 25% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_to_discard<-list()
        species_fold<-list()
        species_to_discard[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/4)],1]
        species_to_discard[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/4)+1):(2*ceiling(length(random_rows)/4))],1]
        species_to_discard[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/4)+1):(3*ceiling(length(random_rows)/4))],1]
        species_to_discard[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/4)+1):length(random_rows)],1]
        species_fold[[1]]<-trait_data[Binomial %ni% as.character(species_to_discard[[1]]$Binomial)]
        species_fold[[2]]<-trait_data[Binomial %ni% as.character(species_to_discard[[2]]$Binomial)]
        species_fold[[3]]<-trait_data[Binomial %ni% as.character(species_to_discard[[3]]$Binomial)]
        species_fold[[4]]<-trait_data[Binomial %ni% as.character(species_to_discard[[4]]$Binomial)]
        rm(species_to_discard)
        
      } # 75% of species in the training-fold and 25% in the validation-fold
      if(c==4){ # 90% of species in the training-fold and 10% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_to_discard<-list()
        species_fold<-list()
        species_to_discard[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/10)],1]
        species_to_discard[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/10)+1):(2*ceiling(length(random_rows)/10))],1]
        species_to_discard[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/10)+1):(3*ceiling(length(random_rows)/10))],1]
        species_to_discard[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/10)+1):(4*ceiling(length(random_rows)/10))],1]
        species_to_discard[[5]]<-trait_data[random_rows[(4*ceiling(length(random_rows)/10)+1):(5*ceiling(length(random_rows)/10))],1]
        species_to_discard[[6]]<-trait_data[random_rows[(5*ceiling(length(random_rows)/10)+1):(6*ceiling(length(random_rows)/10))],1]
        species_to_discard[[7]]<-trait_data[random_rows[(6*ceiling(length(random_rows)/10)+1):(7*ceiling(length(random_rows)/10))],1]
        species_to_discard[[8]]<-trait_data[random_rows[(7*ceiling(length(random_rows)/10)+1):(8*ceiling(length(random_rows)/10))],1]
        species_to_discard[[9]]<-trait_data[random_rows[(8*ceiling(length(random_rows)/10)+1):(9*ceiling(length(random_rows)/10))],1]
        species_to_discard[[10]]<-trait_data[random_rows[(9*ceiling(length(random_rows)/10)+1):length(random_rows)],1]
        species_fold[[1]]<-trait_data[Binomial %ni% as.character(species_to_discard[[1]]$Binomial)]
        species_fold[[2]]<-trait_data[Binomial %ni% as.character(species_to_discard[[2]]$Binomial)]
        species_fold[[3]]<-trait_data[Binomial %ni% as.character(species_to_discard[[3]]$Binomial)]
        species_fold[[4]]<-trait_data[Binomial %ni% as.character(species_to_discard[[4]]$Binomial)]
        species_fold[[5]]<-trait_data[Binomial %ni% as.character(species_to_discard[[5]]$Binomial)]
        species_fold[[6]]<-trait_data[Binomial %ni% as.character(species_to_discard[[6]]$Binomial)]
        species_fold[[7]]<-trait_data[Binomial %ni% as.character(species_to_discard[[7]]$Binomial)]
        species_fold[[8]]<-trait_data[Binomial %ni% as.character(species_to_discard[[8]]$Binomial)]
        species_fold[[9]]<-trait_data[Binomial %ni% as.character(species_to_discard[[9]]$Binomial)]
        species_fold[[10]]<-trait_data[Binomial %ni% as.character(species_to_discard[[10]]$Binomial)]
        rm(species_to_discard)
        
      } # 90% of species in the training-fold and 10% in the validation-fold
      for(d in 1:length(species_fold)){species_fold[[d]]<-droplevels(species_fold[[d]])} # discard non-used levels
      
      # Create empty lists to store some outputs computed across different time periods (see below):
      list1_1<-list() # this will receive the object DiscProb_dt across k-folds
      list2_1<-list() # this will receive the object DiscProb_dt_L95 across k-folds
      list3_1<-list() # this will receive the object DiscProb_dt_U95 across k-folds
      list4_1<-list() # this will receive the object YearPred_dt across k-folds
      
      # Perform computations for each k-fold partitions:
      for(b in 1:Kfold_partitions[c]){
        
        # Filter the species described within the time period of investigation:
        trait_data_subset<-trait_data[Binomial %in% as.character(species_fold[[b]]$Binomial)]
        
        # Load the RData object for the vertebrate group v and time period t:
        FinalOutput<-loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_", VertGroups[v], "_", PercSppModelled[c], ".RData")))
        FinalOutput_L95<-loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_L95_", VertGroups[v], "_", PercSppModelled[c], ".RData")))
        FinalOutput_U95<-loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_U95_", VertGroups[v], "_", PercSppModelled[c], ".RData")))
        
        # Subset one time period and one k-fold partition:
        FinalOutput<-as.data.frame(FinalOutput)
        FinalOutput_L95<-as.data.frame(FinalOutput_L95)
        FinalOutput_U95<-as.data.frame(FinalOutput_U95)
        set_of_models<-FinalOutput[FinalOutput$TimePeriod==TimePeriod[t] & FinalOutput$KFold==b,]
        set_of_models_L95<-FinalOutput_L95[FinalOutput_L95$TimePeriod==TimePeriod[t] & FinalOutput_L95$KFold==b,]
        set_of_models_U95<-FinalOutput_U95[FinalOutput_U95$TimePeriod==TimePeriod[t] & FinalOutput_U95$KFold==b,]
        
        # Replace NA values by 0 in the 'ModelParameters' objects:
        set_of_models[is.na(set_of_models)] <- 0
        set_of_models_L95[is.na(set_of_models_L95)] <- 0
        set_of_models_U95[is.na(set_of_models_U95)] <- 0
        
        # Create an empty dataframe to store the estimated discovery probability (95% CI) in the year 2015 (all species, all AFT models):
        DiscProb_dt<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data_subset)))
        DiscProb_dt_L95<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data_subset)))
        DiscProb_dt_U95<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data_subset)))
        DiscProb_dt[is.na(DiscProb_dt)] = 0
        DiscProb_dt_L95[is.na(DiscProb_dt_L95)] = 0
        DiscProb_dt_U95[is.na(DiscProb_dt_U95)] = 0
        
        # Same as above, but to store the predicted year of discovery (= the year in which discovery probability is 0.50):
        YearPred_dt<-as.data.table(matrix(nrow=nrow(set_of_models), ncol=nrow(trait_data_subset))) 
        YearPred_dt[is.na(YearPred_dt)] = 0
        gc()
        
        # Get the discovery probability (95% CI) in 2015 and estimated year of discovery for all species and AFT models:
        cl<-makePSOCKcluster(detectCores()-1)
        registerDoParallel(cl)
        getDoParWorkers()
        for(i in 1:nrow(set_of_models)){
          
          if(unique(set_of_models$FamilyError)=="exp"){
            TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                        .combine = 'cbind',
                                        .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                          qexp(pct, rate=exp(as.numeric(set_of_models[i,23] +
                                                                          trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                          trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                          trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                          trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                          trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                          trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                          trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                          trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                          trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                          trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                          trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                        }
            ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                      .combine = 'cbind',
                                      .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                        pexp(pct, rate=exp(as.numeric(set_of_models[i,23] +
                                                                        trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                        trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                        trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                        trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                        trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                        trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                        trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                        trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                        trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                        trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                        trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                      }
            ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pexp(pct, rate=exp(as.numeric(set_of_models_L95[i,12] +
                                                                            trait_data_subset$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                            trait_data_subset$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                            trait_data_subset$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                            trait_data_subset$LogAMT[k]*set_of_models_L95[i,4] +
                                                                            trait_data_subset$LogAPP[k]*set_of_models_L95[i,5] +
                                                                            trait_data_subset$LogTS[k]*set_of_models_L95[i,6] +
                                                                            trait_data_subset$LogPS[k]*set_of_models_L95[i,7] +
                                                                            trait_data_subset$LogElevM[k]*set_of_models_L95[i,8] +
                                                                            trait_data_subset$LogPopD[k]*set_of_models_L95[i,9] +
                                                                            trait_data_subset$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                            trait_data_subset$LogRarity[k]*set_of_models_L95[i,11])))
                                          }
            ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pexp(pct, rate=exp(as.numeric(set_of_models_U95[i,12] +
                                                                            trait_data_subset$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                            trait_data_subset$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                            trait_data_subset$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                            trait_data_subset$LogAMT[k]*set_of_models_U95[i,4] +
                                                                            trait_data_subset$LogAPP[k]*set_of_models_U95[i,5] +
                                                                            trait_data_subset$LogTS[k]*set_of_models_U95[i,6] +
                                                                            trait_data_subset$LogPS[k]*set_of_models_U95[i,7] +
                                                                            trait_data_subset$LogElevM[k]*set_of_models_U95[i,8] +
                                                                            trait_data_subset$LogPopD[k]*set_of_models_U95[i,9] +
                                                                            trait_data_subset$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                            trait_data_subset$LogRarity[k]*set_of_models_U95[i,11])))
                                          }
          }
          
          if(unique(set_of_models$FamilyError)=="weibull"){
            TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                        .combine = 'cbind',
                                        .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                          qweibull(pct, shape=set_of_models[i,23], 
                                                   scale=as.numeric(set_of_models[i,25]*
                                                                      exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                            trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                            trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                            trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                            trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                            trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                            trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                            trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                            trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                            trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                            trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                        }
            ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                      .combine = 'cbind',
                                      .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                        pweibull(pct, shape=set_of_models[i,23], 
                                                 scale=as.numeric(set_of_models[i,25]*
                                                                    exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                          trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                          trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                          trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                          trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                          trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                          trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                          trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                          trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                          trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                          trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                      }
            ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pweibull(pct, shape=set_of_models_L95[i,12], 
                                                     scale=as.numeric(set_of_models_L95[i,13]*
                                                                        exp(trait_data_subset$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                              trait_data_subset$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                              trait_data_subset$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                              trait_data_subset$LogAMT[k]*set_of_models_L95[i,4] +
                                                                              trait_data_subset$LogAPP[k]*set_of_models_L95[i,5] +
                                                                              trait_data_subset$LogTS[k]*set_of_models_L95[i,6] +
                                                                              trait_data_subset$LogPS[k]*set_of_models_L95[i,7] +
                                                                              trait_data_subset$LogElevM[k]*set_of_models_L95[i,8] +
                                                                              trait_data_subset$LogPopD[k]*set_of_models_L95[i,9] +
                                                                              trait_data_subset$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                              trait_data_subset$LogRarity[k]*set_of_models_L95[i,11])))
                                          }
            ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pweibull(pct, shape=set_of_models_U95[i,12], 
                                                     scale=as.numeric(set_of_models_U95[i,13]*
                                                                        exp(trait_data_subset$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                              trait_data_subset$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                              trait_data_subset$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                              trait_data_subset$LogAMT[k]*set_of_models_U95[i,4] +
                                                                              trait_data_subset$LogAPP[k]*set_of_models_U95[i,5] +
                                                                              trait_data_subset$LogTS[k]*set_of_models_U95[i,6] +
                                                                              trait_data_subset$LogPS[k]*set_of_models_U95[i,7] +
                                                                              trait_data_subset$LogElevM[k]*set_of_models_U95[i,8] +
                                                                              trait_data_subset$LogPopD[k]*set_of_models_U95[i,9] +
                                                                              trait_data_subset$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                              trait_data_subset$LogRarity[k]*set_of_models_U95[i,11])))
                                          }
            
            # colnames = species & rownames = discovery probability values
            names(TimeToDisc_per_spp)<-trait_data_subset$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
            
            # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
            names(ProbDisc_per_spp)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp)<-pct 
            names(ProbDisc_per_spp_L95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
            names(ProbDisc_per_spp_U95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
          }
          
          if(unique(set_of_models$FamilyError)=="lnorm"){
            TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                        .combine = 'cbind',
                                        .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                          qlnorm(pct, meanlog=as.numeric(set_of_models[i,23]*
                                                                           exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                                 trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                                 trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                                 trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                                 trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                                 trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                                 trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                                 trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                                 trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                                 trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                                 trait_data_subset$LogRarity[k]*set_of_models[i,11],
                                                                               sdlog=set_of_models[i,25])))
                                        }
            ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                      .combine = 'cbind',
                                      .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                        plnorm(pct, meanlog=as.numeric(set_of_models[i,23]*
                                                                         exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                               trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                               trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                               trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                               trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                               trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                               trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                               trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                               trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                               trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                               trait_data_subset$LogRarity[k]*set_of_models[i,11],
                                                                             sdlog=set_of_models[i,25])))
                                      }
            ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            plnorm(pct, meanlog=as.numeric(set_of_models_L95[i,12]*
                                                                             exp(trait_data_subset$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                                   trait_data_subset$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                                   trait_data_subset$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                                   trait_data_subset$LogAMT[k]*set_of_models_L95[i,4] +
                                                                                   trait_data_subset$LogAPP[k]*set_of_models_L95[i,5] +
                                                                                   trait_data_subset$LogTS[k]*set_of_models_L95[i,6] +
                                                                                   trait_data_subset$LogPS[k]*set_of_models_L95[i,7] +
                                                                                   trait_data_subset$LogElevM[k]*set_of_models_L95[i,8] +
                                                                                   trait_data_subset$LogPopD[k]*set_of_models_L95[i,9] +
                                                                                   trait_data_subset$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                                   trait_data_subset$LogRarity[k]*set_of_models_L95[i,11],
                                                                                 sdlog=set_of_models_L95[i,13])))
                                          }
            ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            plnorm(pct, meanlog=as.numeric(set_of_models_U95[i,12]*
                                                                             exp(trait_data_subset$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                                   trait_data_subset$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                                   trait_data_subset$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                                   trait_data_subset$LogAMT[k]*set_of_models_U95[i,4] +
                                                                                   trait_data_subset$LogAPP[k]*set_of_models_U95[i,5] +
                                                                                   trait_data_subset$LogTS[k]*set_of_models_U95[i,6] +
                                                                                   trait_data_subset$LogPS[k]*set_of_models_U95[i,7] +
                                                                                   trait_data_subset$LogElevM[k]*set_of_models_U95[i,8] +
                                                                                   trait_data_subset$LogPopD[k]*set_of_models_U95[i,9] +
                                                                                   trait_data_subset$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                                   trait_data_subset$LogRarity[k]*set_of_models_U95[i,11],
                                                                                 sdlog=set_of_models_U95[i,13])))
                                          }
            
            # colnames = species & rownames = discovery probability values
            names(TimeToDisc_per_spp)<-trait_data_subset$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
            
            # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
            names(ProbDisc_per_spp)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp)<-pct 
            names(ProbDisc_per_spp_L95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
            names(ProbDisc_per_spp_U95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
          }
          
          if(unique(set_of_models$FamilyError)=="llogis"){
            TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                        .combine = 'cbind',
                                        .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                          qllogis(pct, shape=set_of_models[i,23], 
                                                  scale=as.numeric(set_of_models[i,25]*
                                                                     exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                           trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                           trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                           trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                           trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                           trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                           trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                           trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                           trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                           trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                           trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                        }
            ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                      .combine = 'cbind',
                                      .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                        pllogis(pct, shape=set_of_models[i,23], 
                                                scale=as.numeric(set_of_models[i,25]*
                                                                   exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                         trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                         trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                         trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                         trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                         trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                         trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                         trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                         trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                         trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                         trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                      }
            ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pllogis(pct, shape=set_of_models_L95[i,12], 
                                                    scale=as.numeric(set_of_models_L95[i,13]*
                                                                       exp(trait_data_subset$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                             trait_data_subset$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                             trait_data_subset$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                             trait_data_subset$LogAMT[k]*set_of_models_L95[i,4] +
                                                                             trait_data_subset$LogAPP[k]*set_of_models_L95[i,5] +
                                                                             trait_data_subset$LogTS[k]*set_of_models_L95[i,6] +
                                                                             trait_data_subset$LogPS[k]*set_of_models_L95[i,7] +
                                                                             trait_data_subset$LogElevM[k]*set_of_models_L95[i,8] +
                                                                             trait_data_subset$LogPopD[k]*set_of_models_L95[i,9] +
                                                                             trait_data_subset$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                             trait_data_subset$LogRarity[k]*set_of_models_L95[i,11])))
                                          }
            ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pllogis(pct, shape=set_of_models_U95[i,12], 
                                                    scale=as.numeric(set_of_models_U95[i,13]*
                                                                       exp(trait_data_subset$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                             trait_data_subset$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                             trait_data_subset$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                             trait_data_subset$LogAMT[k]*set_of_models_U95[i,4] +
                                                                             trait_data_subset$LogAPP[k]*set_of_models_U95[i,5] +
                                                                             trait_data_subset$LogTS[k]*set_of_models_U95[i,6] +
                                                                             trait_data_subset$LogPS[k]*set_of_models_U95[i,7] +
                                                                             trait_data_subset$LogElevM[k]*set_of_models_U95[i,8] +
                                                                             trait_data_subset$LogPopD[k]*set_of_models_U95[i,9] +
                                                                             trait_data_subset$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                             trait_data_subset$LogRarity[k]*set_of_models_U95[i,11])))
                                          }
            
            # colnames = species & rownames = discovery probability values
            names(TimeToDisc_per_spp)<-trait_data_subset$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
            
            # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
            names(ProbDisc_per_spp)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp)<-pct 
            names(ProbDisc_per_spp_L95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
            names(ProbDisc_per_spp_U95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
          }
          
          if(unique(set_of_models$FamilyError)=="gamma"){
            TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                        .combine = 'cbind',
                                        .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                          qgamma(pct, shape=set_of_models[i,23], 
                                                 rate=as.numeric(set_of_models[i,25]*
                                                                   exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                         trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                         trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                         trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                         trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                         trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                         trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                         trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                         trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                         trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                         trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                        }
            ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                      .combine = 'cbind',
                                      .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                        pgamma(pct, shape=set_of_models[i,23], 
                                               rate=as.numeric(set_of_models[i,25]*
                                                                 exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                       trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                       trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                       trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                       trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                       trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                       trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                       trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                       trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                       trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                       trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                      }
            ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pgamma(pct, shape=set_of_models_L95[i,12], 
                                                   rate=as.numeric(set_of_models_L95[i,13]*
                                                                     exp(trait_data_subset$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                           trait_data_subset$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                           trait_data_subset$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                           trait_data_subset$LogAMT[k]*set_of_models_L95[i,4] +
                                                                           trait_data_subset$LogAPP[k]*set_of_models_L95[i,5] +
                                                                           trait_data_subset$LogTS[k]*set_of_models_L95[i,6] +
                                                                           trait_data_subset$LogPS[k]*set_of_models_L95[i,7] +
                                                                           trait_data_subset$LogElevM[k]*set_of_models_L95[i,8] +
                                                                           trait_data_subset$LogPopD[k]*set_of_models_L95[i,9] +
                                                                           trait_data_subset$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                           trait_data_subset$LogRarity[k]*set_of_models_L95[i,11])))
                                          }
            ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pgamma(pct, shape=set_of_models_U95[i,12], 
                                                   rate=as.numeric(set_of_models_U95[i,13]*
                                                                     exp(trait_data_subset$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                           trait_data_subset$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                           trait_data_subset$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                           trait_data_subset$LogAMT[k]*set_of_models_U95[i,4] +
                                                                           trait_data_subset$LogAPP[k]*set_of_models_U95[i,5] +
                                                                           trait_data_subset$LogTS[k]*set_of_models_U95[i,6] +
                                                                           trait_data_subset$LogPS[k]*set_of_models_U95[i,7] +
                                                                           trait_data_subset$LogElevM[k]*set_of_models_U95[i,8] +
                                                                           trait_data_subset$LogPopD[k]*set_of_models_U95[i,9] +
                                                                           trait_data_subset$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                           trait_data_subset$LogRarity[k]*set_of_models_U95[i,11])))
                                          }
            
            # colnames = species & rownames = discovery probability values
            names(TimeToDisc_per_spp)<-trait_data_subset$Binomial; rownames(TimeToDisc_per_spp)<-(1-pct) # 
            
            # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
            names(ProbDisc_per_spp)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp)<-pct 
            names(ProbDisc_per_spp_L95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
            names(ProbDisc_per_spp_U95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
          }
          
          if(unique(set_of_models$FamilyError)=="gompertz"){
            TimeToDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                        .combine = 'cbind',
                                        .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                          qgompertz(pct, shape=set_of_models[i,23],
                                                    rate=as.numeric(set_of_models[i,25]*
                                                                      exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                            trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                            trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                            trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                            trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                            trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                            trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                            trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                            trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                            trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                            trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                        }
            ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data_subset),
                                      .combine = 'cbind',
                                      .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                        pgompertz(pct, shape=set_of_models[i,23],
                                                  rate=as.numeric(set_of_models[i,25]*
                                                                    exp(trait_data_subset$LogBodySize[k]*set_of_models[i,1] +
                                                                          trait_data_subset$LogTaxPerClade[k]*set_of_models[i,2] +
                                                                          trait_data_subset$LogRangeSize[k]*set_of_models[i,3] +
                                                                          trait_data_subset$LogAMT[k]*set_of_models[i,4] +
                                                                          trait_data_subset$LogAPP[k]*set_of_models[i,5] +
                                                                          trait_data_subset$LogTS[k]*set_of_models[i,6] +
                                                                          trait_data_subset$LogPS[k]*set_of_models[i,7] +
                                                                          trait_data_subset$LogElevM[k]*set_of_models[i,8] +
                                                                          trait_data_subset$LogPopD[k]*set_of_models[i,9] +
                                                                          trait_data_subset$LogTaxPerBiome[k]*set_of_models[i,10] +
                                                                          trait_data_subset$LogRarity[k]*set_of_models[i,11])))
                                      }
            ProbDisc_per_spp_L95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pgompertz(pct, shape=set_of_models_L95[i,12],
                                                      rate=as.numeric(set_of_models_L95[i,13]*
                                                                        exp(trait_data_subset$LogBodySize[k]*set_of_models_L95[i,1] +
                                                                              trait_data_subset$LogTaxPerClade[k]*set_of_models_L95[i,2] +
                                                                              trait_data_subset$LogRangeSize[k]*set_of_models_L95[i,3] +
                                                                              trait_data_subset$LogAMT[k]*set_of_models_L95[i,4] +
                                                                              trait_data_subset$LogAPP[k]*set_of_models_L95[i,5] +
                                                                              trait_data_subset$LogTS[k]*set_of_models_L95[i,6] +
                                                                              trait_data_subset$LogPS[k]*set_of_models_L95[i,7] +
                                                                              trait_data_subset$LogElevM[k]*set_of_models_L95[i,8] +
                                                                              trait_data_subset$LogPopD[k]*set_of_models_L95[i,9] +
                                                                              trait_data_subset$LogTaxPerBiome[k]*set_of_models_L95[i,10] +
                                                                              trait_data_subset$LogRarity[k]*set_of_models_L95[i,11])))
                                          }
            ProbDisc_per_spp_U95<-foreach(k = 1:nrow(trait_data_subset),
                                          .combine = 'cbind',
                                          .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                            pgompertz(pct, shape=set_of_models_U95[i,12],
                                                      rate=as.numeric(set_of_models_U95[i,13]*
                                                                        exp(trait_data_subset$LogBodySize[k]*set_of_models_U95[i,1] +
                                                                              trait_data_subset$LogTaxPerClade[k]*set_of_models_U95[i,2] +
                                                                              trait_data_subset$LogRangeSize[k]*set_of_models_U95[i,3] +
                                                                              trait_data_subset$LogAMT[k]*set_of_models_U95[i,4] +
                                                                              trait_data_subset$LogAPP[k]*set_of_models_U95[i,5] +
                                                                              trait_data_subset$LogTS[k]*set_of_models_U95[i,6] +
                                                                              trait_data_subset$LogPS[k]*set_of_models_U95[i,7] +
                                                                              trait_data_subset$LogElevM[k]*set_of_models_U95[i,8] +
                                                                              trait_data_subset$LogPopD[k]*set_of_models_U95[i,9] +
                                                                              trait_data_subset$LogTaxPerBiome[k]*set_of_models_U95[i,10] +
                                                                              trait_data_subset$LogRarity[k]*set_of_models_U95[i,11])))
                                          }
          }
          
          # colnames = species & rownames = discovery probability values
          names(TimeToDisc_per_spp)<-trait_data_subset$Binomial; rownames(TimeToDisc_per_spp)<-pct # 
          
          # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018)
          names(ProbDisc_per_spp)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp)<-pct 
          names(ProbDisc_per_spp_L95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_L95)<-pct
          names(ProbDisc_per_spp_U95)<-trait_data_subset$Binomial; rownames(ProbDisc_per_spp_U95)<-pct
          
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
          for (k in 1:nrow(trait_data_subset)){data.table::set(DiscProb_dt, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015[k])}
          for (k in 1:nrow(trait_data_subset)){data.table::set(YearPred_dt, i=as.integer(i), j=as.integer(k), value=YearPred[k])}
          for (k in 1:nrow(trait_data_subset)){data.table::set(DiscProb_dt_L95, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015_L95[k])}
          for (k in 1:nrow(trait_data_subset)){data.table::set(DiscProb_dt_U95, i=as.integer(i), j=as.integer(k), value=ProbDisc_2015_U95[k])}
          
        } # end of i for loop (across all AFT models)
        stopCluster(cl)
        
        # Register the time period, kfold,
        DiscProb_dt$TimePeriod<-TimePeriod[t]
        DiscProb_dt_L95$TimePeriod<-TimePeriod[t]
        DiscProb_dt_U95$TimePeriod<-TimePeriod[t]
        YearPred_dt$TimePeriod<-TimePeriod[t]
        
        DiscProb_dt$CV_PartitionSize<-PercSppModelled[c]
        DiscProb_dt_L95$CV_PartitionSize<-PercSppModelled[c]
        DiscProb_dt_U95$CV_PartitionSize<-PercSppModelled[c]
        YearPred_dt$CV_PartitionSize<-PercSppModelled[c]
        
        DiscProb_dt$KFold<-b
        DiscProb_dt_L95$KFold<-b
        DiscProb_dt_U95$KFold<-b
        YearPred_dt$KFold<-b
        
        DiscProb_dt$Class<-VertGroups[v]
        DiscProb_dt_L95$Class<-VertGroups[v]
        DiscProb_dt_U95$Class<-VertGroups[v]
        YearPred_dt$Class<-VertGroups[v]
        
        # Store the results in the empty lists prepared a priori:
        list1_1[[b]]<-DiscProb_dt
        list2_1[[b]]<-DiscProb_dt_L95
        list3_1[[b]]<-DiscProb_dt_U95
        list4_1[[b]]<-YearPred_dt
      
      } # end of b for loop (k-fold partitions)
      
      list1_2[[t]]<-list1_1 # DiscProb_dt
      list2_2[[t]]<-list2_1 # DiscProb_dt_L95
      list3_2[[t]]<-list3_1 # DiscProb_dt_U95
      list4_2[[t]]<-list4_1 # YearPred_dt
      
    } # end of t for loop (time periods)
    
    # Export outputs:
    save(list1_2, file=paste0("ModelValidation/DiscProb_dt_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(list2_2, file=paste0("ModelValidation/DiscProb_dt_L95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(list3_2, file=paste0("ModelValidation/DiscProb_dt_U95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(list4_2, file=paste0("ModelValidation/YearPred_dt_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    
  } # end of c for loop (cross-validation partition sizes)
  
} # end of v for loop (vertebrate groups)

#####

# STEP 4 - GET SPECIES-LEVEL AVG WEIGHTED DISCOVERY METRICS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
##########################################################################################################################
# STEP 4 - GET SPECIES-LEVEL AVG WEIGHTED DISCOVERY METRICS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
rm(list=ls())

# Function to load RData using user-specified object names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create objects to guide for loops:
TimePeriod<-c(1758, seq(from=1760, to=1970, by=10)) # the starting date of each time period
VertGroups<-c("Amphibia", "Reptilia", "Mammalia", "Aves")
PercSppModelled<-c("25T75V", "50T50V", "75T25V", "90T10V") # T = training-fold / V = validation-fold
Kfold_partitions<-c(4, 2, 4, 10) # number of k-fold partitions represented in each levels of 'PercSppModelled'

# Compute the average weighted discovery metrics for each time period and vertebrate group:
for(v in 1:4){ # v-number of vertebrate groups
  
  # Load the response and predictor variables:
  raw_trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  raw_trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>=1759 & raw_trait_data$YearOfDescription<=2014,]
  
  # Get the response variable (time-to-event) standardized between 0 and 1:
  raw_trait_data$Time<-(raw_trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
  raw_trait_data$Censor<-1 # Censor variable (it informs if the event happened)
  
  # Filter the dataset to represent only one vertebrate group:
  raw_trait_data<-raw_trait_data[raw_trait_data$Class==VertGroups[v],]
  raw_trait_data<-droplevels(raw_trait_data) # drop unused levels
  raw_trait_data<-raw_trait_data[complete.cases(raw_trait_data), ] # remove species without data available
  raw_trait_data<-raw_trait_data[, .(Binomial, Class, YearOfDescription)]
  
  # Create empty lists to store some outputs computed across different time periods (see below):
  VertList<-list() # this will receive the object DiscProb_dt
   
  for(c in 1:2){ # c-cross validation partition sizes
    
    # Perform computations for each cross-validation partition:
    for(t in 1:length(TimePeriod)){
      
      # Filter species for a given time period:
      trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>TimePeriod[t],]
      
      # Randomly split species data in training- and validation-data subsets:
      set.seed(12)
      if(c==1){ # 25% of species in the training-fold and 75% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_fold<-list()
        species_fold[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/4)],1]
        species_fold[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/4)+1):(2*ceiling(length(random_rows)/4))],1]
        species_fold[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/4)+1):(3*ceiling(length(random_rows)/4))],1]
        species_fold[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/4)+1):length(random_rows)],1]
        
      } # 25% of species in the training-fold and 75% in the validation-fold
      if(c==2){ # 50% of species in the training-fold and 50% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_fold<-list()
        species_fold[[1]]<-trait_data[random_rows[1:floor(length(random_rows)/2)],1]
        species_fold[[2]]<-trait_data[random_rows[(floor(length(random_rows)/2)+1):length(random_rows)],1]
        
      } # 50% of species in the training-fold and 50% in the validation-fold
      if(c==3){ # 75% of species in the training-fold and 25% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_to_discard<-list()
        species_fold<-list()
        species_to_discard[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/4)],1]
        species_to_discard[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/4)+1):(2*ceiling(length(random_rows)/4))],1]
        species_to_discard[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/4)+1):(3*ceiling(length(random_rows)/4))],1]
        species_to_discard[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/4)+1):length(random_rows)],1]
        species_fold[[1]]<-trait_data[Binomial %ni% as.character(species_to_discard[[1]]$Binomial)]
        species_fold[[2]]<-trait_data[Binomial %ni% as.character(species_to_discard[[2]]$Binomial)]
        species_fold[[3]]<-trait_data[Binomial %ni% as.character(species_to_discard[[3]]$Binomial)]
        species_fold[[4]]<-trait_data[Binomial %ni% as.character(species_to_discard[[4]]$Binomial)]
        rm(species_to_discard)
        
      } # 75% of species in the training-fold and 25% in the validation-fold
      if(c==4){ # 90% of species in the training-fold and 10% in the validation-fold
        
        random_rows<-sample(x=c(1:nrow(trait_data)), size=nrow(trait_data), replace=F)
        species_to_discard<-list()
        species_fold<-list()
        species_to_discard[[1]]<-trait_data[random_rows[1:ceiling(length(random_rows)/10)],1]
        species_to_discard[[2]]<-trait_data[random_rows[(ceiling(length(random_rows)/10)+1):(2*ceiling(length(random_rows)/10))],1]
        species_to_discard[[3]]<-trait_data[random_rows[(2*ceiling(length(random_rows)/10)+1):(3*ceiling(length(random_rows)/10))],1]
        species_to_discard[[4]]<-trait_data[random_rows[(3*ceiling(length(random_rows)/10)+1):(4*ceiling(length(random_rows)/10))],1]
        species_to_discard[[5]]<-trait_data[random_rows[(4*ceiling(length(random_rows)/10)+1):(5*ceiling(length(random_rows)/10))],1]
        species_to_discard[[6]]<-trait_data[random_rows[(5*ceiling(length(random_rows)/10)+1):(6*ceiling(length(random_rows)/10))],1]
        species_to_discard[[7]]<-trait_data[random_rows[(6*ceiling(length(random_rows)/10)+1):(7*ceiling(length(random_rows)/10))],1]
        species_to_discard[[8]]<-trait_data[random_rows[(7*ceiling(length(random_rows)/10)+1):(8*ceiling(length(random_rows)/10))],1]
        species_to_discard[[9]]<-trait_data[random_rows[(8*ceiling(length(random_rows)/10)+1):(9*ceiling(length(random_rows)/10))],1]
        species_to_discard[[10]]<-trait_data[random_rows[(9*ceiling(length(random_rows)/10)+1):length(random_rows)],1]
        species_fold[[1]]<-trait_data[Binomial %ni% as.character(species_to_discard[[1]]$Binomial)]
        species_fold[[2]]<-trait_data[Binomial %ni% as.character(species_to_discard[[2]]$Binomial)]
        species_fold[[3]]<-trait_data[Binomial %ni% as.character(species_to_discard[[3]]$Binomial)]
        species_fold[[4]]<-trait_data[Binomial %ni% as.character(species_to_discard[[4]]$Binomial)]
        species_fold[[5]]<-trait_data[Binomial %ni% as.character(species_to_discard[[5]]$Binomial)]
        species_fold[[6]]<-trait_data[Binomial %ni% as.character(species_to_discard[[6]]$Binomial)]
        species_fold[[7]]<-trait_data[Binomial %ni% as.character(species_to_discard[[7]]$Binomial)]
        species_fold[[8]]<-trait_data[Binomial %ni% as.character(species_to_discard[[8]]$Binomial)]
        species_fold[[9]]<-trait_data[Binomial %ni% as.character(species_to_discard[[9]]$Binomial)]
        species_fold[[10]]<-trait_data[Binomial %ni% as.character(species_to_discard[[10]]$Binomial)]
        rm(species_to_discard)
        
      } # 90% of species in the training-fold and 10% in the validation-fold
      for(d in 1:length(species_fold)){species_fold[[d]]<-droplevels(species_fold[[d]])} # discard non-used levels
      
      # Create empty lists to store some outputs computed across different time periods (see below):
      VertList_sublist<-list() # this will receive the object DiscProb_dt
      
      # Perform computations for each k-fold partitions:
      for(b in 1:Kfold_partitions[c]){
        
        # Load the model parameters and filter for the selected TimePeriod and KFold:
        FinalOutput<-as.data.frame(loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_", VertGroups[v], "_", PercSppModelled[c], ".RData"))))
        FinalOutput_L95<-as.data.frame(loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_L95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))))
        FinalOutput_U95<-as.data.frame(loadRData(file.path(paste0("ModelValidation/ModelParametersAcrossTimeAcrossKfolds_U95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))))
        set_of_models<-FinalOutput[FinalOutput$TimePeriod==TimePeriod[t] & FinalOutput$KFold==b,]
        set_of_models_L95<-FinalOutput_L95[FinalOutput_L95$TimePeriod==TimePeriod[t] & FinalOutput_L95$KFold==b,]
        set_of_models_U95<-FinalOutput_U95[FinalOutput_U95$TimePeriod==TimePeriod[t] & FinalOutput_U95$KFold==b,]
        
        # Load the RData object for the vertebrate group v and time period t:
        DiscProb_dt<-loadRData(paste0("ModelValidation/DiscProb_dt_", VertGroups[v], "_", PercSppModelled[c], ".RData"))[[t]][[b]]
        DiscProb_dt_L95<-loadRData(paste0("ModelValidation/DiscProb_dt_L95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))[[t]][[b]]
        DiscProb_dt_U95<-loadRData(paste0("ModelValidation/DiscProb_dt_U95_", VertGroups[v], "_", PercSppModelled[c], ".RData"))[[t]][[b]]
        YearPred_dt<-loadRData(paste0("ModelValidation/YearPred_dt_", VertGroups[v], "_", PercSppModelled[c], ".RData"))[[t]][[b]]
        
        # Filter the species described within the time period of investigation:
        trait_data_subset<-trait_data[Binomial %in% as.character(species_fold[[b]]$Binomial)]
        
        # Keep only the columns representing species in the dataset:
        cols<-names(DiscProb_dt)[1:(ncol(DiscProb_dt)-4)]
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
        trait_data_subset$DiscProb<-AvgW_DiscProb_dt
        trait_data_subset$DiscProb_L95<-AvgW_DiscProb_dt_L95
        trait_data_subset$DiscProb_U95<-AvgW_DiscProb_dt_U95
        trait_data_subset$YearPred<-AvgW_YearPred_dt
        trait_data_subset$TimePeriod<-TimePeriod[t]
        trait_data_subset$CV_PartitionSize<-PercSppModelled[c]
        trait_data_subset$KFold<-b
        
        # Store the outputs for each vertebrate group:
        VertList_sublist[[b]]<-trait_data_subset
        rm(AvgW_DiscProb_dt, AvgW_DiscProb_dt_L95, AvgW_DiscProb_dt_U95, AvgW_YearPred_dt, trait_data_subset,
           DiscProb_dt, DiscProb_dt_L95, DiscProb_dt_U95, YearPred_dt, set_of_models, cols, spp_columns)
      
        } # end of b for loop (across kfolds)
      
      # Store the outputs for each vertebrate group:
      VertList[[t]]<-VertList_sublist

    } # end of t for loop across time periods
    
    # Export outputs:
    save(VertList, file=paste0("ModelValidation/SpeciesLevelEstimates_ModelValidation_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    
  } # end of c for loop (cross-validation partition sizes)
  
} # end of v for loop (vertebrate groups)

#####

# STEP 5 - TAXON-LEVEL ANALYSIS, GET DISCOVERY METRICS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
##########################################################################################################################
# STEP 5 - TAXON-LEVEL ANALYSIS, GET DISCOVERY METRICS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
rm(list=ls())

# Function to load RData using user-specified object names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create a function to compute the geometric mean for a set of species:
geometric.mean <- function(x,na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm)) }

'%ni%' <- Negate('%in%')

# Create objects to guide for loops:
TimePeriod<-c(1758, seq(from=1760, to=1970, by=10)) # the starting date of each time period
VertGroups<-c("Amphibia", "Reptilia", "Mammalia", "Aves")
PercSppModelled<-c("25T75V", "50T50V", "75T25V", "90T10V") # T = training-fold / V = validation-fold
Kfold_partitions<-c(4, 2, 4, 10) # number of k-fold partitions represented in each levels of 'PercSppModelled'

# Manually set the taxonomic groups of interest:
non_neobatrachia<-c("Leiopelmatidae", "Ascaphidae", "Bombinatoridae", "Discoglossidae", "Alytidae", "Pipidae", 
                    "Rhinophrynidae", "Scaphiopodidae", "Pelodytidae", "Megophryidae","Pelobatidae")
hyloidea_families<-c("Calyptocephalellidae", "Myobatrachidae", "Ceuthomantidae", "Brachycephalidae", "Eleutherodactylidae", "Craugastoridae",
                     "Hemiphractidae", "Hylidae", "Pelodryadidae", "Phyllomedusidae", "Bufonidae", "Dendrobatidae", "Aromobatidae", 
                     "Allophrynidae", "Centrolenidae", "Leptodactylidae", "Ceratophryidae", "Odontophrynidae", "Cycloramphidae", "Alsodidae",
                     "Hylodidae", "Telmatobiidae", "Batrachylidae", "Rhinodermatidae","Strabomantidae")
ranoidea_families<-c("Arthroleptidae", "Brevicipitidae", "Ceratobatrachidae", "Conrauidae", "Dicroglossidae", "Hemisotidae", "Hyperoliidae", 
                     "Mantellidae", "Micrixalidae", "Microhylidae",  "Nyctibatrachidae", "Odontobatrachidae", "Petropedetidae", 
                     "Phrynobatrachidae", "Ptychadenidae", "Pyxicephalidae", "Ranidae", "Ranixalidae", "Rhacophoridae")
neobatrachia_others<-c("Nasikabatrachidae", "Sooglossidae", "Calyptocephalellidae", "Myobatrachidae", "Heleophrynidae", "Limnodynastidae")
cryptodira_families<-c("Chelydridae", "Emydidae", "Geoemydidae", "Platysternidae", "Carettochelyidae", "Dermatemydidae",
                       "Kinosternidae", "Testudinidae", "Trionychidae")
pleurodira_families<-c("Chelidae", "Pelomedusidae", "Podocnemididae")
gekkota_families<-c("Gekkonidae", "Carphodactylidae", "Diplodactylidae", "Phyllodactylidae", "Sphaerodactylidae", "Pygopodidae", "Eublepharidae")
iguania_families<-c("Chamaeleonidae", "Iguanidae", "Agamidae", "Dactyloidae", "Opluridae", "Leiocephalidae", "Liolaemidae", "Phrynosomatidae", "Tropiduridae",
                    "Leiosauridae", "Polychrotidae", "Corytophanidae", "Crotaphytidae", "Hoplocercidae")
scincoidea_families<-c("Gerrhosauridae", "Cordylidae", "Scincidae", "SCINCIDAE", "Xantusiidae")
lacertiformes_families<-c("Teiidae",  "Gymnophthalmidae", "Lacertidae", "Amphisbaenidae", "Cadeidae", "Trogonophiidae", "Bipedidae", "Blanidae", "Rhineuridae")
anguimorpha_families<-c( "Anguidae", "Diploglossidae", "Anniellidae", "Xenosauridae", "Helodermatidae", "Shinisauridae", "Lanthanotidae", "Varanidae")
dibamia_family<-"Dibamidae"
serpentes_families<-c("Typhlopidae", "Colubridae", "Leptotyphlopidae", "Natricidae", "Gerrhopilidae", "Anomalepididae", "Dipsadidae", "Homalopsidae",
                      "Lamprophiidae", "Xenodermatidae", "Pseudoxenodontidae", "Viperidae", "Uropeltidae", "Tropidophiidae", "Elapidae", "Anomochilidae",
                      "Xenophidiidae", "Xenotyphlopidae", "Cylindrophiidae", "Pareatidae", "Boidae", "Xenopeltidae", "Pythonidae", "Acrochordidae",
                      "Aniliidae", "Bolyeridae", "Loxocemidae")
oscines_passeri<-c("Menuridae", "Atrichornithidae", "Ptilonorhynchidae", "Climacteridae", "Maluridae", "Dasyornithidae", "Meliphagidae", 
                   "Pardalotidae", "Acanthizidae", "Orthonychidae", "Pomatostomidae", "Melanocharitidae", "Cnemophilidae", "Callaeidae",
                   "Notiomystidae", "Mohouidae", "Neosittidae", "Campephagidae", "Cinclosomatidae", "Pachycephalidae", "Falcunculidae",
                   "Oreoicidae", "Psophodidae", "Eulacestomatidae", "Paramythiidae", "Vireonidae", "Oriolidae", "Rhagologidae",
                   "Machaerirhynchidae", "Artamidae", "Platysteiridae", "Vangidae", "Pityriasidae", "Aegithinidae", "Malaconotidae",
                   "Dicruridae", "Rhipiduridae", "Laniidae", "Corvidae", "Monarchidae", "Corcoracidae", "Melampittidae", "Ifritidae",
                   "Paradisaeidae", "Eupetidae", "Petroicidae", "Promeropidae", "Dicaeidae", "Nectariniidae", "Irenidae", "Urocynchramidae",
                   "Prunellidae", "Peucedramidae", "Ploceidae", "Estrildidae", "Viduidae", "Passeridae", "Motacillidae", "Fringillidae",
                   "Calcariidae", "Rhodinocichlidae", "Emberizidae", "Passerellidae", "Phaenicophilidae", "Zeledoniidae", "Parulidae",
                   "Icteridae", "Calyptophilidae", "Mitrospingidae", "Cardinalidae", "Thraupidae", "Hyliotidae", "Stenostiridae", "Paridae",
                   "Remizidae", "Nicatoridae", "Alaudidae", "Panuridae", "Macrosphenidae", "Cisticolidae", "Locustellidae", "Donacobiidae",
                   "Bernieridae", "Acrocephalidae", "Pnoepygidae", "Hirundinidae", "Pycnonotidae", "Phylloscopidae", "Scotocercidae",
                   "Aegithalidae", "Sylviidae", "Zosteropidae", "Timaliidae", "Pellorneidae", "Leiothrichidae", "Regulidae", "Dulidae", 
                   "Bombycillidae", "Hypocoliidae", "Ptiliogonatidae", "Mohoidae", "Certhiidae", "Sittidae", "Troglodytidae", "Polioptilidae",
                   "Buphagidae", "Mimidae", "Sturnidae", "Cinclidae", "Muscicapidae", "Turdidae",
                   "Chloropseidae", "Colluricinclidae", "Cracticidae", "Picathartidae", "Reguliidae", "Rhabdornithidae", "Callaeatidae", "Pityriaseidae") # songbirds
suboscines_tyranni<-c("Pittidae", "Eurylaimidae", "Philepittidae", "Calyptomenidae", "Sapayoidae", "Pipridae", "Cotingidae", 
                      "Oxyruncidae", "Onychorhynchidae", "Tityridae", "Platyrinchidae", "Pipritidae", "Tachurisidae", "Pipromorphidae", 
                      "Tyrannidae", "Thamnophilidae", "Melanopareiidae", "Conopophagidae", "Grallariidae", "Rhinocryptidae", "Formicariidae",
                      "Scleruridae", "Dendrocolaptidae", "Furnariidae", "Sapayoaidae")

# Create a list to store the results for all vertebrates:
VertList<-list()

# Compute the average weighted discovery metrics for each time period and vertebrate group:
for(v in 1:4){ # v-number of vertebrate groups
  
  # Load the response and predictor variables:
  raw_trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  raw_trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>=1759 & raw_trait_data$YearOfDescription<=2014,]
  
  # Filter the dataset to represent only one vertebrate group:
  raw_trait_data<-raw_trait_data[raw_trait_data$Class==VertGroups[v],]
  raw_trait_data<-droplevels(raw_trait_data) # drop unused levels
  raw_trait_data<-raw_trait_data[complete.cases(raw_trait_data), ] # remove species without data available
  raw_trait_data$Suborder<-raw_trait_data$Order
  raw_trait_data<-raw_trait_data[, .(Binomial, Class, Order, Suborder, Family, YearOfDescription)]
  trackthis<-which(names(raw_trait_data)=="Suborder")
  
  if(v==1){ 
    raw_trait_data[which(raw_trait_data$Family %in% as.character(ranoidea_families)),trackthis]<-"Neob. (Ranoidea)"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(hyloidea_families)),trackthis]<-"Neob. (Hyloidea)"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(neobatrachia_others)),trackthis]<-"Neobatrachia (others)"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(non_neobatrachia)),trackthis]<-"Non-Neobatrachia"
  } # replace suborder for amphibian species
  if(v==2){ 
    raw_trait_data[which(raw_trait_data$Family %in% as.character(pleurodira_families)),trackthis]<-"Pleurodira"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(cryptodira_families)),trackthis]<-"Cryptodira"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(gekkota_families)),trackthis]<-"Gekkota"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(iguania_families)),trackthis]<-"Iguania"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(scincoidea_families)),trackthis]<-"Scincoidea"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(lacertiformes_families)),trackthis]<-"Lacertiformes"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(anguimorpha_families)),trackthis]<-"Anguimorpha"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(dibamia_family)),trackthis]<-"Dibamoidea"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(serpentes_families)),trackthis]<-"Serpentes"
  } # replace suborder for reptile species
  if(v==4){ 
    raw_trait_data[which(raw_trait_data$Family %in% as.character(oscines_passeri)),trackthis]<-"Oscines"
    raw_trait_data[which(raw_trait_data$Family %in% as.character(suboscines_tyranni)),trackthis]<-"Suboscines"
    # The Passeriformes species from the Acanthisittidae family are neither Oscines nor Suboscines.
    # We used the taxonomic family to define their suborder (they are only two species):
    raw_trait_data[which(raw_trait_data$Family %in% as.character("Acanthisittidae")),1] # which species
    raw_trait_data[which(raw_trait_data$Family %in% as.character("Acanthisittidae")),trackthis]<-"Acanthisittidae"
  } # replace suborder for bird species
  
  PartitionSize_list<-list()
  
  for(c in 1:2){ # c-cross validation partition sizes
    
    TimePeriod_list<-list()
    
    # Perform computations for each cross-validation partition:
    for(t in 1:length(TimePeriod)){
      
      Kfold_list<-list()
      
      # Perform computations for each k-fold partitions:
      for(b in 1:Kfold_partitions[c]){
        
        # Filter species for a given time period:
        taxon_rank_data<-raw_trait_data[raw_trait_data$YearOfDescription>TimePeriod[t],]
        
        # Create empty lists to store the outputs on different taxonomic ranks
        list1<-list()
        list2<-list()
        
        # Perfom computations across different taxonomic ranks:
        for(l in 1:2){ #l taxonomic ranks
          
          # Set the taxonomic rank level (higher-level grouping or family):
          if(l==1){
            taxon_rank_data_subset<-taxon_rank_data[, .(Binomial, Order, Suborder)]
            names(taxon_rank_data_subset)[3]<-"TaxonRank"
          }
          if(l==2){
            taxon_rank_data_subset<-taxon_rank_data[, .(Binomial, Order, Family)]
            names(taxon_rank_data_subset)[3]<-"TaxonRank"
          }
          
          # Read the trait data for the vertebrate group v, time period t, partition size c, k-fold b:
          training_fold_spp<-loadRData(file.path(paste0("ModelValidation/SpeciesLevelEstimates_ModelValidation_",
                                                 VertGroups[v], "_", PercSppModelled[c], ".RData")))[[t]][[b]]
          training_fold_spp<-droplevels(training_fold_spp)
          training_fold_spp<-merge(training_fold_spp, taxon_rank_data_subset, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
          
          # Get the species included in the training- and validation-fold:
          validation_fold_spp<-raw_trait_data[Binomial %ni% as.character(training_fold_spp$Binomial)]
          species_left_out<-validation_fold_spp[validation_fold_spp$YearOfDescription<=TimePeriod[t],] # to be used later
          validation_fold_spp<-validation_fold_spp[validation_fold_spp$YearOfDescription>TimePeriod[t],]
          validation_fold_spp<-droplevels(validation_fold_spp); species_left_out<-droplevels(species_left_out)
          
          # Compute the clade-level metrics of discovery potential:
          spp_data_described<-training_fold_spp %>%
            dplyr::group_by(Class, Order, TaxonRank) %>% 
            dplyr::summarise(PropUnknown= (1 - geometric.mean(DiscProb, na.rm=T)),
                             PropUnknown_L95= (1 - geometric.mean(DiscProb_U95, na.rm=T)), # note that the upper and lower bounds are now inverted
                             PropUnknown_U95= (1 - geometric.mean(DiscProb_L95, na.rm=T))) # note that the upper and lower bounds are now inverted
          spp_data_described<-droplevels(spp_data_described)
          
          # Obs.: Richness need to be recomputed for all species described (without discard early descriptions) to avoid confounding
          # factors related to temporally varying per-assemblage richness:
          species_left_in<-raw_trait_data[Binomial %in% as.character(training_fold_spp$Binomial)]
          all_species_but_validation<-rbind(species_left_in, species_left_out)
          all_species_but_validation<-droplevels(all_species_but_validation)
          if(l==1){all_species_but_validation<-all_species_but_validation[,.(.N), by = .(Suborder)]}
          if(l==2){all_species_but_validation<-all_species_but_validation[,.(.N), by = .(Family)]}
          names(all_species_but_validation)<-c('TaxonRank','Richness')
          spp_data_described<-merge(spp_data_described, all_species_but_validation, by='TaxonRank', all.x=TRUE, allow.cartesian=TRUE)
          
          # Compute the estimated UnknownSR for each level of TaxonRank:
          spp_data_described$Est_UnknownSR<- (spp_data_described$Richness/(1-spp_data_described$PropUnknown)) - spp_data_described$Richness
          spp_data_described$Est_UnknownSR_L95<- (spp_data_described$Richness/(1-spp_data_described$PropUnknown_L95)) - spp_data_described$Richness
          spp_data_described$Est_UnknownSR_U95<- (spp_data_described$Richness/(1-spp_data_described$PropUnknown_U95)) - spp_data_described$Richness
          
          # Compute the observerd discoveries based on the validation-fold:
          if(l==1){spp_data_not_described_yet<-validation_fold_spp[,.(.N), by = .(Suborder)]}
          if(l==2){spp_data_not_described_yet<-validation_fold_spp[,.(.N), by = .(Family)]}
          names(spp_data_not_described_yet)<-c('TaxonRank','Obs_UnknownSR')
          
          # Bind results on estimated and observed discoveries:
          spp_data_described<-merge(spp_data_described, spp_data_not_described_yet, by='TaxonRank', all.x=TRUE, allow.cartesian=TRUE)
          
          # Register the time period, CV partition size, k-fold number, and taxonomic rank level:
          spp_data_described$TimePeriod<-TimePeriod[t]
          spp_data_described$CV_PartitionSize<-PercSppModelled[c]
          spp_data_described$KFold<-b
          if(l==1){
            spp_data_described$TaxonLevel<-"Suborder"
            list1[[t]]<-spp_data_described
          }  
          if(l==2){
            spp_data_described$TaxonLevel<-"Family"
            list2[[t]]<-spp_data_described
          }
          
          rm(spp_data_described, spp_data_not_described_yet, validation_fold_spp, training_fold_spp,
             species_left_in, species_left_out, all_species_but_validation)
          
        } # end of l for loop (taxonomic rank levels)
        
        # Store the results for each vertebrate group:
        list1<-data.table::rbindlist(list1)
        list2<-data.table::rbindlist(list2)
        Kfold_list[[b]]<-rbind(list1, list2)
        rm(list1, list2)
        
      } # end of b for loop (kfolds)
      
      TimePeriod_list[[t]]<-data.table::rbindlist(Kfold_list)
      
    } # end of t for loop (time periods)
    
    PartitionSize_list[[c]]<-data.table::rbindlist(TimePeriod_list)
  
  } # end of c for loop (cross-validation partition size)
  
  VertList[[v]]<-data.table::rbindlist(PartitionSize_list)
  
} # end of v for loop (vertebrate groups)
  
# Unlist all results for terrestrial vertebrates:
VertList<-data.table::rbindlist(VertList)

# Filter taxa with at least 5 species and export the results:
# VertList<-VertList[VertList$Richness>=5,]
fwrite(VertList, "ModelValidation/TaxonLevelEstimates_ModelValidation.csv")

#####

# STEP 6 - ASSEMBLAGE-LEVEL ANALYSIS, GET DISCOVERY METRICS ACROSS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
##########################################################################################################################
# STEP 6 - ASSEMBLAGE-LEVEL ANALYSIS, GET DISCOVERY METRICS ACROSS ACROSS TIME PERIODS, CROSS-VALIDATION PARTITION SIZES, AND K-FOLD SUBDATASETS
rm(list=ls())

# Function to load RData using user-specified object names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Create a function to compute the geometric mean for a set of species:
geometric.mean <- function(x,na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm)) }

'%ni%' <- Negate('%in%')

# Create objects to guide for loops:
TimePeriod<-c(1758, 1970) # the starting date of each time period
#TimePeriod<-c(1758, seq(from=1760, to=1970, by=10)) # the starting date of each time period
VertGroups<-c("Amphibia", "Reptilia", "Mammalia", "Aves")
PercSppModelled<-c("25T75V", "50T50V", "75T25V", "90T10V") # T = training-fold / V = validation-fold
Kfold_partitions<-c(4, 2, 4, 10) # number of k-fold partitions represented in each levels of 'PercSppModelled'
spatial_scale<-c("220km", "440km", "880km")
sub_list<-list() # list of subsampling values (maximum number of subsampled occurrences per species)
sub_list[[1]]<-c(1, 5, 10, 50, 100, 200) # 220km
sub_list[[2]]<-c(1, 5, 10, 50, 100) # 440km
sub_list[[3]]<-c(1, 5, 10) # 880km

# List of directories and files for the observed and subsampled assemblage data:
file_list<-list()
subsampled_list<-list()
file_list[[1]]<-"AssemblageLevelAnalysis/Amphibia/AssemblageData/AssemblageLevelData_Amph_"
file_list[[2]]<-"AssemblageLevelAnalysis/Reptilia/AssemblageData/AssemblageLevelData_Rept_"
file_list[[3]]<-"AssemblageLevelAnalysis/Mammalia/AssemblageData/AssemblageLevelData_Mamm_"
file_list[[4]]<-"AssemblageLevelAnalysis/Aves/AssemblageData/AssemblageLevelData_Aves_"
subsampled_list[[1]]<-"AssemblageLevelAnalysis/Amphibia/AssemblageData/Subsampling/subsampled_"
subsampled_list[[2]]<-"AssemblageLevelAnalysis/Reptilia/AssemblageData/Subsampling/subsampled_"
subsampled_list[[3]]<-"AssemblageLevelAnalysis/Mammalia/AssemblageData/Subsampling/subsampled_"
subsampled_list[[4]]<-"AssemblageLevelAnalysis/Aves/AssemblageData/Subsampling/subsampled_"
csv_scale<-c("220km.csv", "440km.csv", "880km.csv")

# Compute the metrics of discovery potential at the assemblage level:
for (v in 1:4){ # v-number of vertebrate groups
  
  for(c in 1:2){ # c-cross validation partition sizes
    
    # Get raw trait data on species names and year of description (to be used to get accumulated richness across assemblages):
    raw_trait_data<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
    raw_trait_data<-raw_trait_data[raw_trait_data$YearOfDescription>=1759 & raw_trait_data$YearOfDescription<=2014,]
    raw_trait_data<-raw_trait_data[raw_trait_data$Class==VertGroups[v],]
    raw_trait_data<-raw_trait_data[, .(Binomial, YearOfDescription)]
    
    # Create lists to store the outputs for each time period:
    TimePeriod_220km_list<-list()
    TimePeriod_440km_list<-list()
    TimePeriod_880km_list<-list()
    
    # Perform computations for each cross-validation partition:
    for(t in 1:length(TimePeriod)){
      
      # Create lists to store the outputs for each k-fold partition:
      KFold_220km_list<-list()
      KFold_440km_list<-list()
      KFold_880km_list<-list()
      
      # Perform computations for each k-fold partitions:
      for(b in 1:Kfold_partitions[c]){
        
        for (j in 1:3){ # j spatial resolutions (220, 440, 880km)
          
          # Create lists to store the outputs for each subsampling level:
          SubsamplingLevel_220km_list<-list()
          SubsamplingLevel_440km_list<-list()
          SubsamplingLevel_880km_list<-list()
          
          # Read the trait data for the vertebrate group v, time period t, partition size c, k-fold b:
          training_fold_spp<-loadRData(file.path(paste0("ModelValidation/SpeciesLevelEstimates_ModelValidation_",
                                                        VertGroups[v], "_", PercSppModelled[c], ".RData")))[[t]][[b]]
          training_fold_spp<-droplevels(training_fold_spp)
          
          # Get the species included in the training-fold:
          full_assemblage_data<-fread(paste0(file_list[[v]], csv_scale[[j]]))
          full_assemblage_data$Cell_Id<-as.character(full_assemblage_data$Cell_Id)
          training_assemblage_data<-full_assemblage_data[Binomial %in% as.character(training_fold_spp$Binomial)]
          training_assemblage_data<-droplevels(training_assemblage_data)
          
          # Get the species included in the validation-fold:
          validation_fold_spp<-raw_trait_data[Binomial %ni% as.character(training_fold_spp$Binomial)]
          species_left_out<-validation_fold_spp[validation_fold_spp$YearOfDescription<=TimePeriod[t],] # to be used later
          validation_fold_spp<-validation_fold_spp[validation_fold_spp$YearOfDescription>TimePeriod[t],]
          validation_fold_spp<-droplevels(validation_fold_spp); species_left_out<-droplevels(species_left_out)
          validation_assemblage_data<-full_assemblage_data[Binomial %in% as.character(validation_fold_spp$Binomial)]
          validation_assemblage_data<-droplevels(validation_assemblage_data)
          rm(validation_fold_spp)
          
          # Add the species-level discovery metrics to the assemblage data and compute discovery metrics:
          assemblage_with_traits<-merge(training_assemblage_data, training_fold_spp, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
          assemblage_data<-assemblage_with_traits[, .((1 - geometric.mean(DiscProb, na.rm=T)), # PropUnknown per assemblage
                                                      (1 - geometric.mean(DiscProb_U95, na.rm=T)), # PropUnknown_L95, note that the upper and lower bounds are now inverted
                                                      (1 - geometric.mean(DiscProb_L95, na.rm=T))), # PropUnknown_U95, note that the upper and lower bounds are now inverted
                                                  by = .(Cell_Id)] 
          names(assemblage_data)<-c('Cell_Id', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95')
          rm(assemblage_with_traits)
          
          # Obs.: Richness need to be recomputed for all species described (without discard early descriptions) to avoid confounding
          # factors related to temporally varying per-assemblage richness:
          spp_left_out_ass_data<-full_assemblage_data[Binomial %in% as.character(species_left_out$Binomial)]
          all_species_but_validation<-rbind(training_assemblage_data, spp_left_out_ass_data)
          all_species_but_validation<-droplevels(all_species_but_validation)
          all_species_but_validation<-all_species_but_validation[,.(.N), by = .(Cell_Id)]
          names(all_species_but_validation)<-c('Cell_Id','Richness')
          assemblage_data<-merge(assemblage_data, all_species_but_validation, by='Cell_Id', all.x=TRUE, allow.cartesian=TRUE)
          rm(all_species_but_validation, spp_left_out_ass_data, training_assemblage_data, species_left_out)
          
          # Identify those grid cells with at least 5 species (observed):
          cells_to_keep<-assemblage_data[which(assemblage_data$Richness>=5),1]
          assemblage_data<-assemblage_data[Cell_Id %in% as.character(cells_to_keep$Cell_Id)] # remove those assemblages with < 5 species
          
          # Compute the estimated UnknownSR for each assemblage:
          assemblage_data$Est_UnknownSR<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown)) - assemblage_data$Richness
          assemblage_data$Est_UnknownSR_L95<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown_L95)) - assemblage_data$Richness
          assemblage_data$Est_UnknownSR_U95<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown_U95)) - assemblage_data$Richness
          
          # Compute the observerd discoveries based on the validation-fold:
          spp_data_not_described_yet<-validation_assemblage_data[,.(.N), by = .(Cell_Id)]
          names(spp_data_not_described_yet)<-c('Cell_Id','Obs_UnknownSR')
          
          # Register the time period, CV partition size, k-fold number, and taxonomic rank level:
          assemblage_data<-merge(assemblage_data, spp_data_not_described_yet, by='Cell_Id', all.x=TRUE, allow.cartesian=TRUE)
          rm(spp_data_not_described_yet)
          
          # Register the subsampling level, the spatial resolution, and store the results in a list:
          assemblage_data$Subsampling<-"Obs"
          assemblage_data$SpatialScale<-spatial_scale[j]
          assemblage_data$TimePeriod<-TimePeriod[t]
          assemblage_data$CV_PartitionSize<-PercSppModelled[c]
          assemblage_data$KFold<-b
          
          # Store the output for observed assemblages in a list:
          if(j==1){SubsamplingLevel_220km_list[[(length(sub_list[[j]])+1)]]<-assemblage_data}
          if(j==2){SubsamplingLevel_440km_list[[(length(sub_list[[j]])+1)]]<-assemblage_data}
          if(j==3){SubsamplingLevel_880km_list[[(length(sub_list[[j]])+1)]]<-assemblage_data}
          rm(assemblage_data)
          
          # Get the estimates of discovery potential across all levels of the subsampled assemblages: 
          for(k in 1:length(sub_list[[j]])){ # k = sub_list value giving the number of pseudoreplications
            
            # Specify the filename of the subsampled assemblage data and load it:
            myfile<-file.path(paste0(subsampled_list[[v]], spatial_scale[[j]], "_", sub_list[[j]][k], ".Rdata"))
            load(myfile) # load the subsampled assemblage data
            
            # Get the training-data for the subsampled assemblages:
            ssam$Cell_Id<-as.character(ssam$Cell_Id)
            training_ssam_data<-ssam[Binomial %in% as.character(training_fold_spp$Binomial)]
            training_ssam_data<-droplevels(training_ssam_data)
            
            # Get the species included in the validation-fold:
            validation_fold_spp<-raw_trait_data[Binomial %ni% as.character(training_fold_spp$Binomial)]
            species_left_out<-validation_fold_spp[validation_fold_spp$YearOfDescription<=TimePeriod[t],] # to be used later
            validation_fold_spp<-validation_fold_spp[validation_fold_spp$YearOfDescription>TimePeriod[t],]
            validation_fold_spp<-droplevels(validation_fold_spp); species_left_out<-droplevels(species_left_out)
            validation_ssam_data<-ssam[Binomial %in% as.character(validation_fold_spp$Binomial)]
            validation_ssam_data<-droplevels(validation_ssam_data)
            rm(validation_fold_spp)
            
            # Add the species-level discovery metrics to the assemblage data and compute discovery metrics:
            assemblage_with_traits<-merge(training_ssam_data, training_fold_spp, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
            ssam_data<-assemblage_with_traits[, .((1 - geometric.mean(DiscProb, na.rm=T)), # PropUnknown per assemblage
                                                  (1 - geometric.mean(DiscProb_U95, na.rm=T)), # PropUnknown_L95, note that the upper and lower bounds are now inverted
                                                  (1 - geometric.mean(DiscProb_L95, na.rm=T))), # PropUnknown_U95, note that the upper and lower bounds are now inverted
                                                    by = .(Iter, Cell_Id)] 
            names(ssam_data)<-c('Iter', 'Cell_Id', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95')
            rm(assemblage_with_traits)
            
            # Average subsampled assemblage metrics across iterations:
            assemblage_data<-ssam_data[, .(mean(PropUnknown, na.rm=T), # PropUnknown per assemblage
                                   mean(PropUnknown_L95, na.rm=T), # PropUnknown_L95, note that the upper and lower bounds are now inverted
                                   mean(PropUnknown_U95, na.rm=T)), # PropUnknown_U95, note that the upper and lower bounds are now inverted
                                              by = .(Cell_Id)] 
            names(assemblage_data)<-c('Cell_Id', 'PropUnknown', 'PropUnknown_L95', 'PropUnknown_U95')
            
            # Obs.: Richness need to be recomputed for all species described (without discard early descriptions) to avoid confounding
            # factors related to temporally varying per-assemblage richness:
            spp_left_out_ass_data<-ssam[Binomial %in% as.character(species_left_out$Binomial)]
            all_species_but_validation<-rbind(training_ssam_data, spp_left_out_ass_data)
            all_species_but_validation<-droplevels(all_species_but_validation)
            all_species_but_validation<-all_species_but_validation[,.(.N), by = .(Iter, Cell_Id)]
            names(all_species_but_validation)<-c('Iter','Cell_Id','Richness')
            all_species_but_validation<-all_species_but_validation[,.(mean(Richness, na.rm=T)), by = .(Cell_Id)]
            names(all_species_but_validation)<-c('Cell_Id','Richness')
            assemblage_data<-merge(assemblage_data, all_species_but_validation, by='Cell_Id', all.x=TRUE, allow.cartesian=TRUE)
            rm(all_species_but_validation, spp_left_out_ass_data, training_ssam_data, species_left_out)
            
            # Identify those grid cells with at least 5 (observed, not subsampled) species:
            assemblage_data<-assemblage_data[Cell_Id %in% as.character(cells_to_keep$Cell_Id)] # remove those assemblages with < 5 species
            
            # Compute the estimated UnknownSR for each assemblage:
            assemblage_data$Est_UnknownSR<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown)) - assemblage_data$Richness
            assemblage_data$Est_UnknownSR_L95<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown_L95)) - assemblage_data$Richness
            assemblage_data$Est_UnknownSR_U95<- (assemblage_data$Richness/(1-assemblage_data$PropUnknown_U95)) - assemblage_data$Richness
            
            # Compute the observerd discoveries based on the validation-fold:
            spp_data_not_described_yet<-validation_ssam_data[,.(.N), by = .(Iter, Cell_Id)]
            names(spp_data_not_described_yet)<-c('Iter','Cell_Id','Obs_UnknownSR')
            spp_data_not_described_yet<-spp_data_not_described_yet[,.(mean(Obs_UnknownSR, na.rm=T)), by = .(Cell_Id)]
            names(spp_data_not_described_yet)<-c('Cell_Id','Obs_UnknownSR')
            
            # Register the time period, CV partition size, k-fold number, and taxonomic rank level:
            assemblage_data<-merge(assemblage_data, spp_data_not_described_yet, by='Cell_Id', all.x=TRUE, allow.cartesian=TRUE)
            rm(spp_data_not_described_yet)
            
            # Register the subsampling level, the spatial resolution, and store the results in a list:
            assemblage_data$Subsampling<-sub_list[[j]][k]
            assemblage_data$SpatialScale<-spatial_scale[j]
            assemblage_data$TimePeriod<-TimePeriod[t]
            assemblage_data$CV_PartitionSize<-PercSppModelled[c]
            assemblage_data$KFold<-b
            
            # Store the output in a list:
            if(j==1){SubsamplingLevel_220km_list[[k]]<-assemblage_data}
            if(j==2){SubsamplingLevel_440km_list[[k]]<-assemblage_data}
            if(j==3){SubsamplingLevel_880km_list[[k]]<-assemblage_data}
            rm(assemblage_data)
            
          } # end of k loop (subsampling levels)
          rm(cells_to_keep, full_assemblage_data)
          
          # Store the results in a list:
          if(j==1){KFold_220km_list[[b]]<-rbind(data.table::rbindlist(SubsamplingLevel_220km_list))}
          if(j==2){KFold_440km_list[[b]]<-rbind(data.table::rbindlist(SubsamplingLevel_440km_list))}
          if(j==3){KFold_880km_list[[b]]<-rbind(data.table::rbindlist(SubsamplingLevel_880km_list))}
          rm(SubsamplingLevel_220km_list, SubsamplingLevel_440km_list, SubsamplingLevel_880km_list)
          
        } # end of j for loop across spatial resolutions
        
      } # end of b for loop (kfolds)
      
      # Unlist the output to get one single data.table with all subsampling levels and k-fold partitons:
      TimePeriod_220km_list[[t]]<-data.table::rbindlist(KFold_220km_list)
      TimePeriod_440km_list[[t]]<-data.table::rbindlist(KFold_440km_list)
      TimePeriod_880km_list[[t]]<-data.table::rbindlist(KFold_880km_list)
      rm(KFold_220km_list, KFold_440km_list, KFold_880km_list)
      
    } # end of t for loop (time periods)
    
    # Unlist the output to get one single data.table with all subsampling levels, k-fold partitions, and time periods:
    TimePeriod_220km_list<-data.table::rbindlist(TimePeriod_220km_list)
    TimePeriod_440km_list<-data.table::rbindlist(TimePeriod_440km_list)
    TimePeriod_880km_list<-data.table::rbindlist(TimePeriod_880km_list)
    
    # Export outputs:
    save(TimePeriod_220km_list, file=paste0("ModelValidation/AssemblageLevelEstimates_220km_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(TimePeriod_440km_list, file=paste0("ModelValidation/AssemblageLevelEstimates_440km_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    save(TimePeriod_880km_list, file=paste0("ModelValidation/AssemblageLevelEstimates_880km_", VertGroups[v], "_", PercSppModelled[c], ".RData"))
    rm(TimePeriod_220km_list, TimePeriod_440km_list, TimePeriod_880km_list)
    
  } # end of c for loop (cross-validation partition sizes)
  
} # end of v for loop (vertebrate groups)

# Load the RData object for the vertebrate group v and time period t (just for output visualization):
FinalOutput<-loadRData(file.path(paste0("ModelValidation/AssemblageLevelEstimates_880km_", VertGroups[v], "_", PercSppModelled[c], ".RData")))