############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com


# SUPPORTING SCRIPT 9: PLOTTING THE FIGURES OF THE MAIN TEXT
############################################################################################################################
# SUPPORTING SCRIPT 9: PLOTTING THE FIGURES OF THE MAIN TEXT

# Steps in this script:
#  1. Plot figure 1.
#  2. Plot figure 2.
#  3. Plot figure 3.
#  4. Plot figure 4.

# First, clean workspace:
rm(list=ls()); gc()
setwd("C:/Users/apena/OneDrive/Mapping Linnean Shortfall")
setwd("C:/Users/apena/Desktop/NATURE_SCRIPTS")

# Install and load R packages needed to run the analysis:
install.packages("devtools")
install_github("kassambara/easyGgplot2")
needed_packages<-c("data.table", "stringr", "plyr", "readr", "tools", "Hmisc", "viridis", "tidyverse", "ggnewscale", "MuMIn", 
                   "ggpubr", "RColorBrewer", "ggplot2", "easyGgplot2", "foreach", "cowplot", "dplyr" ,"egg", "gridExtra", 
                   "flexsurv", "rgdal", "raster", "tidyr", "maptools", "scatterpie", "lattice", "scales", "grid", "sf")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages,new.packages)

#####

# STEP 1 - PLOT FIGURE 1
##########################################################################################################################
# STEP 1 - PLOT FIGURE 1

# Figure 1 is composed of different pieces of information. 
# Part 1: it includes species-level attributes for four example species.
# Part 2: it includes the histograms of description dates per vertebrate group.
# Part 3: it includes two components: (i) cumulative number of species described and (ii) species discovery curves.

### Figure 1 - Part 1
# Load trait data and reorder levels of tetrapod classes:
SpeciesLevelData<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
SpeciesLevelData$Class<-factor(SpeciesLevelData$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Get the trait values 'line range plots'
# Create a table with the significant example species:
species_list<-c("Brachycephalus guarani", "Sceloporus adleri", "Chaetophractus vellerosus", "Dromaius novaehollandiae")
example_species<-as.data.frame(matrix(nrow=4, ncol=5))
names(example_species)<-c("Binomial", "Class", "RangeSize", "BodySize", "ActivityBioregion")

# Get the standardized value of each predictor variable (converted from 0 to 1):
range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
for(i in 1:4){
  
  # Filter species-level attributes for a given vertebrate group:
  trait_data_subset<-SpeciesLevelData[SpeciesLevelData$Class %in% levels(SpeciesLevelData$Class)[i],]
  trait_data_subset<-droplevels(trait_data_subset)
  
  # Log10 transform the predictors intended for plotting:
  trait_data_subset$LogBodySize<-scale(log10(trait_data_subset$BodySize), center=F, scale=F)
  trait_data_subset$LogRangeSize<-scale(log10(trait_data_subset$RangeSize), center=F, scale=F)
  trait_data_subset$LogActivityBioregion<-scale(log10(trait_data_subset$ActivityBioregion+1), center=F, scale=F)

  # Standardize predictors to vary from 0 to 1:
  trait_data_subset$StdBodySize<-range01(trait_data_subset$LogBodySize)
  trait_data_subset$StdRangeSize<-range01(trait_data_subset$LogRangeSize)
  trait_data_subset$StdActivityBioregion<-range01(trait_data_subset$LogActivityBioregion)
  
  # Store information on the example species:
  example_species[i,1]<-species_list[i] # Binomial
  example_species[i,2]<-levels(SpeciesLevelData$Class)[i] # Taxonomic class
  example_species[i,3]<-trait_data_subset[trait_data_subset$Binomial==species_list[i],21] # StdBodySize
  example_species[i,4]<-trait_data_subset[trait_data_subset$Binomial==species_list[i],22] # StdRangeSize
  example_species[i,5]<-trait_data_subset[trait_data_subset$Binomial==species_list[i],23] # StdActivityBioregion
  rm(trait_data_subset)
}

# Filter the most illustrative  variables for the species set selected:
example_species<-reshape2::melt(example_species, id.vars= c("Class", "Binomial"), variable.name="Label", value.name="Value")
example_species$Label<-factor(example_species$Label,
                              levels = c("ActivityBioregion", "BodySize", "RangeSize"),
                              labels = c("Taxonomic activity", "Body size", "Geographic\n range size"))

# Plot the species-level attributes for example species:
myRangePlots<-list()
myColors <- c('violetred2', 'green3', 'orange',  'slateblue2')
for(i in 1:4){
  species_example_filtered<-example_species[example_species$Binomial==species_list[i],]
  myRangePlots[[i]]<-ggplot(data=species_example_filtered, aes(x=Label, y=Value, ymin=Value, ymax=Value)) +
    geom_linerange(aes(ymin=0, ymax=1), cex=0.5, color="gray50", linetype="dotted") +
    geom_pointrange(colour=myColors[i], fill=myColors[i], shape=108, size=2) + 
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("") + ylab("Standardized attribute value") +
    scale_y_continuous(limits=c(0,1), breaks = c(0, 0.5, 1)) +
    theme(panel.grid.minor = element_blank(), # remove minor gridlines
          panel.grid.major = element_blank(), # remove major gridlines
          panel.background = element_blank(), # white background
          axis.title = element_text(size=12, face="bold", colour="black"), # axis title aesthetics
          axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)), # margin between axis.title and axis.values
          axis.title.y = element_text(margin=margin(t=0, r=0, b=0, l=0)), # margin between axis.title and axis.values
          
          axis.text.x = element_text(size=12, face="bold", colour="black"), # axis text aesthetics
          axis.text.y = element_blank(), # axis text aesthetics
          
          axis.line.x = element_line(colour="black"), # axis lines aesthetitcs
          axis.ticks.y=element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y = element_blank(),
          #strip.text.y = element_text(hjust=1, vjust=0.25, angle=180, size = 10),
          legend.position="none",
          plot.background=element_rect(fill="transparent", colour=NA)) +
    guides(col=guide_legend(ncol=1, byrow=FALSE, reverse=T), fill=guide_legend(ncol=1, byrow=FALSE, reverse=T))
}

# Plot the part 1 of Figure 1, and export for later use:
Multipanel_plot_1<-ggpubr::ggarrange(myRangePlots[[1]] + xlab("Example\nspecies"), myRangePlots[[2]], myRangePlots[[3]], myRangePlots[[4]],
                                     labels=c("A", "", "", ""), font.label=list(size=12, color = "black"), ncol=4, nrow=1) +
  annotate("text", label=bquote('Body size'), x=0.15, y=0.89, color="black", size=3, hjust=0.5) +
  annotate("text", label=bquote('Geographic range size'), x=0.15, y=0.705, color="black", size=3, hjust=0.5) +
  annotate("text", label=bquote('Taxonomic activity'), x=0.123, y=0.51, color="black", size=3, hjust=0.5); Multipanel_plot_1
ggsave("Figure1_Part1.pdf", plot=Multipanel_plot_1, width=10, height=1.5, units="in", bg = "transparent")



### Figure 1 - Part 2
# Load trait data and reorder levels of tetrapod classes:
SpeciesLevelData<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
SpeciesLevelData$Class<-factor(SpeciesLevelData$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Get the median year of description for each tetrapod class:
SpeciesLevelData %>%
  dplyr::group_by(Class) %>% 
  dplyr::summarise(Richness=median(YearOfDescription, na.rm=T))

# Plot histograms of description dates for each vertebrate group:
myPlot<-list()
mySolidColors <- c('violetred2', 'green3', 'orange',  'slateblue2')
for(i in 1:4){
  trait_data_subset<-SpeciesLevelData[Class %in% levels(SpeciesLevelData$Class)[i]]
  myPlot[[i]]<-ggplot2.histogram(data=trait_data_subset, xName='YearOfDescription', groupName='Class',
                                 addMeanLine=TRUE,
                                 meanLineSize=0.75,
                                 meanLineColor=mySolidColors[i],
                                 groupColors=mySolidColors[[i]],
                                 alpha=0.5,
                                 addDensityCurve = TRUE,
                                 removePanelGrid=TRUE,
                                 removePanelBorder=TRUE,
                                 axisLine=c(0.5, "solid", "black"),
                                 showLegend=TRUE,
                                 binwidth=5,
                                 scale="density",
                                 backgroundColor="white") +
    labs(x="", y="") +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=10, colour="black", face="plain"),
          axis.title.x=element_text(size=12, colour="black", face="plain"),
          axis.title.y=element_text(size=12, colour="black", face="plain"),
          plot.margin = rep(unit(0,"null"),4),
          panel.spacing = unit(0,"null"),
          panel.background=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.background=element_rect(fill="transparent", colour=NA)) +
    scale_x_continuous(limits = c(1755,2014)) +
    scale_y_continuous(limits = c(0,0.017))
}



### Figure 1 - Part 3
# Load trait data and reorder levels of tetrapod classes:
SpeciesLevelData<-fread("SpeciesLevelPredictorsSubset.csv", stringsAsFactors=TRUE, encoding="UTF-8")
SpeciesLevelData$Class<-factor(SpeciesLevelData$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Get the species accumulation curves for each group:
ecdf_vertebrates<-as.data.frame(matrix(nrow=1, ncol=3)) 
colnames(ecdf_vertebrates)<-c("YearOfDescription", "PercentageDescribed", "Class")
for(i in 1:4){
  trait_data_subset<-SpeciesLevelData[Class==levels(SpeciesLevelData$Class)[i]]
  n<-sum(!is.na(trait_data_subset$YearOfDescription))
  cumulative_data<-as.data.frame(cbind(YearOfDescription=sort(trait_data_subset$YearOfDescription), PercentageDescribed=(1:n)/n))
  cumulative_data$Class<-levels(SpeciesLevelData$Class)[i]
  ecdf_vertebrates<-rbind(ecdf_vertebrates, cumulative_data)
  rm(cumulative_data, n, trait_data_subset)
}
ecdf_vertebrates<-ecdf_vertebrates[-1,]

# Build discovery curves for each species using the avg. weigthed coefs.
# Note that this is only for visualization purposes, since computations to build the average weighted discovery curve are truly gigantic.
myfiles<-paste("ModelAveraging/AvgModel_Coefs*.csv", sep="")
All_Avg_Coefs<-lapply(Sys.glob(myfiles), fread, h=T)
All_Avg_Coefs<-rbindlist(All_Avg_Coefs)

# For each vertebrate group, estimate the species discovery curves:
for (i in 1:4){
  
  # Get the model coefs and raw data for the vertebrate group i:
  model_coefs<-All_Avg_Coefs[All_Avg_Coefs$Class==levels(SpeciesLevelData$Class)[i],]
  model_coefs<-as.data.frame(model_coefs)
  trait_data<-SpeciesLevelData[Class==levels(SpeciesLevelData$Class)[i]]
  
  # Filter the species data to comprise those described between 1759 and 2014:
  trait_data<-trait_data[trait_data$YearOfDescription>=1759 & trait_data$YearOfDescription<=2014,]
  
  # Load the model coefficients, including information on the best family error distribution:
  my_csv_file<-paste0("ModelAveraging/set_of_models_", unique(trait_data$Class), ".csv")
  set_of_models<-read.csv(my_csv_file)
  
  # Prepare the response variable for time-to-event models:
  trait_data$Time<-(trait_data$YearOfDescription-1758)/(2018-1758) # time to the discovery event
  trait_data$Censor<-1 # Censor variable (it informs if the event happened)
  
  # Log10 transform and standardize predictors (mean 0 and SD = 1):
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
  trait_data<-as.data.frame(trait_data)
  
  # Predict the survival curves for each species based on the full model coefs:
  pct<-seq(0, 1, by=0.02) # this is a coarser resolution for the discovery curve (visualization purposes only)
  ProbDisc_per_spp<-as.data.frame(matrix(nrow=length(pct), ncol=nrow(trait_data)))
  
  if(unique(set_of_models$FamilyError)=="exp"){
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pexp(pct, rate=exp(as.numeric(model_coefs[12,1] +
                                                                trait_data$LogBodySize[k]*model_coefs[1,1] +
                                                                trait_data$LogTaxPerClade[k]*model_coefs[2,1] +
                                                                trait_data$LogRangeSize[k]*model_coefs[3,1] +
                                                                trait_data$LogAMT[k]*model_coefs[4,1] +
                                                                trait_data$LogAPP[k]*model_coefs[5,1] +
                                                                trait_data$LogTS[k]*model_coefs[6,1] +
                                                                trait_data$LogPS[k]*model_coefs[7,1] +
                                                                trait_data$LogElevM[k]*model_coefs[8,1] +
                                                                trait_data$LogPopD[k]*model_coefs[9,1] +
                                                                trait_data$LogTaxPerBiome[k]*model_coefs[10,1] +
                                                                trait_data$LogRarity[k]*model_coefs[11,1])))}

  } # species-level predictions exponential model
  
  if(unique(set_of_models$FamilyError)=="weibull"){
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pweibull(pct, shape=model_coefs[12,1], 
                                         scale=as.numeric(model_coefs[13,1]*
                                                            exp(trait_data$LogBodySize[k]*model_coefs[1,1] +
                                                                  trait_data$LogTaxPerClade[k]*model_coefs[2,1] +
                                                                  trait_data$LogRangeSize[k]*model_coefs[3,1] +
                                                                  trait_data$LogAMT[k]*model_coefs[4,1] +
                                                                  trait_data$LogAPP[k]*model_coefs[5,1] +
                                                                  trait_data$LogTS[k]*model_coefs[6,1] +
                                                                  trait_data$LogPS[k]*model_coefs[7,1] +
                                                                  trait_data$LogElevM[k]*model_coefs[8,1] +
                                                                  trait_data$LogPopD[k]*model_coefs[9,1] +
                                                                  trait_data$LogTaxPerBiome[k]*model_coefs[10,1] +
                                                                  trait_data$LogRarity[k]*model_coefs[11,1])))}
    }
  
  if(unique(set_of_models$FamilyError)=="lnorm"){
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                plnorm(pct, meanlog=as.numeric(model_coefs[12,1]*
                                                                 exp(trait_data$LogBodySize[k]*model_coefs[1,1] +
                                                                       trait_data$LogTaxPerClade[k]*model_coefs[2,1] +
                                                                       trait_data$LogRangeSize[k]*model_coefs[3,1] +
                                                                       trait_data$LogAMT[k]*model_coefs[4,1] +
                                                                       trait_data$LogAPP[k]*model_coefs[5,1] +
                                                                       trait_data$LogTS[k]*model_coefs[6,1] +
                                                                       trait_data$LogPS[k]*model_coefs[7,1] +
                                                                       trait_data$LogElevM[k]*model_coefs[8,1] +
                                                                       trait_data$LogPopD[k]*model_coefs[9,1] +
                                                                       trait_data$LogTaxPerBiome[k]*model_coefs[10,1] +
                                                                       trait_data$LogRarity[k]*model_coefs[11,1],
                                                                     sdlog=model_coefs[13,1])))}
    }
  
  if(unique(set_of_models$FamilyError)=="llogis"){
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pllogis(pct, shape=model_coefs[12,1], 
                                        scale=as.numeric(model_coefs[13,1]*
                                                           exp(trait_data$LogBodySize[k]*model_coefs[1,1] +
                                                                 trait_data$LogTaxPerClade[k]*model_coefs[2,1] +
                                                                 trait_data$LogRangeSize[k]*model_coefs[3,1] +
                                                                 trait_data$LogAMT[k]*model_coefs[4,1] +
                                                                 trait_data$LogAPP[k]*model_coefs[5,1] +
                                                                 trait_data$LogTS[k]*model_coefs[6,1] +
                                                                 trait_data$LogPS[k]*model_coefs[7,1] +
                                                                 trait_data$LogElevM[k]*model_coefs[8,1] +
                                                                 trait_data$LogPopD[k]*model_coefs[9,1] +
                                                                 trait_data$LogTaxPerBiome[k]*model_coefs[10,1] +
                                                                 trait_data$LogRarity[k]*model_coefs[11,1])))}
    }
  
  if(unique(set_of_models$FamilyError)=="gamma"){
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                              .combine = 'cbind',
                              .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                pgamma(pct, shape=model_coefs[12,1], 
                                       rate=as.numeric(model_coefs[13,1]*
                                                         exp(trait_data$LogBodySize[k]*model_coefs[1,1] +
                                                               trait_data$LogTaxPerClade[k]*model_coefs[2,1] +
                                                               trait_data$LogRangeSize[k]*model_coefs[3,1] +
                                                               trait_data$LogAMT[k]*model_coefs[4,1] +
                                                               trait_data$LogAPP[k]*model_coefs[5,1] +
                                                               trait_data$LogTS[k]*model_coefs[6,1] +
                                                               trait_data$LogPS[k]*model_coefs[7,1] +
                                                               trait_data$LogElevM[k]*model_coefs[8,1] +
                                                               trait_data$LogPopD[k]*model_coefs[9,1] +
                                                               trait_data$LogTaxPerBiome[k]*model_coefs[10,1] +
                                                               trait_data$LogRarity[k]*model_coefs[11,1])))}
    }
  
  if(unique(set_of_models$FamilyError)=="gompertz"){
    ProbDisc_per_spp<-foreach(k = 1:nrow(trait_data),
                                .combine = 'cbind',
                                .packages = c("flexsurv"))  %dopar% { # cbind the results of foreach
                                  pgompertz(pct, shape=model_coefs[12,1],
                                            rate=as.numeric(model_coefs[13,1]*
                                                              exp(trait_data$LogBodySize[k]*model_coefs[1,1] +
                                                                    trait_data$LogTaxPerClade[k]*model_coefs[2,1] +
                                                                    trait_data$LogRangeSize[k]*model_coefs[3,1] +
                                                                    trait_data$LogAMT[k]*model_coefs[4,1] +
                                                                    trait_data$LogAPP[k]*model_coefs[5,1] +
                                                                    trait_data$LogTS[k]*model_coefs[6,1] +
                                                                    trait_data$LogPS[k]*model_coefs[7,1] +
                                                                    trait_data$LogElevM[k]*model_coefs[8,1] +
                                                                    trait_data$LogPopD[k]*model_coefs[9,1] +
                                                                    trait_data$LogTaxPerBiome[k]*model_coefs[10,1] +
                                                                    trait_data$LogRarity[k]*model_coefs[11,1])))}
    }
  
  # colnames = species & # rownames = estimated year of discovery scaled from 0 to 1 (0 = 1758; 1 = 2018):
  ProbDisc_per_spp<-as.data.frame(ProbDisc_per_spp)
  names(ProbDisc_per_spp)<-trait_data$Binomial;   rownames(ProbDisc_per_spp)<-pct 
  
  # Add the year of discovery as a new column in the ProbDisc object
  ProbDisc_per_spp$YearPred<-as.numeric(as.character(row.names(ProbDisc_per_spp)))
  ProbDisc_per_reshaped<-reshape2::melt(ProbDisc_per_spp, id.vars= c("YearPred"), variable.name="Species", value.name="ProbDisc") # prepare for plot
  ProbDisc_per_reshaped<-ProbDisc_per_reshaped[,c(3,2,1)]
  Predictions_per_spp<-as.data.table(ProbDisc_per_reshaped)
  rm(ProbDisc_per_reshaped, ProbDisc_per_spp)
  
  # Join table with traits:
  trait_data<-as.data.table(trait_data)
  trait_data<-trait_data[, .(Binomial, Class, YearOfDescription)]
  Predictions_per_spp<-merge(Predictions_per_spp, trait_data, by.x="Species", by.y="Binomial", all.x=T)
  
  # Store discovery curves for each group in separate objects:
  if(i==1){Predictions_per_spp_amph<-Predictions_per_spp}
  if(i==2){Predictions_per_spp_rept<-Predictions_per_spp}
  if(i==3){Predictions_per_spp_mamm<-Predictions_per_spp}
  if(i==4){Predictions_per_spp_aves<-Predictions_per_spp}
  rm(Predictions_per_spp, trait_data, set_of_models, model_coefs)
  
} # end of i for loop (across vertebrate groups)

# Merge the estimates of all vertebrates:
Predictions_per_spp<-rbind(Predictions_per_spp_amph, Predictions_per_spp_rept, Predictions_per_spp_mamm, Predictions_per_spp_aves)
rm(Predictions_per_spp_amph, Predictions_per_spp_rept, Predictions_per_spp_mamm, Predictions_per_spp_aves)

# Set the correct order of levels for the taxonomic class variable:
Predictions_per_spp$Class<-factor(Predictions_per_spp$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))
Predictions_per_spp<-droplevels(Predictions_per_spp)

# Create pre-defined objects storing the colours that will be used for plotting:
mySolidColors<-c('violetred2', 'green3', 'orange', 'slateblue2')
myLightColors<-c('#fde0ef', '#e6f5d0', '#fee0b6', '#d8daeb')

# Recall the binomial name of the example species:
example_species<-c("Brachycephalus guarani", "Sceloporus adleri", "Chaetophractus vellerosus", "Dromaius novaehollandiae")

# Plot the species discovery curves for each vertebrate group:
Survival_plots<-list()
for (i in 1:4){
  
  # Filter one vertebrate group and put species in the order of their discovery probability:
  ecdf_vertebrates_subset<-ecdf_vertebrates[ecdf_vertebrates$Class==levels(Predictions_per_spp$Class)[i],]
  Predictions_per_spp_subset<-Predictions_per_spp[Class %in% as.character(levels(Predictions_per_spp$Class)[i])] 
  Predictions_per_spp_subset<-droplevels(Predictions_per_spp_subset)
  Predictions_per_example_spp<-Predictions_per_spp_subset[Species %in% as.character(example_species[i])]
  
  # Get the average discovery curve:
  AvgCurve<-Predictions_per_spp_subset %>% dplyr::group_by(YearPred) %>% 
    dplyr::summarise(AvgProbDisc=mean(ProbDisc, na.rm=T))
  
  # Increasing order to plot the higher values on top of the lower values:
  Survival_plots[[i]]<-ggplot() +
    # Individual discovery curves for each species:
    geom_smooth(data=Predictions_per_spp_subset, aes(x=((YearPred*260)+1758), y=ProbDisc, group=Species), colour=myLightColors[i], span=.1, size=0.15, fill=NA) +
    
    # Add the line for the cumulative proportion of species descriptions across time (empirical pattern):
    geom_line(data=ecdf_vertebrates_subset, aes(x=YearOfDescription, y=PercentageDescribed), colour="black", size=0.75) +
    
    # Add the average discovery curve across all species (estimated pattern):
    geom_smooth(data=AvgCurve, aes(x=((YearPred*260)+1758), y=AvgProbDisc), colour=mySolidColors[i], span=.2, size=0.75, linetype="longdash", fill=NA) +
    
    # Add the discovery curve of the example species:
    geom_smooth(data=Predictions_per_example_spp, aes(x=((YearPred*260)+1758), y=ProbDisc), colour=mySolidColors[i], span=.2, size=0.15, linetype="solid", fill=NA) +
    
    labs(x="Year of discovery", y="") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour="black"),
          axis.text.x=element_text(size=10, colour="black", face="plain"),
          axis.text.y=element_text(size=10, colour="black", face="plain"),
          axis.title.x=element_text(size=10, colour="black", face="bold"),
          axis.title.y=element_text(size=10, colour="black", face="bold"),
          plot.background=element_rect(fill="white", colour=NA),
          legend.position="none") +
    #geom_vline(xintercept=2015, colour="black", size=1.2, linetype="dotted") +
    scale_x_continuous(limits = c(1759,2019), breaks=c(1800, 1900, 2000)) +
    scale_y_continuous(limits = c(0,1))
  
  rm(Predictions_per_spp_subset, AvgCurve, Predictions_per_example_spp, ecdf_vertebrates_subset)
}

# Assemblage the part 2 and 3 of Figure 1:
Multipanel_plot_2<-ggpubr::ggarrange(myPlot[[1]] + ylab("Density"), myPlot[[2]], myPlot[[3]], myPlot[[4]],
                           Survival_plots[[1]] + ylab("Probability | Proportion"), Survival_plots[[2]], Survival_plots[[3]], Survival_plots[[4]],
                           labels=c("B", "C", "D", "E", "F", "G", "H", "I"), align="hv",
                           font.label=list(size=12, color = "black"), ncol=4, nrow=2)

ggsave("Figure1_Part2_3.png", plot=Multipanel_plot_2, width=10, height=4, units="in", bg = "transparent")

# Replot the part 1 of Figure 1:
Multipanel_plot_1<-ggpubr::ggarrange(myRangePlots[[1]] + xlab("Example\nspecies"), myRangePlots[[2]], myRangePlots[[3]], myRangePlots[[4]],
                                     labels=c("A", "", "", ""), font.label=list(size=12, color = "black"), ncol=4, nrow=1) +
  annotate("text", label=bquote('Body size'), x=0.15, y=0.89, color="black", size=3, hjust=0.5) +
  annotate("text", label=bquote('Geographic range size'), x=0.15, y=0.705, color="black", size=3, hjust=0.5) +
  annotate("text", label=bquote('Taxonomic activity'), x=0.123, y=0.51, color="black", size=3, hjust=0.5); Multipanel_plot_1

# Export for later use:
ggsave("Figure1_Part1.pdf", plot=Multipanel_plot_1, width=10, height=1.5, units="in", bg = "transparent")

# As a note, 'Figure1_Part1.pdf' and 'Figure1_Part2_3.pdf' were edited externally in InkScape to reach the final layout.
# Animal silhuoetes were extracted from http://www.phylopic.org/: Dendrobates-azureus, Ardeosaurus-brevipes, Chrysocyon-brachyurus, Turdus-merula

#####

# STEP 2 - PLOT FIGURE 2
##########################################################################################################################
# STEP 2 - PLOT FIGURE 2
rm(list=ls()); gc()

# Load data on average weighted coefficients for the predictors across time (output of Sensitivity Analysis):
CoefData<-fread("SensitivityAnalysis/TimePeriodOfDiscovery/AverageWeightedCoefsAcrossTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
CoefData$Class<-factor(CoefData$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Set the order that predictors will appear in the plot:
CoefData$CoefLabel<-factor(CoefData$Label,
                             levels = c("LogBodySize", "LogRangeSize", "LogRarity", "LogAPP", "LogAMT", "LogPS", "LogTS", "LogElevM", "LogPopD", "LogTaxPerClade", "LogTaxPerBiome"),
                             labels = c("Body size", "Geographic\n range size", "Range\n Rarity", "Annual\n precipitation", "Annual mean\n temperature", "Temperature\n seasonality",
                                        "Precipitation\n seasonality", "Mean\n elevation", "Human pop.\n density", "Activity per\n family", "Activity per\n bioregion"))

# Make the plts:
Coefs_plots<-list()
for (k in 1:4){
  
  # Get the coefficient data for one vertebrate group:
  my_coef_filtered<-CoefData %>% filter(Class==levels(CoefData$Class)[k])
  my_coef_filtered<-droplevels(my_coef_filtered)
  
  # Interpolate the colour palette to get 22 different levels (= time periods):
  colourCount = length(unique(my_coef_filtered$TimePeriod))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Spectral"))
  
  my_coef_filtered$TimePeriod<-as.factor(my_coef_filtered$TimePeriod)
  
  # Create the plots:
  Coefs_plots[[k]]<-ggplot(data=my_coef_filtered, aes(x=TimePeriod, y=Avg.coef, ymin=Lower_IC, ymax=Upper_IC)) +
      geom_rect(xmin=0.5, xmax=23.5, ymin=-Inf, ymax=Inf, fill="grey95") +
      geom_pointrange(aes(colour=TimePeriod, fill=TimePeriod), shape=20, size=0.2) + 
      #geom_hline(yintercept=0, lty="solid", size=0.25, colour="gray50") +  
      geom_errorbar(aes(ymin=Lower_IC, ymax=Upper_IC, col=TimePeriod), width=0.5, cex=1) +
      facet_wrap(~CoefLabel, strip.position="left", nrow=23, scales = "free_y") +
      scale_colour_manual(values=rev(getPalette(colourCount))) +
      scale_fill_manual(values=rev(getPalette(colourCount))) +
      coord_flip() +  
      xlab("") + ylab("") +
      scale_y_continuous(breaks=c(-1.00, -0.5, 0, 0.5)) +
      theme(panel.grid.minor = element_blank(), # remove minor gridlines
            panel.grid.major = element_blank(), # remove major gridlines
            panel.background = element_blank(), # white background
            axis.title = element_text(size=12), # axis title aesthetics
            axis.line = element_line(colour="black"), # axis lines aesthetitcs
            axis.text.y = element_blank(), # axis text aesthetics
            axis.ticks.y=element_blank(),
            axis.title.y = element_text(margin=margin(t=0, r=0, b=0, l=0)), # margin between axis.title and axis.values
            axis.text.x = element_text(size=10), # axis text aesthetics
            axis.title.x = element_text(margin=margin(t=0, r=0, b=0, l=0)), # margin between axis.title and axis.values
            plot.background=element_rect(fill="transparent", colour=NA),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            legend.position="none")
  
}

# Set equal width for each panel size:
Coefs_plots[[1]]<-egg::set_panel_size(Coefs_plots[[1]] + theme(strip.text.y = element_text(hjust=1, vjust=0.25, angle=180, size = 10)), width=unit(5,"cm"), height=unit(0.9,"cm"))
Coefs_plots[[2]]<-egg::set_panel_size(Coefs_plots[[2]], width=unit(5,"cm"), height=unit(0.9,"cm"))
Coefs_plots[[3]]<-egg::set_panel_size(Coefs_plots[[3]], width=unit(5,"cm"), height=unit(0.9,"cm"))
Coefs_plots[[4]]<-egg::set_panel_size(Coefs_plots[[4]], width=unit(5,"cm"), height=unit(0.9,"cm"))

# Create the multipanel plot and explort it:
Multipanel_plot<-ggpubr::ggarrange(Coefs_plots[[1]], Coefs_plots[[2]], Coefs_plots[[3]], Coefs_plots[[4]],
                                   labels=c("", "", "", ""),
                                   font.label=list(size=12, color = "black"), ncol=4, nrow=1, align = "h") +
  annotate("text", label='Standardized Coefficient (95% CI)', x=0.5, y=0.05, color="black", fontface="bold", size=4.5, hjust=0.5); Multipanel_plot
ggsave("Figure_2.png", plot=Multipanel_plot, width=10, height=6, units="in", bg = "transparent")

# As a note, final layout of Figure 2 was reached in InkScape software.
# Animal silhuoetes were extracted from http://www.phylopic.org/: Dendrobates-azureus, Ardeosaurus-brevipes, Chrysocyon-brachyurus, Turdus-merula

#####

# STEP 3 - PLOT FIGURE 3
##########################################################################################################################
# STEP 3 - PLOT FIGURE 3

# Get the discovery metrics per genus, family, clade and reorder according to the tip.labels:
TaxonLevelPredictions<-fread("TaxonLevelAnalysis/TaxonLevelEstimates_HigherLevelGrouping_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")
TaxonLevelPredictions[TaxonLevelPredictions$Suborder=="Neob. (Hyloidea)",3]<-"Hyloidea"
TaxonLevelPredictions[TaxonLevelPredictions$Suborder=="Neob. (Ranoidea)",3]<-"Ranoidea"

# Create a empty data frame to receive the standardized values of 'PropUnknown' and 'UnknownSR':
StdTaxonLevelPredictions<-as.data.frame(matrix(nrow=1, ncol=ncol(TaxonLevelPredictions), NA))
names(StdTaxonLevelPredictions)<-names(TaxonLevelPredictions)

# Make the 'PropUnknown' to vary between 0 and 1:
for(i in 1:nlevels(TaxonLevelPredictions$Class)){ # i vertebrate groups
  
  vert_filtered<-TaxonLevelPredictions %>% filter(Class==levels(TaxonLevelPredictions$Class)[i])
  vert_filtered<-droplevels(vert_filtered)
  
  track_denominator<-(max(vert_filtered$PropUnknown, na.rm=T) - min(vert_filtered$PropUnknown, na.rm=T))
  track_min<-min(vert_filtered$PropUnknown, na.rm=T)
  vert_filtered$PropUnknown<-(vert_filtered$PropUnknown - track_min) / track_denominator
  vert_filtered$PropUnknown_L95<-(vert_filtered$PropUnknown_L95 - track_min) / track_denominator
  vert_filtered$PropUnknown_U95<-(vert_filtered$PropUnknown_U95 - track_min) / track_denominator
  
  # Bind
  StdTaxonLevelPredictions<-rbind(StdTaxonLevelPredictions, vert_filtered)
}
StdTaxonLevelPredictions<-StdTaxonLevelPredictions[-1,] # remove the first  (empty) row

# Reorder levels for better aesthetics:
StdTaxonLevelPredictions$Class<-factor(StdTaxonLevelPredictions$Class, 
                                       levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"),
                                       labels = c("Amphibians", "Reptiles", "Mammals", "Birds"))

# Make the 'UnknownSR' relative to all terrestrial vertebrates combined:
range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
track_denominator<-sum(StdTaxonLevelPredictions$UnknownSR, na.rm=T)
StdTaxonLevelPredictions$UnknownSR<-(StdTaxonLevelPredictions$UnknownSR/track_denominator)
track_denominator<-sum(StdTaxonLevelPredictions$UnknownSR_L95, na.rm=T)
StdTaxonLevelPredictions$UnknownSR_L95<-(StdTaxonLevelPredictions$UnknownSR_L95/track_denominator)
track_denominator<-sum(StdTaxonLevelPredictions$UnknownSR_U95, na.rm=T)
StdTaxonLevelPredictions$UnknownSR_U95<-(StdTaxonLevelPredictions$UnknownSR_U95/track_denominator)
StdTaxonLevelPredictions %>% dplyr::group_by(Class) %>%
  dplyr::summarize(pies=sum(UnknownSR, na.rm=T))

# Get the most important taxa for future species discovery (Top15 of each vertebrate group):
StdTaxonLevelPredictions<-StdTaxonLevelPredictions[order(StdTaxonLevelPredictions$UnknownSR, decreasing=T),]
DiscovMetrics_Amph<-StdTaxonLevelPredictions[which(StdTaxonLevelPredictions$Class=="Amphibians"),]
DiscovMetrics_Rept<-StdTaxonLevelPredictions[which(StdTaxonLevelPredictions$Class=="Reptiles"),]
DiscovMetrics_Mamm<-StdTaxonLevelPredictions[which(StdTaxonLevelPredictions$Class=="Mammals"),]
DiscovMetrics_Mamm<-DiscovMetrics_Mamm[1:15,]
DiscovMetrics_Aves<-StdTaxonLevelPredictions[which(StdTaxonLevelPredictions$Class=="Birds"),]
DiscovMetrics_Aves<-DiscovMetrics_Aves[1:15,]
StdTaxonLevelPredictions<-rbind(DiscovMetrics_Amph, DiscovMetrics_Rept, DiscovMetrics_Mamm, DiscovMetrics_Aves)
rm(DiscovMetrics_Amph, DiscovMetrics_Rept, DiscovMetrics_Mamm, DiscovMetrics_Aves)
StdTaxonLevelPredictions<-droplevels(StdTaxonLevelPredictions)

# Specify the level of plot:
data<-as.data.frame(StdTaxonLevelPredictions)
data<-droplevels(data)

# Get the minimum and maximum values:
data$UnknownSR_Min<-NA
data$UnknownSR_Max<-NA
for(j in 1:nrow(data)){
  
  # Get minimum UnknownSR
  selected_row<-as.data.frame(data[j, c(8,9,10)])
  data$UnknownSR_Min[j]<-min(selected_row, na.rm=T)
  
  # Get maximum UnknownSR
  selected_row<-as.data.frame(data[j, c(8,9,10)])
  data$UnknownSR_Max[j]<-max(selected_row, na.rm=T)
  
}

# Set a number of 'empty bars' to add at the end of each group (empty bards = graphic space):
empty_bar<-3 
to_add<-data.frame(matrix(NA, empty_bar*nlevels(data$Class), ncol(data)))
colnames(to_add)<-colnames(data)
to_add$Class<-rep(levels(data$Class), each=empty_bar)
to_add<-rbind(to_add, to_add[c(10:12),])

# Merge the datasets with the predicted values and the NA (to make empty bars):
data<-rbind(data, to_add)
data<-data %>% arrange(Class, desc(UnknownSR))
data$id<-seq(1, nrow(data))

# Get the name and the y position of each label
label_data<-data
number_of_bar<-nrow(label_data)
# Id substract 0.5 because the letter must have the angle of the center of the bars.Not extreme right(1) or extreme left (0)
angle<-90 - 360 * (label_data$id-0.5)/number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Prepare a data frame for base lines:
base_data <- data %>% 
  group_by(Class) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# Prepare a data frame for grid (scales):
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Start the circular barplot using only amphibians:
data$NewVarN<-NA # a duplicated and modified column holding UnknownSR
track_column<-which(names(data)=='UnknownSR')
data[data$Class=="Amphibians", 14]<-data[data$Class=="Amphibians", track_column] # get UnknownSR for amphibians only
data$NewVarP<-NA # a duplicated and modified column holding PropUnknown
track_column<-which(names(data)=='PropUnknown')
data[data$Class=="Amphibians", 15]<-data[data$Class=="Amphibians", track_column] # get PropUnknown for amphibians only

# Set automatic breaks for the axes and colorbars:
myBreaks<-c(ceiling(min(data$NewVarP, na.rm=T)*100)/100, 
            floor(max(data$NewVarP, na.rm=T)*100)/100)
myBreaks<-c(myBreaks[1], ((myBreaks[2]+myBreaks[1])/2), myBreaks[2])

# Plot taxon-level predictions for amphibians only:
myPlotA<-ggplot() + 
  
  # Get the graphic space properly: 
  geom_bar(data=data, aes(x=as.factor(id), y=NewVarN, fill=NewVarP), color="black", lwd=0.1, stat="identity") +
  
  # Add thin dotted line segments:
  geom_segment(data=grid_data, aes(x=0, y=0.15, xend=nrow(data), yend=0.15), colour="#BEBEBE", alpha=1, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=0, y=0.1, xend=nrow(data), yend=0.1), colour="#BEBEBE", alpha=1, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=0, y=0.05, xend=nrow(data), yend=0.05), colour="#BEBEBE", alpha=1, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
  
  # Add thin line segments between bars of each vertebrate group:
  geom_segment(data=grid_data, aes(x=end, y=0.15, xend=start, yend=0.15), colour="#BEBEBE", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=0.1, xend=start, yend=0.1), colour="#BEBEBE", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=0.05, xend=start, yend=0.05), colour="#BEBEBE", alpha=1, size=0.3, inherit.aes=FALSE) +
  
  # Add thick line segments in the inner part of the circle to indicate each vertebrate group:
  geom_segment(aes(x=0.5, y=-0.025, xend=6.5, yend=-0.025), alpha=0.5, size=4, color="#BEBEBE", inherit.aes=FALSE) +
  geom_segment(aes(x=9.5, y=-0.025, xend=19.5, yend=-0.025), alpha=0.5, size=4, color="#BEBEBE", inherit.aes=FALSE) +
  geom_segment(aes(x=22.5, y=-0.025, xend=37.5, yend=-0.025), alpha=0.5, size=4, color="#BEBEBE", inherit.aes=FALSE) +
  geom_segment(aes(x=40.5, y=-0.025, xend=55.5, yend=-0.025), alpha=0.5, size=4, color="#BEBEBE", inherit.aes=FALSE) +
  
  # Replot bars on the top of the segments plotted before:
  geom_bar(data=data, aes(x=as.factor(id), y=NewVarN, fill=NewVarP), color="black", lwd=0.1, stat="identity") +
  scale_fill_gradientn(colours=c("#F7F7F7", "#FDE0EF", "#F1B6DA", "#DE77AE", "#C51B7D", "#8E0152"), # shades of pink
                       #na.value = NA, values=c(0, 0.005, 0.025, 0.2, 0.5, 1)) +
                       na.value = NA, values=c(0, 0.1, 0.3, 0.5, 0.7, 1), breaks=myBreaks) +
  geom_text(data=base_data, aes(x = title, y = -18, label=Class), hjust=c(1,1,0,0), 
            colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  geom_text(data=label_data, aes(x=id, y=UnknownSR_Max+0.01, label=Suborder, hjust=hjust), color="black",
            fontface="plain", alpha=0.9, size=2.7, angle=label_data$angle, inherit.aes=FALSE) +
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour="black", alpha=0.8, size=0.6, inherit.aes=FALSE) +
  
  # Add the tick mark values in the y-axis:
  annotate("text", x = rep((max(data$id)+1),4), y = c(0.05, 0.10, 0.15, 0.20), label = c("5%", "10%", "15%", "20%"), color="black", size=2.75 , angle=0, fontface="bold", hjust=1) +
  annotate("text", label=bquote('Percent of total\n future discoveries'), x=59.5, y=0.125, color="black", fontface="bold", size=3, angle=95, hjust=0.5) +
  ylim(-0.2, 0.2) +
  
  # Customize the theme:
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.justification = c(0, 0),
    legend.text=element_text(size=rel(0.8), hjust=0.5),
    legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.3,"cm"),
    legend.background=element_blank(),
    legend.title=element_text(size=10),
    legend.title.align=0.5,
    legend.spacing.x = unit(0.1, 'cm'),
    legend.spacing.y = unit(0.1, 'cm'),
    panel.background=element_blank(),
    plot.background=element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0,5), "cm")) +
  coord_polar()

# Add taxon-level predictions for reptiles:
myPlotA <- myPlotA + new_scale_fill() # add new color ramp
data$NewVarN<-NA # a duplicated and modified column holding UnknownSR
track_column<-which(names(data)=='UnknownSR')
data[data$Class=="Reptiles", 14]<-data[data$Class=="Reptiles", track_column] # get UnknownSR for reptiles only
data$NewVarP<-NA # a duplicated and modified column holding PropUnknown
track_column<-which(names(data)=='PropUnknown')
data[data$Class=="Reptiles", 15]<-data[data$Class=="Reptiles", track_column] # get PropUnknown for reptiles only
myBreaks<-c(ceiling(min(data$NewVarP, na.rm=T)*100)/100, 
            floor(max(data$NewVarP, na.rm=T)*100)/100)
myBreaks<-c(myBreaks[1], ((myBreaks[2]+myBreaks[1])/2), myBreaks[2])
myPlotB<-myPlotA +
  geom_bar(data=data, aes(x=as.factor(id), y=NewVarN, fill=NewVarP), color="black", lwd=0.1, stat="identity") +
  scale_fill_gradientn(colours=c("#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419"), # shades of green
                       na.value = NA, values=c(0, 0.1, 0.3, 0.5, 0.7, 1), breaks=myBreaks)

# Add taxon-level predictions for mammals:
data$NewVarN<-NA # a duplicated and modified column holding UnknownSR
track_column<-which(names(data)=='UnknownSR')
data[data$Class=="Mammals", 14]<-data[data$Class=="Mammals", track_column] # get UnknownSR for reptiles only
data$NewVarP<-NA # a duplicated and modified column holding PropUnknown
track_column<-which(names(data)=='PropUnknown')
data[data$Class=="Mammals", 15]<-data[data$Class=="Mammals", track_column] # get PropUnknown for reptiles only
myBreaks<-c(ceiling(min(data$NewVarP, na.rm=T)*100)/100, 
            floor(max(data$NewVarP, na.rm=T)*100)/100)
myBreaks<-c(myBreaks[1], ((myBreaks[2]+myBreaks[1])/2), myBreaks[2])
myPlotB <- myPlotB + new_scale_fill()
myPlotC<-myPlotB +
  geom_bar(data=data, aes(x=as.factor(id), y=NewVarN, fill=NewVarP), color="black", lwd=0.1, stat="identity") +
  scale_fill_gradientn(colours=c("#F7F7F7", "#FEE0B6", "#FDB863", "#E08214", "#B35806", "#7F3B08"), # shades of orange
                       na.value = NA, values=c(0, 0.1, 0.3, 0.5, 0.7, 1), breaks=myBreaks)

# Add taxon-level predictions for birds:
track_column<-which(names(data)=='UnknownSR')
data[data$Class=="Birds", 14]<-data[data$Class=="Birds", track_column] # get UnknownSR for reptiles only
data$NewVarP<-NA # a duplicated and modified column holding PropUnknown
track_column<-which(names(data)=='PropUnknown')
data[data$Class=="Birds", 15]<-data[data$Class=="Birds", track_column] # get PropUnknown for reptiles only
myBreaks<-c(ceiling(min(data$NewVarP, na.rm=T)*100)/100, 
            floor(max(data$NewVarP, na.rm=T)*100)/100)
myBreaks<-c(myBreaks[1], ((myBreaks[2]+myBreaks[1])/2), myBreaks[2])
myPlotC <- myPlotC + new_scale_fill()
myPlotD<-myPlotC +
  geom_bar(data=data, aes(x=as.factor(id), y=NewVarN, fill=NewVarP), color="black", lwd=0.1, stat="identity") +
  geom_segment(data=data, aes(x=id, y=UnknownSR_Min, xend=id, yend=UnknownSR_Max), colour="black", size=0.25, inherit.aes=FALSE) +
  geom_segment(data=data, aes(x=(id-0.15), y=UnknownSR_Min, xend=(id+0.15), yend=UnknownSR_Min), colour="black", size=0.25, inherit.aes=FALSE) +
  geom_segment(data=data, aes(x=(id-0.2), y=UnknownSR_Max, xend=(id+0.2), yend=UnknownSR_Max), colour="black", size=0.25, inherit.aes=FALSE) +
  scale_fill_gradientn(colours=c("#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B"), # shades of purple
                       na.value = NA, values=c(0, 0.1, 0.3, 0.5, 0.7, 1), breaks=myBreaks)

myPlotD <- myPlotD + new_scale_fill(); myPlotD
myfile_pdf<-paste("Figure3.pdf")
ggsave(myfile_pdf, plot=myPlotD, width=6, height=6, units="in", bg = "transparent")

# As a note, Figure 3 was finished externally, using InkScape software.
# Animal silhuoetes extracted from http://www.phylopic.org/: Dendrobates-azureus, Ardeosaurus-brevipes, Chrysocyon-brachyurus, Turdus-merula

#####

# STEP 4 - PLOT FIGURE 4
##########################################################################################################################
# STEP 4 - PLOT FIGURE 4

# Figure 4 is composed of different pieces of information. 
# Part 1: it includes assemblage-based predictions of future species discovery.
# Part 2: it includes country-based predictions of future species discovery.
# Part 3: it includes a global piechart of future species discovery according to each tetrapod class.

# Before running this script, there are files that need to be downloaded:
# 1. landcover.zip: simplified shapefile with global landcover boundaries (Modified from WWF Ecoregions.
# 2. gridcells_220km.zip: shapefile of 220 x 220 km grid cells.
# 3. hatched_lines.zip: shapefile of hatched lines across the Earth.
# 4. world_limit.zip: shapefile of world limits (boundary box).
# 5. countries.zip: shapefile of country boundaries.

# After downloading, please, define the directory where the zipped files are located:
zip_files_dir<-getwd() #  specify it mannually if different from the working directory

# Create the directory 'Shapefiles' and unzip the shapefiles:
dir.create("Shapefiles", showWarnings = F)
utils::unzip(paste0(zip_files_dir, "/landcover.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/gridcells_220km.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/hatched_lines.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/world_limit.zip"), exdir="Shapefiles")
utils::unzip(paste0(zip_files_dir, "/countries.zip"), exdir="Shapefiles")

### Figure 4 - Part 1
# Load all shapefiles that will be used for part 1 of Figure 4:
terr_cover<-readOGR(dsn='Shapefiles',layer='landcover')
grid_cells<-readOGR(dsn='Shapefiles',layer='gridcells_220km')
hatched_lines<-readOGR(dsn='Shapefiles',layer='hatched_lines')
world_limit<-readOGR(dsn='Shapefiles',layer='world_limit')

# Set the coordinate reference systems that will be used:
wgs84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # WGS84
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection

# Load the assemblage-level estimates of discovery potential at 220 km of spatial resolution:
Vert_220km<-fread("AssemblageLevelAnalysis/AssemblageLevelEstimates_220km.csv", stringsAsFactors=TRUE, encoding="UTF-8")
Vert_220km<-Vert_220km[Vert_220km$Class=="Vertebrates",] # filter for the vertebrates all combined
Vert_220km<-droplevels(Vert_220km) # remove unused levels

# Combine assemblage-level estimates with the attribute table of grid_cells object:
grid_cells@data$n_order<-c(1:nrow(grid_cells@data))
grid_cells@data$Cell_Id220<-as.character(grid_cells@data$Cell_Id220)
Vert_220km$Cell_Id220<-as.character(Vert_220km$Cell_Id220)
grid_cells@data<-merge(grid_cells@data, Vert_220km, by.x="Cell_Id220", by.y="Cell_Id", all.x=T)
grid_cells@data<-grid_cells@data[order(grid_cells@data$n_order),]

# Remove grid cells without estimates (i.e., regions without species, oceans, etc):
grid_cells<-grid_cells[-which(is.na(grid_cells@data$UnknownSR)),] # remove grid cells without species.

# Create separate layers to hold the grid cells in the top5% and top10%:
top10<-grid_cells[grid_cells@data$UnknownSR>=(quantile(grid_cells@data$UnknownSR, probs=0.9, na.rm=T)),1]
top5<-grid_cells[grid_cells@data$UnknownSR>=(quantile(grid_cells@data$UnknownSR,probs=0.95, na.rm=T)),1] 
top10@data$Priority<-1
top5@data$Priority<-1

# Top5% and Top10% grid cells need to be merged whenever possible to remove inner countours:
r<-raster() # create a empty raster to receive data
res(r)<-0.2 # desired raster resolution in decimal degrees
top10_raster<-rasterize(top10, r, field="Priority", background=NA)
top5_raster<-rasterize(top5, r, field="Priority", background=NA)
top10_polygons<-rasterToPolygons(top10_raster, dissolve=T)
top5_polygons<-rasterToPolygons(top5_raster, dissolve=T)
rm(top10, top10_raster, top5, top5_raster)

# Convert the coordinate reference system of spatial objects to an equal area projection:
top10_polygons<-spTransform(top10_polygons, equalareaproj)
top5_polygons<-spTransform(top5_polygons, equalareaproj)
crs(hatched_lines)<-crs(top5_polygons) # to secure identical crs
hatched_top5<-raster::intersect(hatched_lines, top5_polygons) # get the spatial intersection between

# Prepare spatial objects for use in ggplot2:
top10_polygons@data$id<-rownames(top10_polygons@data)
top10_polygons_df<-fortify(top10_polygons, region="id")
top5_polygons@data$id<-rownames(top5_polygons@data)
top5_polygons_df<-fortify(top5_polygons, region="id")
hatched_top5@data$id<-rownames(hatched_top5@data)
hatched_top5_df<-fortify(hatched_top5, region="id")

# Same as above but for the 220 km grid cells
grid_cells<-spTransform(grid_cells, equalareaproj)
grid_cells@data$id<-rownames(grid_cells@data)
grid_cell_df<-fortify(grid_cells, region="id")
grid_cell_df<-join(grid_cell_df, grid_cells@data, by="id")
terr_cover@data$id<-rownames(terr_cover@data)
terr_cover_df<-fortify(terr_cover, region="id")
terr_cover_df<-join(terr_cover_df, terr_cover@data, by="id") # joining tables by "id"
world_limit<-spTransform(world_limit, equalareaproj)
world_limit_df<-fortify(world_limit)

# Plot Figure 4, part 1:
MyPlot_1<-ggplot(data=grid_cell_df) +
  geom_polygon(data=world_limit_df, aes(long, lat, group=group), colour="black", fill="lightcyan", size=0.25) + # plot world limits filled by the ocean
  geom_polygon(data=terr_cover_df, aes(long, lat, group=group), colour=NA, fill="white", size=0.25) +   # plot landcover boundaries
  geom_polygon(data=grid_cell_df, aes(long, lat, group=group, fill=UnknownSR*100)) + # plot grid cells (assemblage-based estimates)
  geom_polygon(data=top10_polygons_df, aes(long, lat, group=group), colour="black", fill=NA, size=0.125) + # add the countour lines for the top10% grid cells
  geom_polygon(data=hatched_top5_df, aes(long, lat, group=group), colour="black", fill=NA, size=0.125) + # add the hatched lines for the top5% grid cells
  geom_polygon(data=terr_cover_df, aes(long, lat, group=group), colour="gray30", fill=NA, size=0.25) + # replot landcover boundaries but without fill colour
  geom_polygon(data=world_limit_df, aes(long, lat, group=group), colour="black", fill=NA, size=0.25) + # replot world boundary box without fill colour
  coord_equal() +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.justification = c(0, 0),
        legend.position=c(0.05, 0.06), legend.direction="horizontal",
        legend.text=element_text(size=rel(0.4), hjust=0.5),
        legend.key.width=unit(0.4,"cm"), legend.key.height=unit(0.175,"cm"),
        legend.background=element_rect(fill="white", size=0.2, linetype='solid', colour='black'),
        legend.title=element_text(size=5),
        legend.title.align=0.5,
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background=element_blank(),
        plot.background=element_blank(),
        plot.margin = rep(unit(-0.03,"null"),4),
        panel.spacing = unit(0,"null")) +
  guides(fill=guide_colorbar(title="Percent of total\nfuture discoveries", title.position = "top", 
                             label=T, label.position="bottom", nbin=15, barwidth=5, 
                             draw.ulim=T, draw.llim=T, frame.colour="black", ticks=T)) +
  scale_fill_gradientn(colours=c('#ffffff', '#f0f0f0', '#d9d9d9', '#bdbdbd', '#969696', '#737373', '#525252', '#252525', '#000000'), # white to black
                       na.value = NA, values=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1), breaks=c(0.0, 0.4, 0.8)) 

ColourBarLegend_1<-ggpubr::as_ggplot((cowplot::get_legend(MyPlot_1))) # get only the legend
MyPlot1_StdLegend<- MyPlot_1 + theme(legend.position = "none") + annotation_custom(ggplotGrob(ColourBarLegend_1), xmin=-19000000, ymin=-8000000)


### Figure 4 - Part 2
# Load all shapefiles that will be used for part 2 of Figure 4 and convert them into equal area projections:
wgs84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # WGS84
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection
country_pol<-readOGR(dsn='Shapefiles',layer='gadm36_cea') # XSP and VAT dropped out
country_pol@data$n_order<-c(1:nrow(country_pol@data))
country_pol<-spTransform(country_pol, equalareaproj)
world_limit<-readOGR(dsn='Shapefiles',layer='world_limit')
world_limit<-spTransform(world_limit, equalareaproj)
world_limit_df<-fortify(world_limit)

# Load the country-level estimates of discovery potential (results are similar for outputs based on 220, 440, or 880 km of spatial resolution):
CountryLevelData<-fread("CountryLevelAnalysis/CountryLevelEstimates_220km_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")
CountryLevelData$Class<-factor(CountryLevelData$Class, levels = c("Amphibia", "Reptilia", "Mammalia", "Aves", "Vertebrates"))

# Convert CountryLevelData dataset into wide format:
UnknownSR_wide<-reshape2::dcast(CountryLevelData, ISO3 ~ Class, value.var="UnknownSR")
UnknownSR_wide_L95<-reshape2::dcast(CountryLevelData, ISO3 ~ Class, value.var="UnknownSR_L95")
UnknownSR_wide_U95<-reshape2::dcast(CountryLevelData, ISO3 ~ Class, value.var="UnknownSR_U95")
names(UnknownSR_wide)<-c("ISO3", "UnknownAmph", "UnknownRept", "UnknownMamm", "UnknownAves", "UnknownVert")
names(UnknownSR_wide_L95)<-c("ISO3", "UnknownSR_L95_Amph", "UnknownSR_L95_Rept", "UnknownSR_L95_Mamm", "UnknownSR_L95_Aves", "UnknownSR_L95_Vert")
names(UnknownSR_wide_U95)<-c("ISO3", "UnknownSR_U95_Amph", "UnknownSR_U95_Rept", "UnknownSR_U95_Mamm", "UnknownSR_U95_Aves", "UnknownSR_U95_Vert")
PropUnknown_wide<-reshape2::dcast(CountryLevelData, ISO3 ~ Class, value.var="PropUnknown")
names(PropUnknown_wide)<-c("ISO3", "PropUnknownAmph", "PropUnknownRept", "PropUnknownMamm", "PropUnknownAves", "PropUnknownVert")

# Bind all wide format datasets:
per_country_results<-cbind(UnknownSR_wide, UnknownSR_wide_L95[,-1], UnknownSR_wide_U95[,-1], PropUnknown_wide[,-1])
rm(UnknownSR_wide, UnknownSR_wide_L95, UnknownSR_wide_U95, PropUnknown_wide, CountryLevelData)

# Make the per country results relative to total (see standardization procedure described in the Methods):
TotalUnknownSR<-sum(per_country_results$UnknownVert, na.rm=T)
TotalUnknownSR_L95<-sum(per_country_results$UnknownSR_L95_Vert, na.rm=T)
TotalUnknownSR_U95<-sum(per_country_results$UnknownSR_U95_Vert, na.rm=T)
per_country_results$UnknownVert<-per_country_results$UnknownVert/TotalUnknownSR
per_country_results$UnknownSR_L95_Vert<-per_country_results$UnknownSR_L95_Vert/TotalUnknownSR_L95
per_country_results$UnknownSR_U95_Vert<-per_country_results$UnknownSR_U95_Vert/TotalUnknownSR_U95
range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
per_country_results$StdPropUnknownVert<-range01(per_country_results$PropUnknownVert)

# Combine country-level estimates with the attribute table of gadm_36 shapefile:
country_pol@data$n_order<-c(1:nrow(country_pol@data))
country_pol@data$GID_0<-as.character(country_pol@data$GID_0)
per_country_results$ISO3<-as.character(per_country_results$ISO3)
country_pol@data<-merge(country_pol@data, per_country_results, by.x="GID_0", by.y="ISO3", all.x=T)
country_pol@data<-country_pol@data[order(country_pol@data$n_order),]

# Prepare spatial objects for use in ggplot2:
country_pol@data$id<-rownames(country_pol@data)
country_pol_df<-fortify(country_pol, region="id")
country_pol_df<-join(country_pol_df, country_pol@data, by="id")

# Create the major plot:
MyPlot_2<-ggplot(data=country_pol_df) +
  geom_polygon(data=world_limit_df, aes(long, lat, group=group), colour="black", fill="lightcyan", size=0.25) + # plot world limits filled by the ocean
  geom_polygon(data=country_pol_df, aes(long, lat, group=group, fill=StdPropUnknownVert), colour="gray40", size=0.25) + # plot countries coloured by UnknownVert
  geom_polygon(data=world_limit_df, aes(long, lat, group=group), colour="black", fill=NA, size=0.25) + # replot world boundary box without fill colour
  coord_equal() +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.justification = c(0, 0),
        legend.position=c(0.05, 0.06), legend.direction="horizontal",
        legend.text=element_text(size=rel(0.4), hjust=0.5),
        legend.key.width=unit(0.4,"cm"), legend.key.height=unit(0.175,"cm"),
        legend.background=element_rect(fill="white", size=0.2, linetype='solid', colour='black'),
        legend.title=element_text(size=5),
        legend.title.align=0.5,
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background=element_blank(),
        plot.background=element_blank(),
        plot.margin = rep(unit(-0.03,"null"),4),
        panel.spacing = unit(0,"null")) +
  guides(fill=guide_colorbar(title="Country-wide prop.\nof unknown species\n(standardized)", title.position = "top", 
                             label=T, label.position="bottom", nbin=15, barwidth=5, 
                             draw.ulim=T, draw.llim=T, frame.colour="black", ticks=T)) +
  scale_fill_gradientn(colours=c('#ffffff', '#f0f0f0', '#d9d9d9', '#bdbdbd', '#969696', '#737373', '#525252', '#252525', '#000000'), # white to black
                       na.value = NA, values=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1), breaks=c(0.0, 0.5, 1)) +
  annotate(geom="rect", size=0.1, xmin=-17200000, xmax=-11300000,  ymin=-2500000, ymax=1000000, colour="black", fill="white")


# Get countries boundaries again and compute their area (in km²):
countries<-st_read("Shapefiles/gadm36_cea.shp",  stringsAsFactors=F) %>% dplyr::mutate(AREA_KM2 = as.numeric(st_area(.)) / 1e6)
countries<-merge(countries, per_country_results[,c(1:16)], by.x="GID_0", by.y="ISO3", all.x=T)

# Get the coordinates of centroids of countries for pie charts placement:
cen <- st_centroid(countries) 
coords <- st_coordinates(cen) 
rads <- sqrt(countries$AREA_KM2 + 1000000)*100 # calculate pie chart radius (for a piechart proportional to the country area, if needed)
dat <- countries %>% st_set_geometry(NULL) %>% 
  cbind(coords, rads) %>% # bind coordinates of centroids and rads size
  filter(AREA_KM2 > 1e4) # remove data from countries too small for visualization
rm(coords, rads, cen, countries)

# Simulate a new country to use as pie size legend:
dat[(nrow(dat)+1),]<-NA
row_to_track<-which(is.na(dat$GID_0)) # track the number of the new row
dat[row_to_track, (3:(ncol(dat)-3))] <- colSums(dat[,(3:(ncol(dat)-3))], na.rm=T) # use global values of UnknownSR for the simulated countries
dat<-rbind(dat, dat[row_to_track,], dat[row_to_track,]) # replicate the new row three times
row_to_track<-which(is.na(dat$GID_0)) # retrack the number of the new rows
dat[row_to_track, 1:2] <- c("EX1", "EX2", "EX3") # give names for the example countries (piecharts)

# Set the three piechart sizes for the legend:
dat$UnknownVert[row_to_track]<-c(0.012, 0.025, 0.075) 
dat$UnknownSR_L95_Vert[row_to_track]<-c(0.012, 0.025, 0.075)
dat$UnknownSR_U95_Vert[row_to_track]<-c(0.012, 0.025, 0.075)

# Set the coordinates of example piecharts:
dat$X[row_to_track]<-c(-16000000, -16000000, -14000000) 
dat$Y[row_to_track]<-c(-800000, -1800000, -1300000)
dat$rads[row_to_track]<-c(10000, 10000, 10000)

# Create a column with constant value to aid with overlaying empty piecharts (95% CI):
dat$NoValue<-0 # to plot empty piecharts

# Plot the pie charts on top of the countries centroids:
MyPieChart <- MyPlot_2 + ggnewscale::new_scale_fill()
MyPieChart <- MyPieChart + theme(legend.position = "none") +
  geom_scatterpie(data=dat, aes(x=X, y=Y, group=GID_0, r=(UnknownVert)*12000000), # 12000000 is a scale factor to expand size of piecharts
                  size=.025, sorted_by_radius=T, cols=colnames(dat[,c(4:7)])) +
  scale_fill_manual(name="", values=c('violetred2', 'slateblue2', 'orange', 'green3')) # amphibians, birds, mammals, reptiles

# Plot empty piecharts representing the 95% upper confidence interval:
MyPieChart2 <- MyPieChart + ggnewscale::new_scale_fill()
MyPieChart2 <- MyPieChart2 + geom_scatterpie(data=dat, aes(x=X, y=Y, group=GID_0, r=(UnknownSR_U95_Vert)*12000000), # 12000000 is a scale factor to expand size of piecharts
                                            linetype="dashed", size=.01, col="grey20", fill=NA, sorted_by_radius=T, cols=colnames(dat[,c(4,22)]))

# Plot empty piecharts representing the 95% lower confidence interval:
MyPieChart3 <- MyPieChart2 + ggnewscale::new_scale_fill()
MyPieChart3 <- MyPieChart3 + geom_scatterpie(data=dat, aes(x=X, y=Y, group=GID_0, r=(UnknownSR_L95_Vert)*12000000),# 12000000 is a scale factor to expand size of piecharts
                  fill=NA, col="grey20", linetype="dashed", size=.01, sorted_by_radius=T, cols=colnames(dat[,c(4,22)]))

# Plot the pie charts on top of the countries centroids:
ColourBarLegend_2<-ggpubr::as_ggplot((cowplot::get_legend(MyPlot_2))) # get only the legend
extent(world_limit) # get the size of the drawing canvas
MyPlot2_StdLegend <- MyPieChart3 + 
  annotation_custom(ggplotGrob(ColourBarLegend_2), xmin=-19000000, ymin=-8000000) +
  annotate("text", label=bquote('Country-wide % of\ntotal future discoveries'),
           x=-14250000, y=450000, color="black", fontface="bold", size=2, hjust=0.5)
ggsave("Figure4_Part2.png", plot=MyPlot2_StdLegend, width=6.7, height=3, units="in", bg = "transparent")


### Figure 4 - Part 3
# Get the global piechart to include in the legend:
global_piechart_df<-colSums(per_country_results[,(2:6)], na.rm=T)
global_piechart_df<-as.data.frame(t(global_piechart_df))
GlobalPieChart<-ggplot(global_piechart_df) +
  geom_scatterpie(data=global_piechart_df, aes(x=0, y=0, r=UnknownVert), size=.025, cols=colnames(global_piechart_df[,c(1:4)])) +
  scale_fill_manual(name="", values=c('violetred2', 'slateblue2', 'orange', 'green3')) + # amphibians, birds, mammals, reptiles
  coord_equal() +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.background=element_blank(),
        legend.position="none")


### Parts of Figure 4 were assemble externally in the InkScape software.
MyPlot1_StdLegend # assemblage-based predictions of future species discovery.
MyPlot2_StdLegend # country-based predictions of future species discovery.
GlobalPieChart # global piechart of future species discovery according to each tetrapod class.
ggsave("Figure4_Part1.pdf", plot=MyPlot1_StdLegend, width=6.7, height=3, units="in", bg = "transparent")
ggsave("Figure4_Part2.pdf", plot=MyPlot2_StdLegend, width=6.7, height=3, units="in", bg = "transparent")
ggsave("Figure4_Part3.pdf", plot=GlobalPieChart, width=2, height=2, units="in", bg = "transparent")
# Animal silhuoetes extracted from http://www.phylopic.org/: Dendrobates-azureus, Ardeosaurus-brevipes, Chrysocyon-brachyurus, Turdus-merula