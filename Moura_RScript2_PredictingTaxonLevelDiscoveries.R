############################################################################################################################
### Supporting Information to ###

# Title: Shortfalls and opportunities in terrestrial vertebrate species discovery
# Authors: Mario R. Moura 1,2,3; Walter Jetz1,2
# Journal: Nature Ecology and Evolution
# 1 Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, USA
# 2 Center for Biodiversity and Global Change, Yale University, New Haven, CT, USA
# 3 Department of Biological Sciences, Federal University of Paraíba, Areia, PB, Brazil
# * Corresponding author: mariormoura@gmail.com


# SUPPORTING SCRIPT 2: PREDICTING TAXON-LEVEL DISCOVERIES
############################################################################################################################
# SUPPORTING SCRIPT 2: PREDICTING TAXON-LEVEL DISCOVERIES

# Steps in this script:
#  1. Load, understand and prepare the dataset for analysis.
#  2. Check the taxonomic heirarchy and define higher-level grouping
#  3. Compute the discovery metrics for each higher-level grouping
#  4. Compute the discovery metrics for each taxonomic family

# First, clean workspace:
rm(list=ls()); gc()

# Install and load R packages needed to run the analysis:
needed_packages<-c("data.table", "dplyr")
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

# Load the species-level estimates needed to perform the taxon-level analysis:
trait_data<-fread("SpeciesLevelData.csv", stringsAsFactors=TRUE, encoding="UTF-8")
trait_data<-trait_data[, .(Binomial, Class, Order, Suborder, Family, Genus, DiscProb, DiscProb_L95, DiscProb_U95)]

# We have 9 columns in the remaining dataset, each of which is explained below:
# Binomial: Sspecies scientific name
# Class: Taxonomic class
# Order: Taxonomic order
# Suborder: Higher-level grouping
# Family: Taxonomic family
# YearOfDescription: Year in which the species was formally described for the first time
# DiscProb: discovery probability in the year 2015
# DiscProb_L95: lower bound of the discovery probability in the year 2015
# DiscProb_U95: upper bound of the discovery probability in the year 2015

# Remove species without data available:
trait_data<-trait_data[complete.cases(trait_data), ]

#####

# STEP 2 - CHECK THE TAXONOMIC HIERARCHY AND DEFINE HIGHER-LEVEL GROUPING
##########################################################################################################################
# STEP 2 - CHECK THE TAXONOMIC HIERARCHY AND DEFINE HIGHER-LEVEL GROUPING
rm(list=setdiff(ls(), c("trait_data")))
setwd("DefineYourDirectory")

# Verify if any taxonomic genus is present in more than one taxonomic family:
families_per_genera<-as.data.frame(trait_data %>%
                                     dplyr::group_by(Genus) %>%
                                     dplyr::summarise(NFamilies=length(unique(Family))))
families_per_genera<-families_per_genera[families_per_genera$NFamilies>=2,]
families_per_genera<-droplevels(families_per_genera); families_per_genera

# Verify if any taxonomic family is present in more than one taxonomic order:
order_per_family<-as.data.frame(trait_data %>%
                                  dplyr::group_by(Family) %>%
                                  dplyr::summarise(NOrders=length(unique(Order))))
order_per_family<-order_per_family[order_per_family$NOrders>=2,]
order_per_family<-droplevels(order_per_family); order_per_family
rm(order_per_family, families_per_genera)

# Higher-level grouping for amphibians.
# Some primitive anuran families:
non_neobatrachia<-c("Leiopelmatidae", "Ascaphidae", "Bombinatoridae", "Discoglossidae", "Alytidae", "Pipidae", 
                    "Rhinophrynidae", "Scaphiopodidae", "Pelodytidae", "Megophryidae","Pelobatidae")

# Within Neobratachia, get the Hyloidea clade (most frog species from the Neartic and Neotropic realms):
hyloidea_families<-c("Calyptocephalellidae", "Myobatrachidae", "Ceuthomantidae", "Brachycephalidae", "Eleutherodactylidae", "Craugastoridae",
                     "Hemiphractidae", "Hylidae", "Pelodryadidae", "Phyllomedusidae", "Bufonidae", "Dendrobatidae", "Aromobatidae", 
                     "Allophrynidae", "Centrolenidae", "Leptodactylidae", "Ceratophryidae", "Odontophrynidae", "Cycloramphidae", "Alsodidae",
                     "Hylodidae", "Telmatobiidae", "Batrachylidae", "Rhinodermatidae","Strabomantidae")

# Within Neobratachia, get the Ranoidea clade (most frog species from the Afrotropical, Paleartic, Indo-Malay, and 256 Australasia realms):
ranoidea_families<-c("Arthroleptidae", "Brevicipitidae", "Ceratobatrachidae", "Conrauidae", "Dicroglossidae", "Hemisotidae", "Hyperoliidae", 
                     "Mantellidae", "Micrixalidae", "Microhylidae",  "Nyctibatrachidae", "Odontobatrachidae", "Petropedetidae", 
                     "Phrynobatrachidae", "Ptychadenidae", "Pyxicephalidae", "Ranidae", "Ranixalidae", "Rhacophoridae")

# Other Neobratachia species not included in Ranoidea or Hyloidea (non-monophyletic set):
neobatrachia_others<-c("Nasikabatrachidae", "Sooglossidae", "Calyptocephalellidae", "Myobatrachidae", "Heleophrynidae", "Limnodynastidae")

# Higher-level grouping for reptiles.
# Get the taxonomic families representing the Cryptodira (hidden-necked turtles) and Pleurodira (side-necked turtles) clades: 
cryptodira_families<-c("Chelydridae", "Emydidae", "Geoemydidae", "Platysternidae", "Carettochelyidae", "Dermatemydidae",
                       "Kinosternidae", "Testudinidae", "Trionychidae")
pleurodira_families<-c("Chelidae", "Pelomedusidae", "Podocnemididae")

# Get the taxonomic families representing the Gekkota clade (gekkos and relatives): 
gekkota_families<-c("Gekkonidae", "Carphodactylidae", "Diplodactylidae", "Phyllodactylidae", "Sphaerodactylidae", "Pygopodidae", "Eublepharidae")

# Iguania clade (iguanas, chameleons, and relatives): 
iguania_families<-c("Chamaeleonidae", "Iguanidae", "Agamidae", "Dactyloidae", "Opluridae", "Leiocephalidae", "Liolaemidae", "Phrynosomatidae", "Tropiduridae",
                    "Leiosauridae", "Polychrotidae", "Corytophanidae", "Crotaphytidae", "Hoplocercidae")

# Scincoidea clade (skins and relatives):
scincoidea_families<-c("Gerrhosauridae", "Cordylidae", "Scincidae", "SCINCIDAE", "Xantusiidae")

# Lacertoidea clade (teiids, lacertids, amphisbaenians, and relatives):
lacertiformes_families<-c("Teiidae",  "Gymnophthalmidae", "Lacertidae", "Amphisbaenidae", "Cadeidae",
                          "Trogonophiidae", "Bipedidae", "Blanidae", "Rhineuridae")

# Anguimorpha clade (glass lizards, monitors, and relatives):
anguimorpha_families<-c( "Anguidae", "Diploglossidae", "Anniellidae", "Xenosauridae", "Helodermatidae", "Shinisauridae", "Lanthanotidae", "Varanidae")

# Dibamids of blink skins (sometimes referred as Dibamoidea clade):
dibamia_family<-"Dibamidae"

# Serpentes clade (snakes):
serpentes_families<-c("Typhlopidae", "Colubridae", "Leptotyphlopidae", "Natricidae", "Gerrhopilidae", "Anomalepididae", "Dipsadidae", "Homalopsidae",
                      "Lamprophiidae", "Xenodermatidae", "Pseudoxenodontidae", "Viperidae", "Uropeltidae", "Tropidophiidae", "Elapidae", "Anomochilidae",
                      "Xenophidiidae", "Xenotyphlopidae", "Cylindrophiidae", "Pareatidae", "Boidae", "Xenopeltidae", "Pythonidae", "Acrochordidae",
                      "Aniliidae", "Bolyeridae", "Loxocemidae")

# Higher-level grouping for birds.
# Split Passeriformes into Oscines (songbirds) and Suboscines (passerine birds):
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

# Set a new higher-level grouping herein named as "Suborder".
# At this point, it behaves like Order, but some more informative divisions will be mannually added:
# In amphibians, Caudata and Gymonphiona were kept as they are.
# In reptiles, Crocodilia was kept as it is.
# In mammals, all orders were kept as they are.
# In birds, only Passeriformes was splitted between Oscines and Suboscines.
trait_data$Suborder<-trait_data$Order
trackthis<-which(names(trait_data)=="Suborder")

# Replace suborder for amphibian species:
trait_data[which(trait_data$Family %in% as.character(ranoidea_families)),trackthis]<-"Neob. (Ranoidea)"
trait_data[which(trait_data$Family %in% as.character(hyloidea_families)),trackthis]<-"Neob. (Hyloidea)"
trait_data[which(trait_data$Family %in% as.character(neobatrachia_others)),trackthis]<-"Neobatrachia (others)"
trait_data[which(trait_data$Family %in% as.character(non_neobatrachia)),trackthis]<-"Non-Neobatrachia"

# Replace suborder for reptile species:
trait_data[which(trait_data$Family %in% as.character(pleurodira_families)),trackthis]<-"Pleurodira"
trait_data[which(trait_data$Family %in% as.character(cryptodira_families)),trackthis]<-"Cryptodira"
trait_data[which(trait_data$Family %in% as.character(gekkota_families)),trackthis]<-"Gekkota"
trait_data[which(trait_data$Family %in% as.character(iguania_families)),trackthis]<-"Iguania"
trait_data[which(trait_data$Family %in% as.character(scincoidea_families)),trackthis]<-"Scincoidea"
trait_data[which(trait_data$Family %in% as.character(lacertiformes_families)),trackthis]<-"Lacertiformes"
trait_data[which(trait_data$Family %in% as.character(anguimorpha_families)),trackthis]<-"Anguimorpha"
trait_data[which(trait_data$Family %in% as.character(dibamia_family)),trackthis]<-"Dibamoidea"
trait_data[which(trait_data$Family %in% as.character(serpentes_families)),trackthis]<-"Serpentes"

# Replace suborder for bird species:
trait_data[which(trait_data$Family %in% as.character(oscines_passeri)),trackthis]<-"Oscines"
trait_data[which(trait_data$Family %in% as.character(suboscines_tyranni)),trackthis]<-"Suboscines"

# The Passeriformes species from the Acanthisittidae family are neither Oscines nor Suboscines.
# We used the taxonomic family to define their suborder (they are only two species):
trait_data[which(trait_data$Family %in% as.character("Acanthisittidae")),1] # which species
trait_data[which(trait_data$Family %in% as.character("Acanthisittidae")),trackthis]<-"Acanthisittidae"

#####

# STEP 3 - COMPUTE THE DISCOVERY METRICS FOR EACH HIGHER-LEVEL GROUPING
##########################################################################################################################
# STEP 3 - COMPUTE THE DISCOVERY METRICS FOR EACH HIGHER-LEVEL GROUPING
rm(list=setdiff(ls(), c("trait_data")))
setwd("DefineYourDirectory")

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

# Reorder levels of taxonomic class (just for organization purposes)
trait_data$Class<-factor(trait_data$Class,
                           levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"),
                           labels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Compute the estimated 'PropUnknown' for each level of Suborder (higher-level grouping):
vert_summary<-trait_data %>%
  dplyr::group_by(Class, Order, Suborder) %>% 
  dplyr::summarise(Richness=length(Suborder),
                   PropUnknown= (1 - geometric.mean(DiscProb, na.rm=T)),
                   PropUnknown_L95= (1 - geometric.mean(DiscProb_U95, na.rm=T)), # note that the upper and lower bounds are now inverted
                   PropUnknown_U95= (1 - geometric.mean(DiscProb_L95, na.rm=T))) # note that the upper and lower bounds are now inverted

# Compute the estimated UnknownSR for each level of Suborder (higher-level grouping):
vert_summary$UnknownSR<- (vert_summary$Richness/(1-vert_summary$PropUnknown)) - vert_summary$Richness
vert_summary$UnknownSR_L95<- (vert_summary$Richness/(1-vert_summary$PropUnknown_L95)) - vert_summary$Richness
vert_summary$UnknownSR_U95<- (vert_summary$Richness/(1-vert_summary$PropUnknown_U95)) - vert_summary$Richness

# Filter the taxon with at least 5 species and export the results:
vert_summary<-vert_summary[vert_summary$Richness>=5,]
dir.create("TaxonLevelAnalysis", showWarnings=F)
fwrite(vert_summary, "TaxonLevelAnalysis/TaxonLevelEstimates_HigherLevelGrouping_WithoutStandardization.csv")

# The sensitivity analysis revealed that 'UnknownSR' and 'PropUnknown' are underestimated.
## To provide more informative metrics of discovery potential, we divided UnknownSR by the total number of estimated discoveries (i.e., 
## sum of UnknownSR across taxa) to provide the estimated percent of total discoveries.
## PropUnknown was standardized to vary between 0 and 1, by first subtracting the minimum observed for each vertebrate class and then
## dividing by the respective range of PropUnknown. The value of 1 indicated the taxon with the highest proportion of unknown species 
## (whatever such number might be), and not necessarily a taxon with 100% of unknown species.

# Create a empty data frame to receive the standardized values of "PropUnknown" and 'UnknownSR':
Std_vert_summary<-as.data.frame(matrix(nrow=1, ncol=ncol(vert_summary), NA))
names(Std_vert_summary)<-names(vert_summary)

# Perform the standardization for each vertebrate class:
for(i in 1:nlevels(vert_summary$Class)){ # i vertebrate groups
  
  vert_summary<-fread("TaxonLevelAnalysis/TaxonLevelEstimates_HigherLevelGrouping_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  vert_filtered<-vert_summary %>% filter(Class==levels(vert_summary$Class)[i])
  vert_filtered<-droplevels(vert_filtered)
  
  # Make the 'UnknownSR' relative to the total:
  track_denominator<-sum(vert_filtered$UnknownSR, na.rm=T)
  vert_filtered$UnknownSR<-vert_filtered$UnknownSR/track_denominator # keep the same denominator
  vert_filtered$UnknownSR_L95<-vert_filtered$UnknownSR_L95/track_denominator # keep the same denominator
  vert_filtered$UnknownSR_U95<-vert_filtered$UnknownSR_U95/track_denominator # keep the same denominator
  
  # Make the 'PropUnknown' to vary between 0 and 1:
  track_denominator<-(max(vert_filtered$PropUnknown, na.rm=T) - min(vert_filtered$PropUnknown, na.rm=T))
  track_min<-min(vert_filtered$PropUnknown, na.rm=T)
  vert_filtered$PropUnknown<-(vert_filtered$PropUnknown - track_min) / track_denominator
  vert_filtered$PropUnknown_L95<-(vert_filtered$PropUnknown_L95 - track_min) / track_denominator
  vert_filtered$PropUnknown_U95<-(vert_filtered$PropUnknown_U95 - track_min) / track_denominator
  
  # Bind
  Std_vert_summary<-rbind(Std_vert_summary, vert_filtered)
}
Std_vert_summary<-Std_vert_summary[-1,] # remove the first  (empty) row

# Export the standardized values of 'PropKnown' and 'UnknownSR' per higher-level grouping:
fwrite(Std_vert_summary, "TaxonLevelAnalysis/TaxonLevelEstimates_HigherLevelGrouping.csv")

#####

# STEP 4 - COMPUTE THE DISCOVERY METRICS FOR EACH TAXONOMIC FAMILY
##########################################################################################################################
# STEP 4 - COMPUTE THE DISCOVERY METRICS FOR EACH TAXONOMIC FAMILY
rm(list=setdiff(ls(), c("trait_data")))
setwd("DefineYourDirectory")

# Same as in Step3, but at the level of taxonomic family.
# Create a function to compute the geometric mean for a set of species:
geometric.mean <- function(x,na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm)) }

# Reorder levels of taxonomic class (just for organization purposes)
trait_data$Class<-factor(trait_data$Class,
                         levels = c("Amphibia", "Reptilia", "Mammalia", "Aves"),
                         labels = c("Amphibia", "Reptilia", "Mammalia", "Aves"))

# Compute the estimated 'PropUnknown' for each level of Suborder (higher-level grouping):
vert_summary<-trait_data %>%
  dplyr::group_by(Class, Order, Family) %>% 
  dplyr::summarise(Richness=length(Family),
                   PropUnknown= (1 - geometric.mean(DiscProb, na.rm=T)),
                   PropUnknown_L95= (1 - geometric.mean(DiscProb_U95, na.rm=T)), # note that the upper and lower bounds are now inverted
                   PropUnknown_U95= (1 - geometric.mean(DiscProb_L95, na.rm=T))) # note that the upper and lower bounds are now inverted

# Compute the estimated UnknownSR for each level of Suborder (higher-level grouping):
vert_summary$UnknownSR<- (vert_summary$Richness/(1-vert_summary$PropUnknown)) - vert_summary$Richness
vert_summary$UnknownSR_L95<- (vert_summary$Richness/(1-vert_summary$PropUnknown_L95)) - vert_summary$Richness
vert_summary$UnknownSR_U95<- (vert_summary$Richness/(1-vert_summary$PropUnknown_U95)) - vert_summary$Richness

# Filter the taxon with at least 5 species and export the results:
vert_summary<-vert_summary[vert_summary$Richness>=5,]
dir.create("TaxonLevelAnalysis", showWarnings=F)
fwrite(vert_summary, "TaxonLevelAnalysis/TaxonLevelEstimates_TaxonomicFamily_WithoutStandardization.csv")

# Create a empty data frame to receive the standardized values of 'PropUnknown' and 'UnknownSR':
Std_vert_summary<-as.data.frame(matrix(nrow=1, ncol=ncol(vert_summary), NA))
names(Std_vert_summary)<-names(vert_summary)

# Perform the standardization for each vertebrate class:
for(i in 1:nlevels(vert_summary$Class)){ # i vertebrate groups
  
  vert_summary<-fread("TaxonLevelAnalysis/TaxonLevelEstimates_TaxonomicFamily_WithoutStandardization.csv", stringsAsFactors=TRUE, encoding="UTF-8")
  vert_filtered<-vert_summary %>% filter(Class==levels(vert_summary$Class)[i])
  vert_filtered<-droplevels(vert_filtered)
  
  # Make the 'UnknownSR' relative to the total:
  track_denominator<-sum(vert_filtered$UnknownSR, na.rm=T)
  vert_filtered$UnknownSR<-vert_filtered$UnknownSR/track_denominator # keep the same denominator
  vert_filtered$UnknownSR_L95<-vert_filtered$UnknownSR_L95/track_denominator # keep the same denominator
  vert_filtered$UnknownSR_U95<-vert_filtered$UnknownSR_U95/track_denominator # keep the same denominator
  
  # Make the 'PropUnknown' to vary between 0 and 1:
  track_denominator<-(max(vert_filtered$PropUnknown, na.rm=T) - min(vert_filtered$PropUnknown, na.rm=T))
  track_min<-min(vert_filtered$PropUnknown, na.rm=T)
  vert_filtered$PropUnknown<-(vert_filtered$PropUnknown - track_min) / track_denominator
  vert_filtered$PropUnknown_L95<-(vert_filtered$PropUnknown_L95 - track_min) / track_denominator
  vert_filtered$PropUnknown_U95<-(vert_filtered$PropUnknown_U95 - track_min) / track_denominator
  
  # Bind
  Std_vert_summary<-rbind(Std_vert_summary, vert_filtered)
}
Std_vert_summary<-Std_vert_summary[-1,] # remove the first  (empty) row

# Export the standardized values of 'PropKnown' and 'UnknownSR' per taxonomic family:
fwrite(Std_vert_summary, "TaxonLevelAnalysis/TaxonLevelEstimates_TaxonomicFamily.csv")

#####