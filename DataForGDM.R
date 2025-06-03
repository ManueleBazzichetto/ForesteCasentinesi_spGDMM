#prepare data to be used for fitting spGDMM

#fit splines
library(splines)
library(splines2)

#for rdist function
library(fields)

#compute Bray-Curtis
library(vegan)
#compute geographical distance
library(geosphere)

#viz
library(ggplot2)
library(ggpubr)

#handle spatial data - not necessarily needed
library(sf)
library(mapview)


#create a copy of EVA_forcas_meta and EVA_forcas_species
#the copies will be modified due to filtering of data for fitting spGDMMs

exists('Forcas_env'); exists('Forcas_sp') #FALSE*2

#environmental data
Forcas_env <- EVA_forcas_meta

#species composition data
Forcas_sp <- EVA_forcas_species

#------------------------------------sanity checks on input data

#check number of unique plots
length(unique(Forcas_sp$PlotObservationID)) #943

#check class of 'PlotObservationID'
class(Forcas_env$PlotObservationID); class(Forcas_sp$PlotObservationID) #integer*2

#coerce PlotObservationID to chr
Forcas_env$PlotObservationID <- as.character(Forcas_env$PlotObservationID)
Forcas_sp$PlotObservationID <- as.character(Forcas_sp$PlotObservationID)

#filter out observations from species data not matching Plot ID in meta data [no problem here]
sum(!unique(Forcas_sp$PlotObservationID) %in% Forcas_env$PlotObservationID) #0
sum(!Forcas_env$PlotObservationID %in% unique(Forcas_sp$PlotObservationID)) #0
sum(is.na(match(unique(Forcas_sp$PlotObservationID), Forcas_env$PlotObservationID))) #0
sum(is.na(match(Forcas_env$PlotObservationID, unique(Forcas_sp$PlotObservationID)))) #0

#------------------------------------check on plot duplicates (species composition and coordinates)

#find and drop plots sampled exactly at the same location, having exactly the same species composition and abundances, 
#but being associated with different sampling years

find_duplicated_plots <- lapply(unique(Forcas_sp$PlotObservationID), function(id) {
  dtf_plot <- Forcas_sp[Forcas_sp$PlotObservationID == id, c('Matched concept corrected', 'Cover %')]
  data_plot <- setNames(object = dtf_plot[['Cover %']], nm = dtf_plot[['Matched concept corrected']])
  #loop over other plot IDs to compare species list and cover
  other_ids <- unique(Forcas_sp$PlotObservationID)
  other_ids <- other_ids[!other_ids %in% id]
  all_compare_res <- sapply(other_ids, function(ot_id) {
    dtf_oth <- Forcas_sp[Forcas_sp$PlotObservationID == ot_id, c('Matched concept corrected', 'Cover %')]
    data_oth <- setNames(object = dtf_oth[['Cover %']], nm = dtf_oth[['Matched concept corrected']])
    if((isTRUE(length(setdiff(names(data_plot), names(data_oth))) == 0) && isTRUE(length(setdiff(names(data_oth), names(data_plot))) == 0))) {
      data_oth <- data_oth[names(data_plot)]
      comp_res <- all.equal(data_oth, data_plot)
      if(isTRUE(comp_res)) comp_res <- TRUE else comp_res <- FALSE
      return(comp_res)
    } else {
      return(FALSE)
    }
  })
  doppelganger <- other_ids[which(all_compare_res)]
})

names(find_duplicated_plots) <- unique(Forcas_sp$PlotObservationID)

#Filter(function(x) length(x) > 0, find_duplicated_plots)
 
duply_plots <- unlist(find_duplicated_plots)

#check duplicates - notice that duplicate pairs are themselves duplicated (dupl A vs. B == dupl B vs. a) - so half of duply_plots can be dropped
duply_plots <- duply_plots[c(1, 3, 5)]

#case 1) same year
cbind(Forcas_sp[Forcas_sp$PlotObservationID == names(duply_plots)[1], c('Matched concept corrected', 'Cover %')], 
      Forcas_sp[Forcas_sp$PlotObservationID == duply_plots[1], c('Matched concept corrected', 'Cover %')])
cbind(Forcas_env[Forcas_env$PlotObservationID == names(duply_plots)[1], c('Year')], 
      Forcas_env[Forcas_env$PlotObservationID == duply_plots[1], c('Year')])

#case 2) same year
cbind(Forcas_sp[Forcas_sp$PlotObservationID == names(duply_plots)[2], c('Matched concept corrected', 'Cover %')], 
      Forcas_sp[Forcas_sp$PlotObservationID == duply_plots[2], c('Matched concept corrected', 'Cover %')])
cbind(Forcas_env[Forcas_env$PlotObservationID == names(duply_plots)[2], c('Year')], 
      Forcas_env[Forcas_env$PlotObservationID == duply_plots[2], c('Year')])

#case 3) same year
cbind(Forcas_sp[Forcas_sp$PlotObservationID == names(duply_plots)[3], c('Matched concept corrected', 'Cover %')], 
      Forcas_sp[Forcas_sp$PlotObservationID == duply_plots[3], c('Matched concept corrected', 'Cover %')])
cbind(Forcas_env[Forcas_env$PlotObservationID == names(duply_plots)[3], c('Year')], 
      Forcas_env[Forcas_env$PlotObservationID == duply_plots[3], c('Year')])

#drop duplicates from the datasets

#keep names(duply_plots) as PlotID
Forcas_env <- Forcas_env[!Forcas_env$PlotObservationID %in% duply_plots, ]

Forcas_sp <- Forcas_sp[Forcas_sp$PlotObservationID %in% Forcas_env$PlotObservationID, ]

#check
names(duply_plots) %in% Forcas_env$PlotObservationID
duply_plots %in% Forcas_env$PlotObservationID

names(duply_plots) %in% Forcas_sp$PlotObservationID
duply_plots %in% Forcas_sp$PlotObservationID

#check if some plots have exactly the same coordinates as this could create a problem with the kernel function
#basically, the 'kernel matrix' is going to be singular
#notice that ReSurvey plots likely have duplicated coordinates

length(which(duplicated(Forcas_env[c('lon', 'lat')]))) #172

#unique duplicated locations
duply_loc <- unique(Forcas_env[which(duplicated(Forcas_env[c('lon', 'lat')])), c('lon', 'lat')])

#coerce to list
duply_loc <- as.list(duply_loc)

#find PlotObservationID with duplicated coordinates - some coordinates are associated with more than 2 locations: 25, 26, 27, 28, 29, 30
duply_loc_plotid <- Map(function(x, y) {
  Forcas_env[Forcas_env$lon == x & Forcas_env$lat == y, 'PlotObservationID']
}, x = duply_loc$lon, y = duply_loc$lat)

lengths(duply_loc_plotid)

#1 m2 plots will be dropped out from the analyses at a later stage
Forcas_env[Forcas_env$PlotObservationID %in% duply_loc_plotid[[25]], ] #plot size 1 - Fangacci
Forcas_env[Forcas_env$PlotObservationID %in% duply_loc_plotid[[26]], ] #plot size 1 - Fangacci
Forcas_env[Forcas_env$PlotObservationID %in% duply_loc_plotid[[27]], ] #plot size 1 - Fangacci
Forcas_env[Forcas_env$PlotObservationID %in% duply_loc_plotid[[28]], ] #plot size 1 - Fonte Sodo Conti
Forcas_env[Forcas_env$PlotObservationID %in% duply_loc_plotid[[29]], ] #plot size 1 - Casone della Burraia
Forcas_env[Forcas_env$PlotObservationID %in% duply_loc_plotid[[30]], ] #plot size 300 - Sasso Fratino

#exclude elements of duply_loc_plotid from 25 to 29 as these will be dropped out when filtering for plot size
duply_loc_plotid <- duply_loc_plotid[-c(25, 26, 27, 28, 29)]

#sample 1 random plot from duplicates - i is always going to have length > 1
#I should end up with a vector with length equal to length(duply_loc_plotid) + 18
set.seed(74859)

duply_loc_drop <- unlist(lapply(duply_loc_plotid, function(i) {
  if(length(i) == 2) {
    res <- sample(x = i, size = 1, replace = F)
    return(res)
    } else {
      res <- sample(x = i, size = (length(i) - 1), replace = F) 
      return(res)
      }
  }))

#drop plots
Forcas_env <- Forcas_env[!Forcas_env$PlotObservationID %in% duply_loc_drop, ]

#------------------------------------check how many habitat types are included in Forcas_env

#check frequency of habitat types
with(Forcas_env, table(Expert_system))

#assignment of habitat category follows Excel file in: https://1drv.ms/f/s!AhheMx-31qY5jqIVbdGFkojTLiCkuA?e=vYZBZd
#name of the file: EUNIS-ESy-2025-01-18_legendfromESYfile_excel

Hab_types <- unique(Forcas_env$Expert_system)

Hab_types <- c("~" = "unclassified", "T" = "Forest", "T19" = "Temp.sub.thermoph.decid.forest", "T1F" = "Ravine.forest", "T3M" = "Coniferous.plant.non.site.nat.trees",
  "R1A" = "Semi.dry.peren.calcareous.grassland", "T13" = "Temp.hardwood.riparian.forest", "T1E" = "Carpinus.Quercus.mesic.decid.forest",
  "R" = "Grasslands", "T31" = "Temp.mountain.Picea.forest", "V11" = "Intensive.unmixed.crops", "R16" = "Peren.rocky.grassland.Central.Europe",
  "R55" = "Lowland.moist.wet.tallherb.fringe", "Sa" = "Scrub", "R14" = "Peren.rocky.grassland.Italian.Peninsula",
  "R22" = "Low.medium.altit.hay.meadow", "R36" = "Moist.wet.mesotroph.to.eutrop.pasture", "S53" = "Spartium.junceum.scrub",
  "R1F" = "Medit.annual.rich.dry.grassland", "T32" = "Temp.mountain.Abies.forest", "S35" = "Temp.submedit.thorn.scrub",
  "U37" = "Temperate.lowland.to.mont.base.rich.inland.cliff", "S31,S35" = "Mix.of.scrub.31_35", "U27" = "Temp.lowland.to.mont.base.rich.scree",
  "T33" = "Medit.mount.Abies.forest", "S22" = "Alpine.subalpine.ericoid.heath", "S33" = "Lowland.to.mont.temp.submed.genistoid.scrub",
  "V" = "Man.made.habitats", "R44" = "Arctic.alpine.calc.grassland", "R1M" = "Lowland.to.mont.dry.mesic.grassland.Nardus.str.",
  "T11" = "Temp.Salix.Populus.riparian.forest", "U38" = "Medit.base.rich.inland.cliff", "S35,S53" = "Mix.of.scrub.35_53",
  "S31" = "Lowland.to.mont.temp.submed.Juniperus.scrub", "R57" = "Herbaceous.forest.clearing.veget",
  "U23" = "Temp.lowland.to.mont.siliceous.scree", "R52" = "Forest.fringe.acidic.nutrient.poor.soils",
  "R56" = " Mont.to.subalp.moist.wet.tall.herb.fern.fringe", "Qb" = "Wetlands", "Q52" = "Small.helophyte.bed",
  "R35" = "Moist.wet.mesotroph.to.eutroph.hay.meadow", "Pa" = "Brackish.water.veget.", "T17" = "Fagus.forest.on.non.acid.soils") 

Hab_types.df <- data.frame(Hab_nm = unname(Hab_types),
                           Hab_code = names(Hab_types),
                           Hab_simpl = substr(names(Hab_types), start = 1, stop = 1))

#create key to go from Hab_code to simplified classification of habitats - higher hierarchical level
Hab_key <- unique(Hab_types.df$Hab_simpl)

Hab_key <- setNames(object = c('Uncl.', 'Forest', 'Grassland', 'Man.made.hab', 'Scrub', 'Inland.sparse.veg', 'Wetlands', 'Surface.waters'),
                    nm = Hab_key)

#add forest vs. non-forest column
Hab_types.df$For_NonFor <- ifelse(Hab_types.df$Hab_simpl == "T", 'Forest', 'NonForest')


Hab_types.df$Hab_simpl <- unname(Hab_key[Hab_types.df$Hab_simpl])

#join Hab_simpl and For_NonFor to Forcas_env: these columns will be used below to filter out habitats not used in the analyses

all(Hab_types.df$Hab_code %in% Forcas_env$Expert_system) #T
all(unique(Forcas_env$Expert_system) %in% Hab_types.df$Hab_code) #T

colnames(Hab_types.df)[c(1, 3, 4)] %in% colnames(Forcas_env) #FALSE*3

Forcas_env <- dplyr::left_join(x = Forcas_env, y = Hab_types.df, by = c('Expert_system' = 'Hab_code'))

#------------------------------------drop observations with missing plot size or coordinates

#drop observations with missing coordinates
sum(!complete.cases(Forcas_env[c('lon', 'lat')])) #0 missing
anyNA(Forcas_env[c('lon', 'lat')]) #FALSE

#drop observations with missing plot size
sum(is.na(Forcas_env$Plot_size)) #208 (!)

Forcas_env <- Forcas_env[!is.na(Forcas_env$Plot_size), ]

#------------------------------------drop habitat types to be excluded from the analyses

#check how many observations are associated with habitats to be excluded
sum(Forcas_env$Hab_simpl %in% c('Uncl.', 'Man.made.hab', 'Wetlands', 'Surface.waters')) #190

Forcas_env <- Forcas_env[!Forcas_env$Hab_simpl %in% c('Uncl.', 'Man.made.hab', 'Wetlands', 'Surface.waters'), ]

#------------------------------------remove plots with size < than 10 m2 for non-forests and 100 m2 for forests

#this stage should drop put plots with duplicated coordinates - check for duplicated coordinates after filtering out obs with a too small plot size

#check range of plot size
range(Forcas_env$Plot_size)
hist(Forcas_env$Plot_size)
table(Forcas_env$Plot_size)

Obs_to_keep <- ifelse((Forcas_env$For_NonFor == "Forest" & Forcas_env$Plot_size >= 100) | (Forcas_env$For_NonFor == "NonForest" & Forcas_env$Plot_size >= 10), TRUE, FALSE)

#observations to keep
sum(Obs_to_keep) #374

#View(Forcas_env[!Obs_to_keep, c(5, 48)])

#drop observations
Forcas_env <- Forcas_env[Obs_to_keep, ]

#check duplicated coordinates
sum(duplicated(Forcas_env[c('lon', 'lat')])) #0

#------------------------------------remove plots that were resurveyed

#check how many resurveyed plots are included in the dataset
unique(Forcas_env$`ReSurvey plot (Y/N)`) #"N" / "Y" / ""

sum(is.na(Forcas_env$`ReSurvey plot (Y/N)`)) #0
sum(Forcas_env$`ReSurvey plot (Y/N)` == "Y") #75

#drop resurveyed plots
Forcas_env <- Forcas_env[Forcas_env$`ReSurvey plot (Y/N)` != "Y", ]

#------------------------------------match species data

Forcas_sp <- Forcas_sp[Forcas_sp$PlotObservationID %in% Forcas_env$PlotObservationID, ]

#fix issue of some plots having multiple times the same species with different cover

#function for combining different cover for the same species at the same layer
combine.cover <- function(x) {
  while (length(x) > 1) {
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
    }
  return(x)
  }


#find plots with multiple entries for the same species
duply_species <- sapply(unique(Forcas_sp$PlotObservationID), function(id) {
  any(duplicated(Forcas_sp[Forcas_sp$PlotObservationID == id, 'Matched concept corrected']))
}) 

#how many cases? 
sum(duply_species) #11

duply_sp_id <- names(which(duply_species)) 

duply_species_in_plot <- lapply(duply_sp_id, function(id) {
  plot_with_dup <- Forcas_sp[Forcas_sp$PlotObservationID == id, 'Matched concept corrected']
  dup_sp <- plot_with_dup[which(duplicated(plot_with_dup))]
  return(dup_sp)
})

names(duply_species_in_plot) <- duply_sp_id

#check cases - same species for different layer
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[1] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[1]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[2] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[2]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[3] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[3]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[4] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[4]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[5] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[5]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[6] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[6]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[7] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[7]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[8] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[8]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[9] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[9]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[10] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[10]], ]
Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[11] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[11]], ]

#example
duply_cov.ex <- Forcas_sp[Forcas_sp$PlotObservationID == duply_sp_id[1] & Forcas_sp$`Matched concept corrected` %in% duply_species_in_plot[[1]], 'Cover %']

combine.cover(duply_cov.ex) #88.36

#combine cover - the for loop should drop sum(lengths(duply_species_in_plot)) observations (22)
for(i in names(duply_species_in_plot)) {
  sp_to_comb <- duply_species_in_plot[[i]]
  for(j in sp_to_comb) {
    sp_position <- which(Forcas_sp$PlotObservationID == i & Forcas_sp$`Matched concept corrected` == j)
    cov_to_comb <- Forcas_sp[sp_position, 'Cover %']
    Forcas_sp[sp_position[1], 'Cover %'] <- combine.cover(cov_to_comb)
    #drop duplicates
    Forcas_sp <- Forcas_sp[-sp_position[-1], ]
  }
  }

rm(sp_to_comb, sp_position, cov_to_comb)
rm(i, j)

#check
sum(sapply(unique(Forcas_sp$PlotObservationID), function(id) {
  any(duplicated(Forcas_sp[Forcas_sp$PlotObservationID == id, 'Matched concept corrected']))
  })) #0 right!


#------------------------------------subset and re-format Forcas_sp from long to wide

#select columns needed for the analyses
Forcas_sp <- Forcas_sp[c('PlotObservationID', 'Matched concept corrected', 'Cover %')]

#rename columns
colnames(Forcas_sp) <- c('ID', 'Species_name', 'Cover_perc') 

#check range of cover values
range(Forcas_sp$Cover_perc) #0.1 - 98.00 (no NA)
anyNA(Forcas_sp$Cover_perc) #F
length(unique(Forcas_sp$Species_name)) #561

#from wide to long format
fcas_spec_mat <- as.data.frame(tidyr::pivot_wider(data = Forcas_sp, names_from = 'Species_name', values_from = 'Cover_perc', values_fill = 0))

#assign PlotID to row.names
row.names(fcas_spec_mat) <- fcas_spec_mat$ID

#drop ID column
fcas_spec_mat$ID <- NULL


#------------------------------------create locality and environmental matrices

fcas_loc_mat <- Forcas_env[c('PlotObservationID', 'lon', 'lat')]

colnames(fcas_loc_mat)[1] <- 'ID'

fcas_env_use <- Forcas_env[c('tmean_avg', 'prcp_avg', 'slope', 'east', 'TCW')]

cor(fcas_env_use)


#check if PlotID order is the same between species, location and env mat
identical(row.names(fcas_spec_mat), fcas_loc_mat$ID) #T


#------------------------------------generate input data for spGDMM

#------------------Dissimilarity in species composition - Z

#save number of sites
N_sites <- nrow(fcas_loc_mat)

#compute Bray-Curtis distance and check proportions of 1s and 0s
bc_distmat <- as.matrix(vegan::vegdist(fcas_spec_mat, 'bray'))

#get vector of observed dissimilarity
Obs_Z <- bc_distmat[upper.tri(bc_distmat)] #length equal to (N_sites*(N_sites - 1))/2

Smp_size <- length(Obs_Z)

#range
range(Obs_Z) #0.015 - 1 [no 0 dissimilarity]

#check perfect and non-perfect dissimilarity
Perf_Z <- which(Obs_Z == 1)

NonPerf_Z <- which(Obs_Z != 1)
N_nonperfZ <- length(NonPerf_Z) #37085

#how many 1s
N_perfZ <- length(Perf_Z) #7466
#prop of 1s
N_perfZ/Smp_size #17%
mean(Obs_Z == 1) #17%

#check
(N_nonperfZ + N_perfZ) == Smp_size #T


#------------------Geographical distance matrix

#compute distance matrix - distance is computed assuming spherical Earth and scaled to km
geo_distmat <- distm(cbind(fcas_loc_mat$lon, fcas_loc_mat$lat))/1000

#get vector of geographical distances
Site_geodist <- geo_distmat[upper.tri(geo_distmat)]

#check range
range(Site_geodist) #0.03 - 41.6 km

#------------------Exploratory analysis: how does Z changes along geo distance

Z_geodist.df <- data.frame(Z = Obs_Z, GeoDist = Site_geodist)

Z_along_geodist <- ggplot(Z_geodist.df, aes(x = GeoDist, y = Z)) +
  geom_hex() +
  geom_smooth(se = F, col = 'red', lwd = 2) +
  scale_fill_viridis_c() +
  ylab('Observed Beta-dissmilarity') + xlab('Geographical distance') + 
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 14))

#change in mean and variance of Z given distance

#bin distances in categories - 5km step
Z_geodist.df$GeoDist_cat <- with(Z_geodist.df, cut(GeoDist, breaks = c(seq(0, 35, by = 5), 42), labels = F)) 

#compute mean and variance of log-Beta dissimilarity (scale of latent Z) - these are conditional mean and variance | distance
Aggr_Z_geodist.df <- data.frame(Mean_Z = with(Z_geodist.df, tapply(log(Z), INDEX = GeoDist_cat, mean)),
                                Var_Z = with(Z_geodist.df, tapply(log(Z), INDEX = GeoDist_cat, var)),
                                Dist_cat = paste0('Cat_', names(with(Z_geodist.df, tapply(log(Z), INDEX = GeoDist_cat, var)))))

#plot mean and variance of Z as a function of distance
ggplot(Aggr_Z_geodist.df, aes(x = Dist_cat, y = Mean_Z, group = 1)) +
  geom_line(lwd = 1) +
  geom_point(size = 4, col = 'green') +
  ylab('Mean log-Beta dissimilarity') + xlab('Distance class') +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 14))

ggplot(Aggr_Z_geodist.df, aes(x = Dist_cat, y = Var_Z, group = 1)) +
  geom_line(lwd = 1) +
  geom_point(size = 4, col = 'red') +
  ylab('Variance log-Beta dissimilarity') + xlab('Distance class') +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 14))

#plot marginal mean against marginal variance
#create .1 bins of Z and compute mean and variance within them
#this should replicate figure A.1 in White's supp material

Z_geodist.df$Z_bins <- with(Z_geodist.df, cut(Z, breaks = seq(0, 1, by = .1)))

Aggr_Z.df <- data.frame(Mean_Z = as.numeric(with(Z_geodist.df, tapply(Z, INDEX = Z_bins, mean))),
                        Var_Z = as.numeric(with(Z_geodist.df, tapply(log(Z), INDEX = Z_bins, var))),
                        Z_bin = names(with(Z_geodist.df, tapply(Z, INDEX = Z_bins, mean))))

ggplot(Aggr_Z.df, aes(x = Mean_Z, y = Var_Z, group = 1)) +
  geom_line(lwd = 1) +
  geom_point(size = 4, col = 'red') +
  ylab('Variance log-Beta dissimilarity') + xlab('Mean Beta dissimilarity') +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 14))


#------------------Generate regressors - splines

X_mat <- fcas_env_use

#define number of knots and degrees for building splines
Spl_deg <- 3
Spl_knots <- 1
Spl_df <- Spl_deg + Spl_knots

#formula to be used to generate basis functions
Spl_formula <- as.formula(paste("~ 0 +", paste(
  paste("iSpline(`", colnames(X_mat),"`,degree=", Spl_deg - 1 ,",df = ", Spl_df,
        " ,intercept = TRUE)", sep = ""), collapse = "+")))

#generate model.matrix of basis functions - number of basis per predictor depends on Spl_df
I_spl_basisfun <- model.matrix(Spl_formula, data = X_mat)

#create design matrix for spGDMM - note that the function computes Euclidean distances between sites' env values
X_for_GDM <- do.call(cbind, lapply(seq_len(ncol(I_spl_basisfun)), function(i) {
  dist_temp <- rdist(I_spl_basisfun[, i]) #compute Euclidean distance between pairs of observations
  vec_dist <- dist_temp[upper.tri(dist_temp)] #extract upper triangle of dist_temp matrix
  return(vec_dist)
  })) 

#check nrow(X_for_GDM) - (nrow(I_spl_basisfun)*(nrow(I_spl_basisfun) - 1))/2
nrow(X_for_GDM) #44551

#compute basis functions for geographical distances
Basis_geodist <- iSpline(Site_geodist, degree = (Spl_deg - 1), df = Spl_df, intercept = TRUE)

#cbind X_for_GDM and Basis_geodist
X_for_GDM <- cbind(X_for_GDM, Basis_geodist)

#delete Basis_geodist
rm(Basis_geodist)

#save ncol of X_for_GDM
N_col_XforGDM <- ncol(X_for_GDM)

#modify colnames of X_for_GDM
colnames(X_for_GDM) <- c(paste(rep(colnames(X_mat), each = Spl_df), "I", rep(1:Spl_df, times = ncol(X_mat)), sep = ""),
                         paste("dist", "I", 1:Spl_df, sep = ""))

#------------------Generate indices for pairing sites

#create a matrix including repetitions of each element in 1:nrow(bc_distmat) per col, each col has length equal to nrow(bc_distmat)
#1 - 2 - 3 ..
#1 - 2 - 3 ..
#.. .. ..  ..
#1 - 2 - 3 ..
tmp <- matrix(rep(1:nrow(bc_distmat), each = nrow(bc_distmat)), nrow = nrow(bc_distmat))
#extract upper triangle
col_ind <- tmp[upper.tri(tmp)] # 2 - 3, 3 - 4, 4, 4 - ..


#same as above, but now 1:nrow(bc_distmat) is repeated nrow(bc_distmat) and put in a matrix
#1 - 1 - 1 ..
#2 - 2 - 2 ..
#.. ..  ..
#nrow(bc_distmat)
tmp <- matrix(rep(1:nrow(bc_distmat), times = nrow(bc_distmat)), nrow = nrow(bc_distmat))
#extract upper triangle
row_ind <- tmp[upper.tri(tmp)] # 1 - 1, 2 - 1, 2, 3 - ..

#------------------Generate 'Kernel' matrix for Gaussian Process

#See appendix of White et al. 2023 for further info on parameters

#rho = 1/phi
#phi controls how quickly covariance among sites decays with distance
#here 10 is used for scaling - this parameter can be estimated in the model
Rho_fix <- max(geo_distmat)/10

#compute Kernel matrix
R_spat <- exp(-geo_distmat/Rho_fix)

#check if R_spat is postive-definite - it must be invertible!
all(eigen(R_spat)$val > 0) #TRUE

#Cholesky decomposition
R_chol <- t(chol(R_spat))

#take inverse of R_spat
R_inv <- solve(R_spat)

#check how covariance decays with distance
Rspat_offd_temp <- R_spat[upper.tri(R_spat)]

plot(Site_geodist, Rspat_offd_temp)

#------------------------------------Export input data to project for fitting spGDMM

#check date of object creation for reference
save(Obs_Z, X_for_GDM, N_col_XforGDM, row_ind, col_ind, Smp_size, N_sites, R_inv, Site_geodist,
     file = '/MOTIVATE/GDM_ForesteCasentinesi/spGDMM_fold/Data_for_spGDMM.RData')

