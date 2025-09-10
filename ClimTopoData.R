#extract climate, topographic and disturbance data to fit GDM

library(sf)
library(mapview)
library(easyclimate) #installed on Jan 9th 2025 - version 0.2.2
library(terra)
library(geosphere)
library(GGally)
library(ggplot2)
library(ggpubr)


#select column of EVA_forcas_meta that will be needed for the analyses

#check if some releve numbers are duplicated
sum(duplicated(EVA_forcas_meta$`TV2 relevé number`)) #0
sum(duplicated(EVA_forcas_meta$PlotObservationID)) #0

#subset columns of EVA_forcas_meta
EVA_forcas_meta <- EVA_forcas_meta[c('PlotObservationID', 'PlotID', 'TV2 relevé number', 'Date of recording', 'Relevé area (m²)', 'Altitude',
                                     'Locality', 'Expert System', 'ReSurvey plot (Y/N)', 'For EVA (Y/N)', 'Plot shape', 'Manipulate (y/n)', 'Type of manipulation',
                                     'lon', 'lat', 'Year')]

#rename some columns
colnames(EVA_forcas_meta)[c(3, 4, 5, 8, 11)] <- c('OriginalID', 'Date_recording', 'Plot_size', 'Expert_system', 'Plot_shape')

#check class of ID col
class(EVA_forcas_meta$OriginalID)

#coerce ID to character in EVA_forcas_meta
EVA_forcas_meta$OriginalID <- as.character(EVA_forcas_meta$OriginalID)

#check on missing sampling years
anyNA(EVA_forcas_meta$Date_recording) #F

#check duplicates
sum(duplicated(EVA_forcas_meta$OriginalID)) #0
 
#----pre-processing before extraction of climate data

#create EVA_forcas_climate, which will include climatic data
EVA_forcas_climate <- EVA_forcas_meta[c('OriginalID', 'lon', 'lat', 'Year')]

#rename OriginalID to ID to make it more general
colnames(EVA_forcas_climate)[1] <- 'ID'

#check range of years covered
range(EVA_forcas_climate$Year) #1985 2022

#we are unsure whether sampling year is real or rather an approximation for multiple sampling year
#(e.g., a group of plots was sampled over multiple years, yet a single year was indicated)
#for this reason, we download climatic data for a temporal window of the form (Year - Span/2):(Year + Span/2)
#span is the length of the temporal window - 1. If we choose span = 6, and Year = 2000, then climatic years are downloaded from 1997 to 2003
#we then aggregate climatic data over that period
#this avoids binning observations per sampling periods and assign observations to bins

#example
(2000 - 3):(2000 + 3)

#----------------------------------------extract climatic information for each plot

#see: https://verughub.github.io/easyclimate/reference/get_daily_climate.html
#https://verughub.github.io/easyclimate/articles/points-df-mat-sf.html
#https://verughub.github.io/easyclimate/articles/polygons-raster.html

#check_server() allows checking server is working correctly - note that get_daily_climate() has an arg check by default TRUE that does the same
check_server(climatic_var = "Tmin", year = 1950, verbose = T)
check_server(climatic_var = "Tmin", year = 1950, verbose = F)

#test example with 1 single plot (also previously ran using text_ex as an sf object)
#as a test, I select a plot that was sampled in 2022, so that the function has to adjust the period for download

test_ex <- EVA_forcas_climate[858, ]

test_span <- 6

#coords must be provided in longlat format
#temp: Celsius degrees
#prc: mm

#the plot was sampled in 2022, so the period for download should go from 2016:2022
test_ex_data <- get_daily_climate(coords = test_ex, climatic_var = c('Tmax', 'Tmin', 'Prcp'),
                                  period = 2016:2022,
                                  output = 'df')

#add Tmean manually
test_ex_data$Tmean <- with(test_ex_data, (Tmin + Tmax)/2)

#add year_climate and new_date manually
test_ex_data$new_date <- format(as.Date(test_ex_data$date, format = "%Y-%m-%d"), format = "%m/%d/%Y")
test_ex_data$year_cl <- format(as.Date(test_ex_data$date, format = "%Y-%m-%d"), format = "%Y")

#check on period
unique(test_ex_data$year_cl)

#function to download data and aggregate daily records over selected period
#I wrapped 'get_daily_climate' within try() to avoid errors with data download breaking the process
#Note that if Year > 2019, then the function needs to adjust the period of download of the data, as 2022 is the last year for
#which climate data are available
#Additionally, note that lower boundary is not checked as well - the function will simply stop if (sampl_yr - span/2) < 1950

get_avg_climate <- function(x, span, cl_lim = 2022L) {
  require(easyclimate)
  require(ClimInd)
  #extract year of sampling
  sampl_yr <- x[['Year']]
  #check that year is covered by climatic data
  if(isTRUE(sampl_yr > cl_lim)) stop('Year larger than cl_lim')
  #set half span
  half_span <- (span/2)
  #if span is an odd number, adjust it with floor
  if(isTRUE((span %% 2) != 0)) {
    half_span <- floor(half_span)
    message(paste0('Span adjusted to: ', half_span))
  }
  #check lower boundary
  if(isTRUE((sampl_yr - half_span) < 1950L)) stop('Lower limit is before 1950')
  #check upper boundary and compute period over which data are downloaded
  up_lim <- (sampl_yr + half_span)
  if(isTRUE(up_lim > cl_lim)) {
    #compute number of years above up_lim
    diff_yrs <- (up_lim - cl_lim)
    #exceeding years are assigned to lower boundary
    lw_lim <- sampl_yr - (half_span + diff_yrs)
    #up_lim is cl_lim
    up_lim <- cl_lim
    #set period
    prd_dwn <- (lw_lim):(up_lim)
  } else {
    lw_lim <- (sampl_yr - half_span)
    prd_dwn <- lw_lim:up_lim
  }
  #download data
  cl_data <- try(get_daily_climate(coords = x, climatic_var = c('Tmax', 'Tmin', 'Prcp'), period = prd_dwn, output = 'df', version = 4))
  #if cl_data inherits 'try-error' returns NA -> this was taken from Advanced R 1st Ed
  if(inherits(cl_data, what = "try-error")) {
    return(NA)
  }
  cl_data <- cl_data[c('date', 'Tmax', 'Tmin', 'Prcp')]
  #add mean temperature
  cl_data$Tmean <- (cl_data[["Tmin"]] + cl_data[["Tmax"]])/2
  #compute a set of climatic indices (bioclimatic variables following WorldClim + info on growing season length and precipitation)
  cl_indices <- compute_climind(dtf = cl_data)
  #re-attach information on the plot
  cl_output <- data.frame(x[c('ID', 'Year', 'lon', 'lat')], period_from = lw_lim, period_to = up_lim, cl_indices)
  return(cl_output)
}

#these checks were done before re-running analyses with EVA from the beginning and before dropping out plots sampled before 1985
#quick sanity check - the function should not work if the period < 1950
#get_avg_climate(x = EVA_forcas_climate[4, ], span = test_span) #good

#Year = 1954
#get_avg_climate(x = EVA_forcas_climate[928, ], span = test_span) #this goes
#get_avg_climate(x = EVA_forcas_climate[928, ], span = 9) #this goes
#get_avg_climate(x = EVA_forcas_climate[928, ], span = 10) #good - this returns an error

#compare against test_ex_data

test_ex_data2 <- get_avg_climate(x = test_ex, span = test_span)

#all looks pretty good - climatic and bioclimatic vars
all(all.equal(mean(tapply(test_ex_data$Tmean, test_ex_data$year_cl, mean)), test_ex_data2$tmean_avg), 
    all.equal(mean(tapply(test_ex_data$Tmax, test_ex_data$year_cl, mean)), test_ex_data2$tmax_avg), 
    all.equal(mean(tapply(test_ex_data$Tmin, test_ex_data$year_cl, mean)), test_ex_data2$tmin_avg), 
    all.equal(mean(tapply(test_ex_data$Prcp, test_ex_data$year_cl, sum)), test_ex_data2$prcp_avg), 
    all.equal(mean(bio4(data = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date))), test_ex_data2$bio4_avg),
    all.equal(mean(bio5(data = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date), tmax = setNames(test_ex_data$Tmax, nm = test_ex_data$new_date))), test_ex_data2$bio5_avg),
    all.equal(mean(bio6(data = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date), tmin = setNames(test_ex_data$Tmin, nm = test_ex_data$new_date))), test_ex_data2$bio6_avg),
    all.equal(mean(bio8(pr = setNames(test_ex_data$Prcp, nm = test_ex_data$new_date), taverage = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date))), test_ex_data2$bio8_avg),
    all.equal(mean(bio10(data = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date))), test_ex_data2$bio10_avg),
    all.equal(mean(bio14(data = setNames(test_ex_data$Prcp, nm = test_ex_data$new_date))), test_ex_data2$bio14_avg),
    all.equal(mean(bio18(pr = setNames(test_ex_data$Prcp, nm = test_ex_data$new_date), taverage = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date))), test_ex_data2$bio18_avg),
    all.equal(mean(gsl(data = setNames(test_ex_data$Tmean, nm = test_ex_data$new_date))), test_ex_data2$gsl_avg),
    all.equal(mean(gsr(data = setNames(test_ex_data$Prcp, nm = test_ex_data$new_date))), test_ex_data2$gsr_avg)) # TRUE (!)

#run it for all plots

#download process started at:
Sys.time() #"2025-02-17 15:22:28 CET"

#get climate data for all sampling plots
FCas_climate_data <- lapply(seq_len(nrow(EVA_forcas_climate)), function(id_r) {
  res <- get_avg_climate(x = EVA_forcas_climate[id_r, ], span = 6L)
  if(is.data.frame(res)) {
    message(paste('Plot ', res[['ID']], ' done!'))
  } else {
    message(paste('Plot returned NA'))
  }
  return(res)
})


#debugonce(get_avg_climate)

#and ended at:
Sys.time() #"2025-02-19 13:56:13 CET" (approx.)


#check if in some of the data.frames there are NAs
sum(sapply(Filter(is.data.frame, FCas_climate_data), anyNA)) #0

#check for which plots it was not possible to retrieve data
EVA_forcas_climate$ID[which(sapply(FCas_climate_data, is.logical))]
length(EVA_forcas_climate$ID[which(sapply(FCas_climate_data, is.logical))]) #13

#quick view of some climatic parameters
plot(tmean_avg ~ Year, data = do.call(rbind, Filter(is.data.frame, FCas_climate_data)))
plot(prcp_avg ~ Year, data = do.call(rbind, Filter(is.data.frame, FCas_climate_data)))

#save a copy of FCas_climate_data for future recovery of the data
#this version has NAs for 13 plots
save(FCas_climate_data, file = 'climatic_data_per_plot.RData')

#name FCas_climate_data with plot ID
names(FCas_climate_data) <- EVA_forcas_climate$ID

#extract plot ID for which I miss climatic data
plot_missing_climate <- names(Filter(is.logical, FCas_climate_data)) 

#subset EVA_forcas_climate to include only plots for which I miss data
plot_missing_climate <- EVA_forcas_climate[EVA_forcas_climate$ID %in% plot_missing_climate, ]  

#try to retrieve climatic data for missing plots
Fcas_cl_missing <- lapply(seq_len(nrow(plot_missing_climate)), function(id_r) {
  res <- get_avg_climate(x = plot_missing_climate[id_r, ], span = 6L)
  return(res)
  })


#check how many NAs
Filter(is.logical, Fcas_cl_missing) #0

#check how many NAs in data.frame
sum(sapply(Fcas_cl_missing, anyNA)) #

#assign names
names(Fcas_cl_missing) <- plot_missing_climate$ID

#assign missing elements to metad_climate
FCas_climate_data[names(Fcas_cl_missing)]

for(nm in names(Fcas_cl_missing)) {
  FCas_climate_data[[nm]] <- Fcas_cl_missing[[nm]]
}

rm(nm)

#sanity check - list names match with plot ID
sum(sapply(names(FCas_climate_data), function(nm) nm == FCas_climate_data[[nm]][['ID']]))

#rbind the list to get all info in a data.frame
Fcas_climate.df <- do.call(rbind, FCas_climate_data)

#---------quick overview on correlation among climatic parameters

#correlation matrix
cor(Fcas_climate.df[7:27])

#plot of correlation matrix
ggcorr(data = Fcas_climate.df[7:27], geom = 'text')

#biplot of PCA
biplot(prcomp(Fcas_climate.df[7:27], center = T, scale. = T))

#create PCA object
Fcas_clim_pca <- prcomp(Fcas_climate.df[7:27], center = T, scale. = T)

#check variance explained by the 2 axis
summary(Fcas_clim_pca) #73%

#get correlation between variables and first 2 axes
Fcas_pca_cors <- do.call(cbind, lapply(seq_len(ncol(Fcas_clim_pca$rotation)), function(i) Fcas_clim_pca$rotation[, i] * Fcas_clim_pca$sdev[i]))

#select correlations related to the first 3 axes
Fcas_pca_cors <- as.data.frame(Fcas_pca_cors[, c(1, 2, 3)])

#add column with parameter name
Fcas_pca_cors$Param_name <- factor(rownames(Fcas_pca_cors), levels = rownames(Fcas_pca_cors))

#make it a long dataframe
Fcas_pca_cors <- data.frame(Cor = c(Fcas_pca_cors$V1, Fcas_pca_cors$V2, Fcas_pca_cors$V3),
                            PC_ax = rep(c("PC1", "PC2", "PC3"), each = nrow(Fcas_pca_cors)),
                            Param_name = rep(Fcas_pca_cors$Param_name, times = 3))

ggplot(Fcas_pca_cors, aes(x = Param_name, y = Cor, col = PC_ax)) +
  geom_hline(yintercept = c(-0.5, 0.5), col = 'red') +
  geom_hline(yintercept = 0, col = 'black') +
  geom_point(size = 3) +
  facet_wrap(~ PC_ax) +
  theme_pubr() +
  theme(legend.position = 'top', axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))

#get name of bioclim parameters showing high correlation (> |r = .7|) with PC1/2/3

cl_indices_per_ax <- lapply(unique(Fcas_pca_cors$PC_ax), function(ax) {
  dtf <- Fcas_pca_cors[Fcas_pca_cors$PC_ax == ax, ]
  cl_pr_id <- which(abs(dtf[['Cor']]) > .7)
  pr_nm <- as.character(dtf[cl_pr_id, "Param_name"])
  return(pr_nm)
  })

names(cl_indices_per_ax) <- unique(Fcas_pca_cors$PC_ax)

#PC1 - all neg. correlations: high values of PC1 -> colder conditions; low values of PC1 -> warmer conditions
Fcas_pca_cors[Fcas_pca_cors$Param_name %in% cl_indices_per_ax$PC1 & Fcas_pca_cors$PC_ax == "PC1", ]

#PC2 - all pos. correlations: high values of PC2 -> wetter conditions; low values of PC2 -> drier conditions 
Fcas_pca_cors[Fcas_pca_cors$Param_name %in% cl_indices_per_ax$PC2 & Fcas_pca_cors$PC_ax == "PC2", ]

#PC3 - all pos. correlations: high values of PC3 -> higher prcp wettest month, coldest quarter and seasonality + min temp coldest month
Fcas_pca_cors[Fcas_pca_cors$Param_name %in% cl_indices_per_ax$PC3 & Fcas_pca_cors$PC_ax == "PC3", ]

#----------------------------------------extract topographic information for each plot

#import topographic layers - Geomorph90m - Amatulli et al. 2020
#for info on dates of download and other things, see README in the TopographicData folder

slope_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/slope_90M_n30e000/slope_90M_n40e010.tif")
north_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/northness_90M_n30e000/northness_90M_n40e010.tif")
east_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/eastness_90M_n30e000/eastness_90M_n40e010.tif")
tri_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/tri_90M_n30e000/tri_90M_n40e010.tif")

terra::compareGeom(slope_tile, north_tile)
terra::compareGeom(north_tile, east_tile)
terra::compareGeom(east_tile, tri_tile)

topo_stack <- c(slope_tile, north_tile, east_tile, tri_tile)

names(topo_stack) <- c('slope', 'north', 'east', 'tri')

#extract value of topographic points at each location

topo_values <- extract(x = topo_stack, y = EVA_forcas_meta[c('lon', 'lat')])

#get rid of ID col

topo_values$ID <- NULL

#attach plot ID back

topo_values <- data.frame(ID = EVA_forcas_meta[['OriginalID']], topo_values)

#check if we get the same using SpatVect (yes)
#spvec_obj <- vect(as.matrix(EVA_forcas_meta[c('lon', 'lat')]), atts = EVA_forcas_meta[c('OriginalID')], crs = "+proj=longlat +datum=WGS84")
#plot(spvec_obj)
#spvec_topo_vals <- extract(x = topo_stack, y = spvec_obj)
#all.equal(spvec_topo_vals$slope, topo_values$slope) #T
#rm(spvec_obj, spvec_topo_vals)

#check correlation between topo_values
cor(topo_values[-1]) #slope and tri are essentially the same variable
plot(topo_values$slope, topo_values$tri)

#check range of variables
lapply(topo_values[-1], range) #looks correct

#check NAs
anyNA(topo_values) #F


#----------------------------------------join climatic and topographic data with EVA metadata

#check duplicated IDs in all datasets
sum(duplicated(Fcas_climate.df$ID)) #0
sum(duplicated(topo_values$ID)) #0
sum(duplicated(EVA_forcas_meta$OriginalID)) #0

#join climatic data
class(EVA_forcas_meta$OriginalID); class(Fcas_climate.df$ID)

#check correspondence of plot IDs
setdiff(EVA_forcas_meta$OriginalID, Fcas_climate.df$ID) #0
setdiff(Fcas_climate.df$ID, EVA_forcas_meta$OriginalID) #0

#check correspondence of Year
identical(EVA_forcas_meta$Year, Fcas_climate.df$Year) #T

EVA_forcas_meta <- dplyr::left_join(x = EVA_forcas_meta, y = Fcas_climate.df[!colnames(Fcas_climate.df) %in% c('lon', 'lat', 'Year')],
                                    by = c('OriginalID' = 'ID'))

class(EVA_forcas_meta)

#join topographic values
class(topo_values$ID)

#check correspondence of plot IDs
setdiff(EVA_forcas_meta$OriginalID, topo_values$ID) #0
setdiff(topo_values$ID, EVA_forcas_meta$OriginalID) #0

EVA_forcas_meta <- dplyr::left_join(x = EVA_forcas_meta, y = topo_values, by = c('OriginalID' = 'ID'))

#coerce Altitude to double
EVA_forcas_meta$Altitude <- as.double(EVA_forcas_meta$Altitude)

#----------------------------------------exploratory analyses

#check variability in Habitat type
unique(EVA_forcas_meta$Expert_system) #tilde symbol means unclassified

sort(table(EVA_forcas_meta$Expert_system), decreasing = T) #180 unclassified habitats

#table of habitat types
Hab_table <- as.data.frame(table(EVA_forcas_meta$Expert_system))

#re-order var1
Hab_table$Var1 <- factor(Hab_table$Var1, levels = as.character(Hab_table[order(Hab_table$Freq, decreasing = T), 'Var1']))

#keep only habitats with at least 20 plots
Hab_table.hr <- Hab_table[Hab_table$Freq >= 20, ]

#coarse attribution of the habitat
Hab_table.hr$Hab_typ <- ifelse(substr(as.character(Hab_table.hr$Var1), 1, 1) == "R", 'Grassland',
                               ifelse(substr(as.character(Hab_table.hr$Var1), 1, 1) == "T", 'Forest', 'Unclassified'))


Hab_typ_plot <- ggplot(Hab_table.hr, aes(x = Var1, y = Freq, fill = Hab_typ)) +
  geom_col() +
  scale_fill_manual(values = c('Forest' = 'darkgreen', 'Grassland' = 'lightgreen', 'Unclassified' = 'grey'), name = 'Habitat Type') +
  ylab('Count') + xlab(NULL) +
  theme_pubr() +
  theme(axis.title = element_text(size = 16), axis.text.x.bottom = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))


#check range of altitude, slope and tri
Slope_plot <- ggplot(EVA_forcas_meta, aes(y = slope)) +
  geom_boxplot() +
  ggtitle('Slope') + ylab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 18),
        title = element_text(size = 20))

#note that Altitude has 4 NAs
Alt_plot <- ggplot(EVA_forcas_meta, aes(y = Altitude)) +
  geom_boxplot() +
  ggtitle('Altitude') + ylab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 18),
        title = element_text(size = 20))

Tri_plot <- ggplot(EVA_forcas_meta, aes(y = tri)) +
  geom_boxplot() +
  ggtitle('Topogr. Roughness Index') + ylab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 18),
        title = element_text(size = 20))

#check range of some bioclimatic parameters
Tmean_plot <- ggplot(EVA_forcas_meta, aes(y = tmean_avg)) +
  geom_boxplot() +
  ggtitle('Mean Temp.') + ylab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 18),
        title = element_text(size = 20))

Precip_plot <- ggplot(EVA_forcas_meta, aes(y = prcp_avg)) +
  geom_boxplot() +
  ggtitle('Mean Annual Precip.') + ylab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks.x.bottom = element_blank(),
        axis.text.y.left = element_text(size = 18),
        title = element_text(size = 20))


ggarrange(Alt_plot, Slope_plot, Tmean_plot, Precip_plot, nrow = 2, ncol = 2)

#----------------------------------------import disturbance data

#import time-series of TCW data - values represent median during the vegetative season (excluding clouds, snow and shadow)
#this version of the data substitutes a previous version with weird offset in TCW values since 2013 (so, from Landsat 8)
Fcas_dist.df <- read.csv(file = '/MOTIVATE/GDM_ForesteCasentinesi/DisturbanceData/TCW_LandTrendr_per_plot_per_year.csv', header = T,
                         sep = ',', dec = '.')

#check number of plots
length(unique(Fcas_dist.df$PlotID)) #943

all(unique(Fcas_dist.df$PlotID) %in% EVA_forcas_meta$PlotObservationID) #T

setdiff(unique(Fcas_dist.df$PlotID), EVA_forcas_meta$PlotObservationID) #0
setdiff(EVA_forcas_meta$PlotObservationID, unique(Fcas_dist.df$PlotID)) #0

#check temporal extent of TCW
range(Fcas_dist.df$year) #1984-2024

943*length(seq(1984, 2024, 1)) #38663

all(table(Fcas_dist.df$PlotID) == 41) #T

#rename cols
colnames(Fcas_dist.df)[3] <- 'Year'

#coerce Year to factor
Fcas_dist.df$Year_fct <- factor(Fcas_dist.df$Year, levels = sort(unique(Fcas_dist.df$Year)))

#check NA
sum(is.na(Fcas_dist.df$TCW)) #0

#exploratory
ggplot(Fcas_dist.df, aes(x = Year_fct, y = TCW)) +
  geom_point() +
  xlab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1))

#check spatial pattern TCW
explo_tcw_trend <- dplyr::left_join(x = Fcas_dist.df, y = EVA_forcas_meta[c('PlotID', 'lon', 'lat')], by = "PlotID")

#mapview(st_as_sf(x = explo_tcw_trend, coords = c('lon', 'lat'), crs = 4326))

ggplot(explo_tcw_trend, aes(x = lon, y = lat, col = TCW)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_pubr()

#add Year chr EVA_forcas_meta
EVA_forcas_meta$Year_chr <- as.character(EVA_forcas_meta$Year)

#add Year chr Fcas_dist.df
Fcas_dist.df$Year_chr <- as.character(Fcas_dist.df$Year)

#join disturbance data to EVA_forcas_meta

#check - also checked if some values (from head() and tail()) of TCW for years in EVA_forcas_meta would match those in Fcas_dist.df, which includes whole time-series, and they match
#prova <- dplyr::left_join(x = EVA_forcas_meta, y = Fcas_dist.df[c('PlotID', 'TCW', 'Year_chr')], by = c('PlotID', 'Year_chr')) 
#identical(prova[-ncol(prova)], EVA_forcas_meta) #TRUE

EVA_forcas_meta <- dplyr::left_join(x = EVA_forcas_meta, y = Fcas_dist.df[c('PlotID', 'TCW', 'Year_chr')], by = c('PlotID', 'Year_chr'))

#----------------------------------------import aspect and compute solar radiation

#import topographic layers - Geomorph90m - Amatulli et al. 2020
#for info on dates of download and other things, see README in the TopographicData folder

#assuming aspect is in degrees as it goes from 0 to about 360
aspect_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/aspect_90M_n40e010.tif")

#compare aspect_tile with the other topographic layers
compareGeom(aspect_tile, tri_tile) #T

#extract aspect values at EVA_forcas_meta coordinates
aspect_values <- terra::extract(aspect_tile, y = EVA_forcas_meta[c('lon', 'lat')])

anyNA(aspect_values) #F
range(aspect_values$aspect_90M_n40e010) #0.2685165 359.5325623

#drop ID, modify name of aspect values' column and add plot ID
aspect_values <- data.frame(OriginalID = EVA_forcas_meta[['OriginalID']], aspect = aspect_values$aspect_90M_n40e010)

#join to EVA_forcas_meta
#prova <- dplyr::left_join(x = EVA_forcas_meta, y = aspect_values, by = 'OriginalID')
#identical(prova[-length(prova)], EVA_forcas_meta) #T
EVA_forcas_meta <- dplyr::left_join(x = EVA_forcas_meta, y = aspect_values, by = 'OriginalID')

#use slope, aspect and latitude to compute Solar Radiation

#import function SolarRad from: https://github.com/fmsabatini/SolarRad/blob/main/SolarRad.R

SolarRad <- function(Slope, Aspect, Latitude, Exp=FALSE){
  if(Slope==0) Aspect=999
  folded_aspect <- 180-abs(Aspect-180)
  rlatitude <- Latitude*pi/180                #transform to radians
  rslope    <- Slope*pi/180
  rfolded_aspect <- folded_aspect*pi/180
  SolarRad <- -1.467+1.582*cos(rlatitude)*cos(rslope) - 1.5*cos(rfolded_aspect)*sin(rslope)*sin(rlatitude) -
    0.262*sin(rlatitude)*sin(rslope) + 0.607*sin(rfolded_aspect)*sin(rslope)
  if(Exp==T) SolarRad <- exp(SolarRad)
  
  folded_aspectNW_SE <- abs(180-abs(Aspect-225))  #transform to radians centered on the SW/NE
  rfolded_aspectNW_SW <- folded_aspectNW_SE*pi/180
  HeatLoad <- -1.467+1.582*cos(rlatitude)*cos(rslope) - 1.5*cos(folded_aspectNW_SE)*sin(rslope)*sin(rlatitude) -
    0.262*sin(rlatitude)*sin(rslope) + 0.607*sin(folded_aspectNW_SE)*sin(rslope)
  if(Exp==T) HeatLoad <- exp(HeatLoad)
  return(data.frame(SolarRad, HeatLoad))
}

SolarRadiation_df <- do.call(rbind, lapply(1:nrow(EVA_forcas_meta), function(i) {
  res <- SolarRad(Slope = EVA_forcas_meta$slope[i], Aspect = EVA_forcas_meta$aspect[i], Latitude = EVA_forcas_meta$lat[i], Exp = F)
  res$OriginalID <- EVA_forcas_meta$OriginalID[i]
  return(res)
  }))

anyNA(SolarRadiation_df) #F

#join to EVA_forcas_meta
#prova <- dplyr::left_join(x = EVA_forcas_meta, y = SolarRadiation_df, by = 'OriginalID')
#identical(prova[-c(length(prova)-1, length(prova))], EVA_forcas_meta) #T
EVA_forcas_meta <- dplyr::left_join(x = EVA_forcas_meta, y = SolarRadiation_df, by = 'OriginalID')

#check correlation with topographic variables
#in the end it's probably better if we use heat load, which is nearly orthogonal with other predictors
#and it has been used in other studies as a proxy of microclimate conditions: https://onlinelibrary.wiley.com/doi/full/10.1111/jvs.13305
cor(EVA_forcas_meta[c('slope', 'north', 'east', 'tri', 'SolarRad', 'HeatLoad')])

#----------------------------------------import environmental raster layers and bring them all to the same resolution and extent of the coarsest layer

#the idea is to create a stack of raster layers to be used to map the RGB of estimated warping functions across the study area
#to do that, we need to bring all layers to the same extent, resolution and so on
#at the moment, the environmental layer(s) with the coarsest resolution is temperature (1 km at the Equator)
#so the topographic and TCW layers have to be modified to match that same resolution, extent and so on
#however, as the RGB (as well as the maps of spatial random effects) will be used for visualization only
#we'll force the layers to have a resolution that allows showing spatial patterns, such as 250 m
#this means that the climatic layers will be forced to have a better resolution than they actually have
#important note: the resolution of the layer cannot be too fine as this would imply sampling GP realizations
#from an high dimension mvnorm

#import TCW layers
TCW_tile_1984 <- rast("TCWData/TCW_1984_corrext.tif")
TCW_tile_2020 <- rast("TCWData/TCW_2020_corrext.tif")

#import eastness, slope and aspect
slope_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/slope_90M_n30e000/slope_90M_n40e010.tif")
east_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/eastness_90M_n30e000/eastness_90M_n40e010.tif")
aspect_tile <- rast(x = "C:/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/TopoLayers/aspect_90M_n40e010.tif")

#import temperature and precipitation (and apply scale to report them to correct range of values)
temp_tile_1984 <- rast(x = "C:/MOTIVATE/GDM_EuropeanEcoregions/EasyClimateData/v4data_from_ftpzilla/tavg/DownscaledTavg1984YearlyAvg_cogeo.tif")
temp_tile_1984 <- temp_tile_1984/100
temp_tile_2020 <- rast(x = "C:/MOTIVATE/GDM_EuropeanEcoregions/EasyClimateData/v4data_from_ftpzilla/tavg/DownscaledTavg2020YearlyAvg_cogeo.tif")
temp_tile_2020 <- temp_tile_2020/100

prcp_tile_1984 <- rast(x = "C:/MOTIVATE/GDM_EuropeanEcoregions/EasyClimateData/v4data_from_ftpzilla/prec/DownscaledPrcp1984YearlySum_cogeo.tif")
prcp_tile_1984 <- prcp_tile_1984/100
prcp_tile_2020 <- rast(x = "C:/MOTIVATE/GDM_EuropeanEcoregions/EasyClimateData/v4data_from_ftpzilla/prec/DownscaledPrcp2020YearlySum_cogeo.tif")
prcp_tile_2020 <- prcp_tile_2020/100

#crop all layers at the extent of FCas_shp
#check if ext() returns a SpatExtent if x = sf object
class(ext(FCas_shp)) #Yes

TCW_tile_1984 <- crop(TCW_tile_1984, y = FCas_shp)
TCW_tile_2020 <- crop(TCW_tile_2020, y = FCas_shp)

slope_tile <- crop(slope_tile, y = FCas_shp)
east_tile <- crop(east_tile, y = FCas_shp)
aspect_tile <- crop(aspect_tile, y = FCas_shp)

temp_tile_1984 <- crop(temp_tile_1984, y = FCas_shp)
temp_tile_2020 <- crop(temp_tile_2020, y = FCas_shp)

prcp_tile_1984 <- crop(prcp_tile_1984, y = FCas_shp)
prcp_tile_2020 <- crop(prcp_tile_2020, y = FCas_shp)

#----use aspect and slope to compute the heat load layer
heatload_tile_df <- as.data.frame(c(aspect_tile, slope_tile), xy = T)

nrow(heatload_tile_df) == ncell(aspect_tile)

#update col names
colnames(heatload_tile_df)[c(3, 4)] <- c('aspect', 'slope')

#check NAs
anyNA(heatload_tile_df[c(3, 4)]) #F

#compute heat load on whole layer
heatload_tile_df$heatload <- sapply(seq_len(nrow(heatload_tile_df)), function(i) {
  asp_val <- heatload_tile_df[i, 'aspect']
  sl_val <- heatload_tile_df[i, 'slope']
  lat_val <- heatload_tile_df[i, 'y']
  res <- SolarRad(Slope = sl_val, Aspect = asp_val, Latitude = lat_val, Exp = F)
  res <- res[['HeatLoad']] 
  return(res)
})

#create heat load layer
heatload_tile <- terra::rasterize(x = heatload_tile_df[c(1, 2)], y = aspect_tile, values = heatload_tile_df$heatload)

#check geometry comparison
compareGeom(heatload_tile, aspect_tile) #T

#----create template layer

#create template layer in EPSG:32632
temp_layer <- rast(crs = 'epsg:32632', extent = ext(st_transform(FCas_shp, crs = 32632)))

#set resolution to 400 m, so that sampling from mvrnorm is feasible
res(temp_layer) <- c(400, 400)

#give it fake values values (not strictly needed)
temp_layer[] <- 1

#project temp_layer to epsg:4326
temp_layer <- project(x = temp_layer, y = 'epsg:4326', method = 'bilinear')

#check what's the new resolution (from longlat to m)
temp_layer_df <- as.data.frame(temp_layer, xy = T)

#min longitude and mean lat: 329.8845
distm(rbind(c(min(temp_layer_df$x), mean(temp_layer_df$y)), c(min(temp_layer_df$x) + res(temp_layer)[1], mean(temp_layer_df$y))))

#max and mean lat: 329.8845
distm(rbind(c(max(temp_layer_df$x), mean(temp_layer_df$y)), c(max(temp_layer_df$x) + res(temp_layer)[1], mean(temp_layer_df$y))))

#mean lat at mean long (distortion should not affect lat at all over long): 456.0014
distm(rbind(c(mean(temp_layer_df$x), mean(temp_layer_df$y)), c(mean(temp_layer_df$x), mean(temp_layer_df$y) + res(temp_layer)[1])))

#check how the ney layer overlaps with FCas perimeter
plot(temp_layer)
plot(st_geometry(FCas_shp), add = T)

#----report all layers to the same resolution as temp_layer

#small differences
lapply(list(TCW_tile_1984, TCW_tile_2020, slope_tile, east_tile, heatload_tile, temp_tile_1984, prcp_tile_1984, temp_layer), ext)

#check on output of using resample
#compareGeom(resample(x = TCW_tile_1984, y = temp_tile_1984, method = 'average'),
#           temp_tile_1984)

#I'm resampling all layers, including the climatic ones, but it could be useful to adopt another strategy for the climatic layers,
#whose resolution is way larger than that of temp_layer [TO CHECK!]
env_stack <- rast(lapply(list(temp_tile_1984, temp_tile_2020, prcp_tile_1984, prcp_tile_2020,
                              TCW_tile_1984, TCW_tile_2020, slope_tile, heatload_tile), function(i, lyr = temp_layer) {
  res <- terra::resample(i, y = lyr, method = 'average')
  return(res)
}))

names(env_stack) <- c("Tavg_1984", "Tavg_2020", "Prcp_1984", "Prcp_2020",
                      "TCW_1984", "TCW_2020", "slope", "HeatLoad")

#combine with temperature and precipitation layer
#env_stack <- c(temp_tile_1984, temp_tile_2020, prcp_tile_1984, prcp_tile_2020, env_stack)

#names(env_stack)[1:4] <- c("temp_1984", "temp_2020", "prcp_1984", "prcp_2020")

#mask the layers to overlap the spatial extent of the Foreste Casentinesi
env_stack <- mask(env_stack, mask = vect(FCas_shp), touches = T)

#coerce the env_stack to a data.frame for further calculations - see DataForGDM
env_stack_df <- as.data.frame(env_stack, xy = T, na.rm = T)


#---------------------------------------saves data for presentations

#these files are not anymore exported at this stage - see DataForGDM

#last save: 26/02/2025 -> BECAUSE workshop and test spGDMM
#save(EVA_forcas_meta, file = 'C://MOTIVATE/Talks_and_presentations/Data/EVA_forcas_data.RData')
#save(EVA_forcas_meta, file = 'C://spGDMM_repo/spGDMM-code-main/fcas_test_data/EVA_forcas_data.RData')

#last save: 26/02/2025 -> test spGDMM
#save(EVA_forcas_species, file = 'C://spGDMM_repo/spGDMM-code-main/fcas_test_data/EVA_forcas_species.RData')


