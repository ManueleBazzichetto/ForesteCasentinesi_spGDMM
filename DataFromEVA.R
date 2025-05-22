##import EVA dataset and select data of the Foreste Casentinesi

library(data.table)
library(sf)
library(mapview)

#import EVA-formatted data from the beginning,
#so that we can then generalize the whole workflow to other regions

#this shapefile of the Foreste Casentinesi can be used to filter observations that will be used for the analyses
FCas_shp <- st_read(dsn = '/MOTIVATE/GDM_ForesteCasentinesi/Shape_FCas/Casentinesi_Boundary/WDPA_WDOECM_Mar2024_Public_64511_shp-polygons.shp')

nrow(FCas_shp) #1

mapview(FCas_shp)

#this EVA version was downloaded on Feb, 3rd, 2025. See the README in MOTIVATE/EVA_dataset for further info

EVA_full <- fread(file = 'C://MOTIVATE/EVA_dataset/200_MOTIVATE20240412_taxa_concepts_unified.csv', data.table = FALSE)

#this EVA metadata version was downloaded on Feb, 3rd, 2025. See the README in MOTIVATE/EVA_dataset for further info

EVA_meta <- fread(file = 'C://MOTIVATE/EVA_dataset/200_MOTIVATE20240412_header_notJUICE_with_precise_coordinates.csv', data.table = FALSE)

#check
#sapply(list(EVA_full, EVA_meta), class) #both are data.frames

#number of (unique) plots

length(unique(EVA_full$PlotObservationID)) #2385767

length(unique(EVA_meta$PlotObservationID)) #2387117

#number of resurveyed plots
sum(EVA_meta$`ReSurvey plot (Y/N)` == 'Y') #425310

#check if Italy is in
'Italy' %in% EVA_meta$Country

#check how many plots there are in Italy
sum(EVA_meta$Country == 'Italy') #99191

#filter records in Italy
EVA_it <- EVA_meta[EVA_meta$Country == 'Italy', ]

#get rid of EVA_meta
rm(EVA_meta)

#quickly check plot location - there's a plot in the Antarctic region - this will be filtered out later
mapview(st_as_sf(EVA_it[complete.cases(EVA_it[c('Longitude', 'Latitude')]), ], coords = c('Longitude', 'Latitude'), crs = 4326))

#check locality names
unique(EVA_it$Locality)
grep(pattern = 'Foreste', x = EVA_it$Locality, value = T) #"Foreste Casentinesi - Romagna"
unique(grep(pattern = 'Foreste Casentinesi', x = EVA_it$Locality, value = T)) #"Foreste Casentinesi - Romagna" / "Foreste Casentinesi - Toscana"  

#check how many data points for 'Foreste Casentinesi - Romagna/Toscana'
nrow(EVA_it[EVA_it$Locality %in% c("Foreste Casentinesi - Romagna", "Foreste Casentinesi - Toscana"), ]) #318

#crop plots overlapping with the shapefile of the Foreste Casentinesi

#drop rows with missing coords
EVA_it <- EVA_it[complete.cases(EVA_it[c('Longitude', 'Latitude')]), ]

#coerce EVA_it to sf obj
EVA_it.sp <- st_as_sf(x = EVA_it, coords = c('Longitude', 'Latitude'), crs = 4326)

#crop EVA_it.sp - note that this will extract plots within the bbox (st_crop != mask)
EVA_forcas_meta <- st_crop(x = EVA_it.sp, y = FCas_shp) 

#check - ok
setdiff(colnames(st_drop_geometry(EVA_forcas_meta)), colnames(EVA_it))
setdiff(colnames(EVA_it), colnames(st_drop_geometry(EVA_forcas_meta)))
identical(colnames(st_drop_geometry(EVA_forcas_meta)), colnames(EVA_it[!colnames(EVA_it) %in% c("Longitude", "Latitude")]))
identical(st_drop_geometry(EVA_forcas_meta), EVA_it[EVA_it$PlotObservationID %in% EVA_forcas_meta$PlotObservationID, !colnames(EVA_it) %in% c("Longitude", "Latitude")]) #T

#check how many plots were resurveyed 
sum(EVA_forcas_meta$`ReSurvey plot (Y/N)` == "Y") #174

#check duplicates
sum(duplicated(EVA_forcas_meta$PlotObservationID)) #good

#drop geometry of EVA_forcas_meta
EVA_forcas_meta$lon <- st_coordinates(EVA_forcas_meta)[, 1]
EVA_forcas_meta$lat <- st_coordinates(EVA_forcas_meta)[, 2]
EVA_forcas_meta <- st_drop_geometry(EVA_forcas_meta)

#----check on missing sampling years

#create SamplYear col - note that there are cases of empty 'Date of recording' ('')
unique(EVA_forcas_meta$`Date of recording`)
sum(EVA_forcas_meta$`Date of recording` == '') #46

#have a look at data for which sampling year is missing
EVA_forcas_meta[EVA_forcas_meta$`Date of recording` == '', ]

#check how substr handles these cases
substr(EVA_forcas_meta$`Date of recording`[EVA_forcas_meta$`Date of recording` == ''], start = 7, stop = 10)
#coercing them to integer returns NAs
as.integer('') #good, we can drop these cases after having created the Year col and coerced it to integer

#add Year column
EVA_forcas_meta$Year <- substr(x = EVA_forcas_meta$`Date of recording`, start = 7, stop = 10)

#count number of years
EVA_forcas_yrs <- as.data.frame(with(EVA_forcas_meta, table(Year)))

#check distribution of years in the dataset
ggplot(EVA_forcas_yrs, aes(x = Year, y = Freq)) +
  geom_col() +
  ylab('Count') + xlab(NULL) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 45, size = 14, vjust = 1, hjust = 1))

#coerce Year to integer
EVA_forcas_meta$Year <- as.integer(EVA_forcas_meta$Year)

#check NAs
sum(is.na(EVA_forcas_meta$Year)) #good, still 46
 
#exclude NAs
EVA_forcas_meta <- EVA_forcas_meta[!is.na(EVA_forcas_meta$Year), ]

#exclude plots sampled before 1985 because these will not be used for the analyses
sum(EVA_forcas_meta$Year >= 1985) #943/1101 - 85% are kept

#filter out observations sampled before 1985
EVA_forcas_meta <- EVA_forcas_meta[EVA_forcas_meta$Year >= 1985, ]

#subset EVA_full to only keep observations of plots needed for the analysis
EVA_forcas_species <- EVA_full[EVA_full$PlotObservationID %in% EVA_forcas_meta$PlotObservationID, ]

#check
length(unique(EVA_forcas_species$PlotObservationID)) #943

#remove EVA_full to save memory
rm(EVA_full)

#-----removed from the script (old procedure for matching data as in EVA and VegBank)

#crop area of metad_historical to retrieve EVA plots corresponding to Foreste Casentinesi
#bbox_metad_hist <- st_bbox(st_as_sf(metad_historical, coords = c('lon', 'lat'), crs = 4326))

#check if plot ID were maintained between EVA and data as originally contributed
#sum(is.na(match(metad_historical$ID, as.character(EVA_forcas_meta$`TV2 relevé number`)))) #1
#which(is.na(match(metad_historical$ID, as.character(EVA_forcas_meta$`TV2 relevé number`)))) #571
#sum(as.character(EVA_forcas_meta$`TV2 relevé number`) == metad_historical$ID[571]) #0

#check if year and coordinates are the same between TV2 relevé number and ID
#there should not be records without Year in metad_historical
#compare_datasets <- lapply(unique(metad_historical[!metad_historical$ID %in% "7608", 'ID']), function(id) {
#  plot_info <- metad_historical[metad_historical$ID == id, c('lon', 'lat', 'Year')]
#  plot_eva_coords <- EVA_forcas_meta[as.character(EVA_forcas_meta$`TV2 relevé number`) == id, c('lon', 'lat', 'Date of recording')]
#  #check coords
#  check_coords <- all.equal(unlist(plot_info[c(1, 2)]), unlist(plot_eva_coords[c(1, 2)]))
#  #check year
#  check_year <- all.equal(plot_info[[3]], as.integer(substr(plot_eva_coords[[3]], start = 7, stop = 10)))
#  res <- list(ck_crds = check_coords, ck_yr = check_year)
#  return(res)
#  })

#compare crds
#compare_crds <- lapply(compare_datasets, function(.) .[['ck_crds']])
#all(unlist(Filter(is.logical, compare_crds))) #T
#compare_crds[which(sapply(compare_crds, function(x) !is.logical(x)))]
#compare year
#compare_yr <- lapply(compare_datasets, function(.) .[['ck_yr']])
#all(unlist(Filter(is.logical, compare_yr))) #T

#subset EVA_forcas_meta to only keep plot in metad_historical (?)
#1097 plots in metad_historical -> 1096 plots in EVA_forcas_meta as 1 plot in metad_historical is missing in EVA_forcas_meta
#sum(duplicated(metad_historical$ID)) #0
#subset EVA_forcas_meta
#EVA_forcas_meta <- EVA_forcas_meta[EVA_forcas_meta$`TV2 relevé number` %in% metad_historical$ID, ]


#--------------------------------funny fact about assignment
#plot_eva_coords <- EVA_forcas_meta[as.character(EVA_forcas_meta$`TV2 relevé number`) == '7024', c('lon', 'lat', 'Date of recording')]
#plot_eva_coords[[3]] <- as.integer(substr(plot_eva_coords[[3]], start = 7, stop = 10))

#plot_eva_coords <- EVA_forcas_meta[as.character(EVA_forcas_meta$`TV2 relevé number`) == '7024', c('lon', 'lat', 'Date of recording')]
#plot_eva_coords[1, 3] <- as.integer(substr(plot_eva_coords[[3]], start = 7, stop = 10))

#plot_eva_coords <- EVA_forcas_meta[as.character(EVA_forcas_meta$`TV2 relevé number`) == '7024', c('lon', 'lat', 'Date of recording')]
#plot_eva_coords$`Date of recording` <- as.integer(substr(plot_eva_coords[[3]], start = 7, stop = 10))

#plot_eva_coords <- EVA_forcas_meta[as.character(EVA_forcas_meta$`TV2 relevé number`) == '7024', c('lon', 'lat', 'Date of recording')]
#plot_eva_coords[1, 3] <- as.numeric(substr(plot_eva_coords[[3]], start = 7, stop = 10))










