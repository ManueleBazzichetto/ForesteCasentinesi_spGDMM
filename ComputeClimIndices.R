#This script tests functionality of the ClimInd R package

#installed on Jan 13th 2025
install.packages('ClimInd')
library(ClimInd)

#notice that data_all (example dataset) has dates formatted with 2-digits years. This creates a problem for years before 1968. I discussed about the issue
#with the package maintainer (see mail to Fergus Reig on 24th May 2023), and he created a new data_all version to run examples. This new data version
#has years formatted as %Y (4-digits). He said that using this format should give the right outputs

#load example dataset - this data_all version has years formatted as %y
#data("data_all")

#load new data_all version from: https://gitlab.com/indecis-eu/indecis
load(file = 'data_all.RData')

#data_all is a list
data_all$tn[which(!is.na(data_all$tn))]

#each element of the list is the value of a climatic parameter for a given day - the Date format is rather weird

#example - computing bio10
ClimInd::bio10(data = data_all$tg) #results are per year!


#so the procedure to use data from easyclimate would be to:

#extract daily climate data

#compute Tmean

#create list with input for ClimInd - this includes creating a list of vectors that are named with dates formatted as in ClimInd

#--------indices to compute

#Basically, I can compute any index that is based on daily tmin, tmax, tmean and prec
#see https://web.archive.org/web/20190714191708/https://www.worldclim.org/bioclim

#bio10: Mean Temperature of Warmest Quarter
bio10(data = data_all$tg)
#bio11: Mean Temperature of Coldest Quarter
bio11(data = data_all$tg)
#bio13: Precipitation of Wettest Month
bio13(data = data_all$rr)
#bio14: Precipitation of Driest Month
bio14(data = data_all$rr)
#bio15: Precipitation Seasonality (Coefficient of Variation)
bio15(data = data_all$rr)
#bio16: Precipitation of Wettest Quarter
bio16(data = data_all$rr)
#bio17: Precipitation of Driest Quarter
bio17(data = data_all$rr)
#bio18: Precipitation of Warmest Quarter
bio18(pr = data_all$rr, taverage = data_all$tg)
#bio19: Precipitation of Coldest Quarter
bio19(pr = data_all$rr, taverage = data_all$tg)
#bio4: Temperature Seasonality (standard deviation *100)
bio4(data = data_all$tg)
#bio5: Max Temperature of Warmest Month
bio5(data = data_all$tg, tmax = data_all$tx)
#bio6: Min Temperature of Coldest Month
bio6(data = data_all$tg, tmin = data_all$tn)
#bio7: Temperature Annual Range (BIO5-BIO6)
bio7(data = data_all$tg, tmin = data_all$tn, tmax = data_all$tx)
#bio8: Mean Temperature of Wettest Quarter
bio8(pr = data_all$rr, taverage = data_all$tg)
#bio9: Mean Temperature of Driest Quarter
bio9(pr = data_all$rr, taverage = data_all$tg)
#gsl: Growing season length
gsl(data = data_all$tg)
#gsr: Growing season precipitation
gsr(data = data_all$rr)


compute_climind <- function(dtf) {
  #re-format date to comply with formatting in ClimInd
  old_date <- as.Date(dtf[["date"]], format = "%Y-%m-%d")
  new_date <- format(old_date, format = "%m/%d/%Y")
  #extract year climatic data
  year_climate <- format(old_date, format = '%Y')
  #generate input data for ClimInd
  input_climind <- lapply(dtf[c('Tmax', 'Tmin', 'Prcp', 'Tmean')], function(cl, nd = new_date) {
    climvar <- setNames(cl, nm = nd)
    return(climvar)
  })
  #compute average of bioclimatic variables over selected period
  res <- data.frame(bio4_avg = mean(bio4(data = input_climind$Tmean)),
                    bio5_avg = mean(bio5(data = input_climind$Tmean, tmax = input_climind$Tmax)),
                    bio6_avg = mean(bio6(data = input_climind$Tmean, tmin = input_climind$Tmin)),
                    bio7_avg = mean(bio7(data = input_climind$Tmean, tmin = input_climind$Tmin, tmax = input_climind$Tmax)),
                    bio8_avg = mean(bio8(pr = input_climind$Prcp, taverage = input_climind$Tmean)),
                    bio9_avg = mean(bio9(pr = input_climind$Prcp, taverage = input_climind$Tmean)),
                    bio10_avg = mean(bio10(data = input_climind$Tmean)),
                    bio11_avg = mean(bio11(data = input_climind$Tmean)),
                    bio13_avg = mean(bio13(data = input_climind$Prcp)),
                    bio14_avg = mean(bio14(data = input_climind$Prcp)),
                    bio15_avg = mean(bio15(data = input_climind$Prcp)),
                    bio16_avg = mean(bio16(data = input_climind$Prcp)),
                    bio17_avg = mean(bio17(data = input_climind$Prcp)),
                    bio18_avg = mean(bio18(pr = input_climind$Prcp, taverage = input_climind$Tmean)),
                    bio19_avg = mean(bio19(pr = input_climind$Prcp, taverage = input_climind$Tmean)),
                    gsl_avg = mean(gsl(data = input_climind$Tmean)),
                    gsr_avg = mean(gsr(data = input_climind$Prcp)),
                    tmean_avg = mean(tapply(X = input_climind$Tmean, INDEX = year_climate, mean)),
                    tmax_avg = mean(tapply(X = input_climind$Tmax, INDEX = year_climate, mean)),
                    tmin_avg = mean(tapply(X = input_climind$Tmin, INDEX = year_climate, mean)),
                    prcp_avg = mean(tapply(X = input_climind$Prcp, INDEX = year_climate, sum)))
  return(res)
  }








