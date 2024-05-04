library(dataRetrieval)
library(dplyr)
library(sjmisc)
library(plyr)
#library(leaflet)
#library(mapview)

setwd('D:/17_TOVA/PCR-GLOBWB2-post-process')

# Code '00060' corresponds to discharge
pCode <- "00060"

# State that I want to extract discharge

State_US <- c('AL', 'AK', 'AZ', 'AR', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'HI', 'ID', 'IL', 'IN',
 'IA', 'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 'NV', 'NH', 'NJ', 
 'NM', 'NY', 'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 
 'WA', 'WV', 'WI', 'WY')

#State_US <- c('TX', 'UT', 'VT', 'VA', 'WA', 'WV', 'WI', 'WY')

SurfaceArea_all_gages <- list()
Lat_all_gages <- list()
Lon_all_gages <- list()


metainfo <- data.frame(matrix(nrow = 0, ncol = length(c('iD', 'area', 'lat', 'lon'))))
colnames(metainfo) <- c('iD', 'area_km2', 'lat', 'lon')

# Loop over the states
for (state in State_US){

  # Read all stations in 'state' 
  flow_df <- readNWISdata(stateCd=state, parameterCd=pCode,
                          service="site", seriesCatalogOutput=TRUE)

  # Selection of the stations that matches with the criteria below
  Flow <- filter(flow_df, parm_cd %in% pCode) %>% filter(count_nu > 1000)

  # Find and remove replicate of 'site_no' (it seems that duplicates are not present in the 
  # SurfaceArea list)
  id_duplicated <- which(duplicated(Flow[,'site_no']))
  if (is_empty(id_duplicated) == FALSE){
    Flow <- Flow[-c(id_duplicated),]
  }

  # Read the area of all currently selected gages (km2)
  SurfaceArea <- readNWISsite(Flow[,'site_no'])[['drain_area_va']] * 2.590
  
  # Remove gages that have missing basin's area
  id_NA <- which(is.na(SurfaceArea), arr.ind = TRUE)
  if (is_empty(id_NA) == FALSE){
    Flow <- Flow[-c(id_NA),]
    SurfaceArea <- SurfaceArea[-c(id_NA)]
  }

  # Remove gages that are smaller than 2,000 knm2
  id_small <- which(SurfaceArea<2000, arr.ind = TRUE)
  if (is_empty(id_small) == FALSE){
    Flow <- Flow[-c(id_small),]
    SurfaceArea <- SurfaceArea[-c(id_small)]
  }

  # Remove gages with startDate after 2015-01-01 and endDate before 1950-01-01
  id_too_new <- which(Flow[,'begin_date'] > '2015-01-01')
  if (is_empty(id_too_new) == FALSE){
    Flow <- Flow[-c(id_too_new),]
    SurfaceArea <- SurfaceArea[-c(id_too_new)]
  }

  id_too_old <- which(Flow[,'end_date'] < '1951-01-01')
  if (is_empty(id_too_old) == FALSE){
    Flow <- Flow[-c(id_too_old),]
    SurfaceArea <- SurfaceArea[-c(id_too_old)]
  }

  if (is_empty(Flow[,'site_no']) == FALSE){
    # Remove gages that have less than 20 years of data
    flow_data <- readNWISdv(siteNumbers = Flow[,'site_no'], parameterCd=pCode,
                            startDate = '1951-01-01', endDate = '2015-12-31')

    num_of_years <- flow_data %>%
      group_by(site_no) %>%
      summarize(duration_years = as.numeric(difftime(max(Date), min(Date), units = "days") / 365.25))

    
    # Get rid of the gages with less than 20 years of data between 1951 and 2015
    id_sites_without_enough_data <- which(num_of_years['duration_years']<20)
    if (is_empty(id_sites_without_enough_data) == FALSE){
      Flow <- Flow[-c(id_sites_without_enough_data),]
      SurfaceArea <- SurfaceArea[-c(id_sites_without_enough_data)]
    }

    # Get rid of the gages not in 'num_of_years' array'
    difference <- mapply(setdiff,Flow[,'site_no'], num_of_years['site_no'])
    non_empty_char_cols <- which(sapply(difference, function(x) is.character(x) && length(x) > 0))
    if (is_empty(non_empty_char_cols) == FALSE){
      Flow <- Flow[-c(non_empty_char_cols),]
      SurfaceArea <- SurfaceArea[-c(non_empty_char_cols)]
    }
    sdf
    # Loop over the selected gages and pick the daily streamflow time series
    id_no_unique_flow_label <- list()
    for (i in seq(dim(Flow)[1])) {

      begin_date = Flow[i,'begin_date']
      end_date   = Flow[i,'end_date']
        
      Q <- readNWISdata(sites=Flow[i,'site_no'], service="dv", 
                      parameterCd=pCode, 
                      startDate=begin_date,endDate=end_date)
      
      # Need to check if the time series is not split into multiple period
      # For the sake of simplicity, we disregard gages that have there historical records
      # split into multiple columns.
      if ('X_00060_00003' %in% names(Q)){
        Q <- subset(Q, select = c('dateTime', 'X_00060_00003'))
        Q <- Q %>% rename('Dates' = 'dateTime')
        Q <- Q %>% rename('Flow_cfs'='X_00060_00003')

      
        write.csv(Q, file = paste('data/USGS/USGS_',Flow[i,'site_no'],'.csv', sep=''), row.names=FALSE)
      } else {
        id_no_unique_flow_label <- rbind(id_no_unique_flow_label, i)
      }
      
        
    }

    # Update the list of basin and surface area list.
    if (is_empty(id_no_unique_flow_label) == FALSE){
      Flow <- Flow[-c(unlist(id_no_unique_flow_label)),]
      SurfaceArea <- SurfaceArea[-c(unlist(id_no_unique_flow_label))]
    }

    # Get latitude and longitude of the selected gages
    latitude <- readNWISsite(Flow[,'site_no'])[['dec_lat_va']]
    longitude <- readNWISsite(Flow[,'site_no'])[['dec_long_va']]

    # Attributes for the current state
    additional_basins = data.frame(
      iD = c(Flow[,'site_no']),
      area_km2 = c(SurfaceArea),
      lat = c(latitude),
      lon = c(longitude)
    )
    
    # Append the metadata/attribute to the other states
    metainfo <- rbind(metainfo, additional_basins)
  }

}

non_numerics<-adply(1:ncol(metainfo),1,function(x)print(is.numeric(metainfo[,x])))
quote_val<-as.numeric(array(non_numerics[which(!non_numerics$V1),1]))

write.csv(metainfo, file = paste('data/USGS/MetaData_USGS.csv', sep=''), 
          row.names = FALSE, quote = quote_val)

