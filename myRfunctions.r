#
# Collect all the useful functions
#
# Author: Zhiyong Wu (zhiyong319@gmail.com)

#### Collection of data ####
# time
date_start <- c(  1, 1, 32, 60,  91, 121, 152, 182, 213, 244, 274, 305, 335)
date_end   <- c(365,31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
period <- sprintf("%02d", 0:12) 
period2 <-c('Annual','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

# US states
# state.abb  # state abbreviations 
# state.name # state full names 

#### Collection of functions ####
state.name2abb <- function(state_vector) {
  for (i in 1:length(state.name)) {
    state_vector[state_vector==state.name[i]] <- state.abb[i]
  }
  return(state_vector)
}

MinMaxScaling <- function(x) {
  minvalue <- min(x, na.rm = T)
  maxvalue <- max(x, na.rm = T)
  # maxvalue <- quantile(x, 0.995, na.rm = T)
  return((x-minvalue)/(maxvalue-minvalue))
}

CoeffOfDivergence <- function(x,y) {
  ii <- !is.na(x) & !is.na(y)
  xx <- x[ii]
  yy <- y[ii]
  
  return(sqrt(mean(((xx-yy)/(xx+yy))^2)))
}

# Create a function to move files from one directory to another using file.rename
# In the process check if destination directory exists and if not, create one
move.file <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

match_rows <- function(row_match,row_raw) {
  indicies <- vector(length=length(row_raw))
  for (i in 1:length(row_raw)) {
    indicies[i] <- which(row_match==row_raw[i])
  }
                     
  return(indicies)
}

### wide table to narrow
dataframe_reshape <- function(df) {
  # the 1st column of input df is kept
  # original df: time, column1, column2, ...
  # reshaped df: time, ColumnName, value
  
  num_column_reshaped <- ncol(df)-1
  ColumnName <- colnames(df)[2:ncol(df)]
  
  df2 <- data.frame(matrix(nrow = nrow(df)*num_column_reshaped, ncol=3))
  colnames(df2) <- c('time','ColumnName','value')
  
  df2$time <- rep(df$time, each=num_column_reshaped)
  df2$ColumnName <- rep(ColumnName, nrow(df))
  df2$value <- c(t(as.matrix(df[,2:ncol(df)])))
    
  # for (i in 1:nrow(df)) {
  #   df2[((i-1)*num_column_reshaped+1):i*num_column_reshaped,3] <- df[i,2:ncol(df)]
  # }
  
  return(df2)
}

### calculate daily averages using hourly data
Hourly2Daily <- function(hourlyData) {
  # hourlyData [time, site1, site2, ...]
  # the hourly time starts from 00 and ends with 23 without missing rows
  
  Sys.setenv(TZ='GMT')
  date_start <- as.Date(min(hourlyData$time),format = "%Y-%m-%d")  # start date
  date_end <- as.Date(max(hourlyData$time),format = "%Y-%m-%d") # end date
  seq_dates <- seq(as.Date(date_start), as.Date(date_end), by="days")
  
  dailyData <- data.frame(matrix(nrow = length(seq_dates), ncol = ncol(hourlyData)))
  colnames(dailyData) <- colnames(hourlyData)
  dailyData$time <- seq_dates
  
  for (iday in 1:length(seq_dates)) {
    dailyData[iday,2:ncol(dailyData)] <- colMeans(hourlyData[((iday-1)*24+1):(iday*24), 2:ncol(hourlyData)],na.rm = T)
  }
  
  return(dailyData)
  
}

##########################################################################################################
##########################################################################################################
#####-------------------   START OF FUNCTION: LAMB_LATLON_TO_IJ         ------------------------------####
# Cacluation of i and j index of a specified latitude and longitude for a given lambert conformal
# projection.
# Input:
#        reflat  -- Latitude of a reference point of the grid. Typically lower left corner of the domain
#        reflon  -- Longitude of a reference point of the grid. Typically lower left corner of the domain
#        iref    -- I index of the reference point.
#        jref    -- J index of the reference point.
#        truelat1-- first true latitude of the projection
#        truelat2-- second true latitude of the projection
#        stdlon  -- center longitude of the LC projection
#        delx    -- grid spacing in meters
#        grdlat  -- Latitude of point of interest (usually an observation location)
#        grdlon  -- Longitude of point of interest (usually an observation location)
#        radius  -- Radius of the earth in meters. Default is current value for WRF and MPAS
#
# Output: List with i and j index values in fractional form

lamb_latlon_to_ij <- function(reflat, reflon, iref, jref, truelat1, truelat2,
                              stdlon, delx, grdlat, grdlon, radius= 6370000.0) {
  
  
  pi     <- 4.0*atan(1.0)
  pi2    <- pi/2.0
  pi4    <- pi/4.0
  d2r    <- pi/180.0
  r2d    <- 180.0/pi
  omega4 <- 4.0*pi/86400.0
  
  if(truelat1 == truelat2) {
    gcon <- sin( abs(truelat1) * d2r)
  } else {
    gcon <- (log(sin((90.0- abs(truelat1))* d2r))-log(sin((90.0- abs(truelat2))* d2r)))/
      (log(tan((90.0- abs(truelat1))*0.5* d2r))-log(tan((90.0- abs(truelat2))*0.5* d2r)))
  }
  
  ogcon  <- 1.0/gcon
  ahem   <- abs( truelat1 / truelat1 )
  deg    <- (90.0- abs(truelat1)) * d2r
  cn1    <- sin(deg)
  cn2    <- radius * cn1 * ogcon
  deg    <- deg * 0.5
  cn3    <- tan(deg)
  deg    <- (90.0- abs(reflat)) * 0.5 * d2r
  cn4    <- tan(deg)
  rih    <- cn2 *(cn4 / cn3)**gcon
  deg    <- (reflon-stdlon) * d2r * gcon
  
  xih    <- rih * sin(deg) - delx/2
  yih    <- (-1 * rih * cos(deg) * ahem) - delx/2
  
  deg    <- (90.0 - (grdlat * ahem) ) * 0.5 * d2r
  cn4    <- tan(deg)
  
  rrih   <- cn2 *  (cn4/cn3)**gcon
  check  <- 180.0-stdlon
  alnfix <- stdlon + check
  alon   <- grdlon + check
  
  if (alon<0.0)  { alon <- alon+360.0 }
  if (alon>360.0){ alon <- alon-360.0 }
  
  deg    <- (alon-alnfix) * gcon * d2r
  XI     <- rrih * sin(deg)
  XJ     <- -1* rrih * cos(deg) * ahem
  
  grdi   <- iref+(XI-xih)/(delx)
  grdj   <- jref+(XJ-yih)/(delx)
  
  return(list(i=grdi, j=grdj))
  
}
#####--------------------------	  END OF FUNCTION: LAMB_LATLON_TO_IJ  --------------------------------####
##########################################################################################################

##########################################################################################################
#####-------------------   START OF FUNCTION: POLARS_LATLON_TO_IJ       ------------------------------####
# Cacluation of i and j index of a specified latitude and longitude for a given lambert conformal
# projection. Note: Calculations derived from Obsgrid Code
# Input:
#        reflat  -- Latitude of a reference point of the grid. Typically lower left corner of the domain
#        reflon  -- Longitude of a reference point of the grid. Typically lower left corner of the domain
#        truelat1-- first true latitude of the projection
#        stdlon  -- center longitude of the LC projection
#        delx    -- grid spacing in meters
#        lato    -- Latitude of point of interest (usually an observation location)
#        lono    -- Longitude of point of interest (usually an observation location)
#        radius  -- Radius of the earth in meters. Default is current value for WRF and MPAS
#
# Output: List with i and j index values in fractional form

polars_latlon_to_ij <- function(lat1, lon1, delx, truelat1, stdlon, 
                                lato, lono, radius= 6370000.0) {
  
  rad_per_deg <- pi/180
  deg_per_rad <- 180/pi
  
  reflon      <- stdlon + 90.
  
  hemi   <- 1
  if(truelat1 < 0) {
    hemi <- -1
  }
  
  rebydx    <- radius / delx
  scale_top <- 1. + hemi * sin(truelat1 * rad_per_deg)
  
  ala1      <- lat1 * rad_per_deg
  rsw       <- rebydx*cos(ala1)*scale_top/(1.+hemi*sin(ala1))
  
  alo1      <- (lon1 - reflon) * rad_per_deg
  polei     <- 1. - rsw * cos(alo1)
  polej     <- 1. - hemi * rsw * sin(alo1)
  
  
  ala       <- lato * rad_per_deg
  rm        <- rebydx * cos(ala) * scale_top/(1. + hemi *sin(ala))
  alo       <- (lono - reflon) * rad_per_deg
  grdi      <- polei + rm * cos(alo)
  grdj      <- polej + hemi * rm * sin(alo)
  
  
  return(list(i=round(grdi), j=round(grdj)))
  
}
#####--------------------------	  END OF FUNCTION: POLARS_LATLON_TO_IJ  ------------------------------####
##########################################################################################################

##########################################################################################################
##########################################################################################################
#####-------------------   START OF FUNCTION: MERCAT_LATLON_TO_IJ       ------------------------------####
# Cacluation of i and j index of a specified latitude and longitude for a given mercator projection.
# Note: Calculations derived from Obsgrid Code
# Input:
#        reflat  -- Latitude of a reference point of the grid. Typically lower left corner of the domain
#        reflon  -- Longitude of a reference point of the grid. Typically lower left corner of the domain
#        lat     -- Latitude array (2D) of domain gridpoints
#        lon     -- Longitude array (2D) of domain gridpoints
#        dx      -- Grid spacing in meters
#        stdlt1  -- first true latitude of the projection
#        lato    -- Latitude of point of interest (usually an observation location)
#        lono    -- Longitude of point of interest (usually an observation location)
#        radius  -- Earth radius in model (m) default is set to WRF.
#
# Output: List with i and j index values in fractional form

mercat_latlon_to_ij <- function(reflat, reflon, lat, lon, dx, stdlt1, lato, lono, radius= 6370000.0) {
  
  # Calcuation used in Obsgrid for Mercator. Saved here for possible use in the future.
  rad_per_deg <- pi/180
  deg_per_rad <- 180/pi
  
  lon0 <- reflon
  if(lon0 < 0) {
    lon0 <- 360+lon0
  } 
  
  lon360 <- lono
  if(lon360 < 0) {
    lon360 <- 360+lon360
  } 
  
  clain  <- cos(rad_per_deg*stdlt1)
  dlon   <- dx / (radius * clain)
  
  rsw <- 0.
  if(reflat != 0) {
    rsw <- (log(tan(0.5*((reflat+90.)*rad_per_deg))))/dlon
  }
  
  
  deltalon <- lon360 - lon0
  deltalat <- lato - reflat
  
  grdi <- 1. + (deltalon/(dlon*deg_per_rad))
  grdj <- 1. + (log(tan(0.5*((lato + 90.) * rad_per_deg)))) / dlon - rsw
  
  # More direct calculation that finds closest grid point using min distance
  # between site lat-lon and all grid point lat-lons.
  d    <- sqrt( (lato-lat)^2 + (lono-lon)^2 )
  ind  <- which(d == min(d), arr.ind =T)
  grdi <-ind[1]
  grdj <-ind[2]
  
  return(list(i=round(grdi), j=round(grdj)))
  
}
#####--------------------------	  END OF FUNCTION: MERCAT_LATLON_TO_IJ  --------------------------------####
##########################################################################################################

##########################################################################################################
#####-------------------   START OF FUNCTION: LATLON_LATLON_TO_IJ         ------------------------------####
# Cacluation of i and j index of a specified latitude and longitude for a given lat-lon projection.
# 
# Input:
#        lat     -- Latitude array (2D) of domain gridpoints
#        lon     -- Longitude array (2D) of domain gridpoints
#        lato    -- Latitude of point of interest (usually an observation location)
#        lono    -- Longitude of point of interest (usually an observation location)
#
# Output: List with i and j index values in fractional form

latlon_latlon_to_ij <- function(lat, lon, lato, lono) {
  
  # More direct calculation that finds closest grid point using min distance
  # between site lat-lon and all grid point lat-lons.
  d    <- sqrt( (lato-lat)^2 + (lono-lon)^2 )
  ind  <- which(d == min(d), arr.ind =T)
  grdi <-ind[1]
  grdj <-ind[2]
  
  return(list(i=grdi, j=grdj))
  
}
#####--------------------------	  END OF FUNCTION: LATLON_LATLON_TO_IJ  --------------------------------####

##########################################################################################################
#####--------------------------   START OF FUNCTION: CONE         ------------------------------------####
#  Cone factor function
# Input:
#        true1  -- first true latitude of the projection
#        true2  -- second true latitude of the projection
#         hemi  -- Hemisphere of the projection (N=1, S=-1)
#
# Output: Cone factor

cone <-function(true1,true2,hemi=1) {
  
  true1  <-33
  true2  <-45
  hemi   <- 1   
  
  pi180  <-pi/180
  
  cone   <- log10(cos(true1 * pi180)) - log10(cos(true2 * pi180))
  cone   <- cone / ( log10( tan((45.0 - hemi*true1/2.0) * pi180)) 
                     -  log10( tan((45.0 - hemi*true2/2.0) * pi180)) )
  
  # if met_tru1 and met_tru2 are the same
  if(true1 == true2) {
    cone <- hemi * sin(true1*pi180)
  }
  return(cone)
  
}
#####--------------------------	  END OF FUNCTION: CONE           ------------------------------------####
##########################################################################################################

# revise from mcip_surface() in MET_model.read.R
##########################################################################################################
#####--------------------------   START OF FUNCTION: MCIP_projection  ------------------------------------####
#  Open MCIP output and extract grid information and surface met data for comparision to MADIS obs
# Input:
#       file   -- Model output file name. Full path and file name
#
# Output: Multi-level list of model projection data.
#
#  projection <-list(mproj=mproj, lat=lat, lon=lon, lat1=lat1, lon1=lon1, nx=nx, ny=ny, dx=dx,
#                    truelat1=truelat1, truelat2=truelat2, standlon=standlon, conef=conef)

mcip_projection <-function(gfile) {
  
  # For MCIP, users must put a generic GRIDCRO2D file in the same directory
  # as METCRO2D files. Lines below piece together directory location from "file"
  # and set grid file name "gfile" to generic GRIDCRO2D.
  # If user does not provide a GRIDCRO2D file in correct location, error message
  # informs user of the mistep.
  
  #gfile    <-paste(paste0(dirs[1:(ndirs-1)], collapse="/"),"/GRIDCRO2D",sep="")
  if(file.exists(gfile)){
    writeLines(paste("Reading supplied GRIDCRO2D for grid information:",gfile))
    fg      <- nc_open(gfile)
    head    <- ncatt_get(fg, varid=0, attname="EXEC_ID" )$value
    tproj   <- ncatt_get(fg, varid=0, attname="GDTYP" )$value
  }
  else {
    writeLines(paste("Terminating AMET."))
    writeLines(paste("User must supply a single GRIDCRO2D as follows:",gfile))
    quit(save="no")
  }   
  
  # Note MCIP grid different from WRF grid numbers
  # MCIP 1-lat/lon, 2-Lambert, 3-Mercator, 6- Polar Stereo
  # WRF  1-Lamb     2-Polar    3-Merc      4- Lat-lon
  if(tproj != 1 & tproj !=2 & tproj !=3 & tproj !=6) { 
    writeLines(paste("Grid projection not compatible. Has to be 1, 2, 3 or 6."))
    writeLines(paste("Terminating AMET. Projection has to be lat-lon, lambert, Merc or PolarS."))
    quit(save="no")
  }
  # Mapping MCIP Proj numbering to WRF 
  if(tproj == 1) { mproj <- 4 }
  if(tproj == 2) { mproj <- 1 }
  if(tproj == 3) { mproj <- 3 }
  if(tproj == 6) { mproj <- 2 }
  
  if(mproj == 1){
    lat1    <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(1,1,1,1))
    lon1    <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(1,1,1,1))
    lat     <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(-1,-1,1,1))
    lon     <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(-1,-1,1,1))
    nx      <- ncatt_get(fg, varid=0, attname="NCOLS" )$value
    ny      <- ncatt_get(fg, varid=0, attname="NROWS" )$value
    dx      <- ncatt_get(fg, varid=0, attname="XCELL" )$value
    truelat1<- ncatt_get(fg, varid=0, attname="P_ALP" )$value
    truelat2<- ncatt_get(fg, varid=0, attname="P_BET" )$value
    standlon<- ncatt_get(fg, varid=0, attname="P_GAM" )$value
  }
  if(mproj == 2){
    lat1    <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(1,1,1,1))
    lon1    <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(1,1,1,1))
    lat     <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(-1,-1,1,1))
    lon     <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(-1,-1,1,1))
    nx      <- ncatt_get(fg, varid=0, attname="NCOLS" )$value
    ny      <- ncatt_get(fg, varid=0, attname="NROWS" )$value
    dx      <- ncatt_get(fg, varid=0, attname="XCELL" )$value
    truelat1<- ncatt_get(fg, varid=0, attname="P_ALP" )$value
    truelat2<- ncatt_get(fg, varid=0, attname="P_BET" )$value
    standlon<- ncatt_get(fg, varid=0, attname="P_GAM" )$value
  }
  if(mproj == 3){
    lat1    <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(1,1,1,1))
    lon1    <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(1,1,1,1))
    lat     <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(-1,-1,1,1))
    lon     <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(-1,-1,1,1))
    nx      <- ncatt_get(fg, varid=0, attname="NCOLS" )$value
    ny      <- ncatt_get(fg, varid=0, attname="NROWS" )$value
    dx      <- ncatt_get(fg, varid=0, attname="XCELL" )$value
    truelat1<- ncatt_get(fg, varid=0, attname="P_ALP" )$value
    truelat2<- ncatt_get(fg, varid=0, attname="P_BET" )$value
    standlon<- ncatt_get(fg, varid=0, attname="P_GAM" )$value
  }
  
  if(mproj == 4){
    lat1    <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(1,1,1,1))
    lon1    <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(1,1,1,1))
    lat     <- ncvar_get(fg, varid="LAT",start=c(1,1,1,1),count=c(-1,-1,1,1))
    lon     <- ncvar_get(fg, varid="LON",start=c(1,1,1,1),count=c(-1,-1,1,1))
    nx      <- ncatt_get(fg, varid=0, attname="NCOLS" )$value
    ny      <- ncatt_get(fg, varid=0, attname="NROWS" )$value
    dx      <- ncatt_get(fg, varid=0, attname="XCELL" )$value
    truelat1<- ncatt_get(fg, varid=0, attname="P_ALP" )$value
    truelat2<- ncatt_get(fg, varid=0, attname="P_BET" )$value
    standlon<- ncatt_get(fg, varid=0, attname="P_GAM" )$value
  }
  nc_close(fg)
  
  # MCIP wind dir is adjusted to true north from grid relative north in WRF
  # A readjustment back to WRF is required before U and V are caculated.
  # This adjustment is added to MCIP WD10 in time loop below.
  conef    <- cone(truelat1,truelat2)
  wdadjust <- (standlon-lon)*conef
  
  projection <-list(mproj=mproj, lat=lat, lon=lon, lat1=lat1, lon1=lon1, nx=nx, ny=ny, dx=dx,
                    truelat1=truelat1, truelat2=truelat2, standlon=standlon, conef=conef)
  
  return(projection)
  
}
#####--------------------------	  END OF FUNCTION: MCIP_projection     -----------------------------------####
##########################################################################################################
