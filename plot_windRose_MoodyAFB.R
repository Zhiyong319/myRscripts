# set the environment
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(openair)

# read raw data
rawdat2015 <-read.csv('74781013857_2015.csv')
rawdat2016 <-read.csv('74781013857_2016.csv')
rawdat2017 <-read.csv('74781013857_2017.csv')
rawdat2018 <-read.csv('74781013857_2018.csv')
rawdat2019 <-read.csv('74781013857_2019.csv')
cols <- c("DATE","WND","REPORT_TYPE" )
rawdat <- rbind(rawdat2015[,cols],rawdat2016[,cols],rawdat2017[,cols],rawdat2018[,cols],rawdat2019[,cols])

rawdat$datetime <- as.POSIXct(paste(substr(rawdat$DATE,1,10),substr(rawdat$DATE,12,19)), origin="1970-01-01", tz="UTC")
rawdat$year  <- as.numeric(substr(rawdat$DATE,1,4))
rawdat$month <- as.numeric(substr(rawdat$DATE,6,7))
rawdat$day   <- as.numeric(substr(rawdat$DATE,9,10))
rawdat$hour  <- as.numeric(substr(rawdat$DATE,12,13))
rawdat$wdir  <- as.numeric(substr(rawdat$WND,1,3))
rawdat$wdir_qc  <- as.numeric(substr(rawdat$WND,5,5))
rawdat$wspd  <- as.numeric(substr(rawdat$WND,9,12))/10  # m/s
rawdat$wspd_qc  <- as.numeric(substr(rawdat$WND,14,14))

# filter data by qc flag
rawdat <- subset(rawdat, REPORT_TYPE=="FM-15" & (wdir_qc==1 | wdir_qc==5) & (wspd_qc==1 | wspd_qc==5))

# set season flags
season_bins <- seq(as.POSIXct("2014-12-01 00:00:00", origin="1970-01-01", tz="UTC"),
                   as.POSIXct("2019-12-01 00:00:00", origin="1970-01-01", tz="UTC"), by = "quarter")
rawdat$season_flag <- as.numeric(cut(rawdat$datetime, season_bins))
seasonNames <- c('(01) Spring 2015','(02) Summer 2015','(03) Fall 2015','(04) Winter 2015',
                 '(05) Spring 2016','(06) Summer 2016','(07) Fall 2016','(08) Winter 2016',
                 '(09) Spring 2017','(10) Summer 2017','(11) Fall 2017','(12) Winter 2017',
                 '(13) Spring 2018','(14) Summer 2018','(15) Fall 2018','(16) Winter 2018',
                 '(17) Spring 2019','(18) Summer 2019','(19) Fall 2019','(20) Winter 2019')
rawdat$seasonName <- seasonNames[rawdat$season_flag]

seasonofMonth <- c('(01) Spring', '(01) Spring', 
                   '(02) Summer', '(02) Summer','(02) Summer',
                   '(03) Fall', '(03) Fall','(03) Fall',
                   '(04) Winter', '(04) Winter', '(04) Winter',
                   '(01) Spring')
rawdat$seasonName2 <- seasonofMonth[rawdat$month]

# windrose for each season of each year
png('MoodyAFB_seasonal_windroses.png',width=7.2, height=7.2*1.2, units="in",res=1200)
windRose(rawdat, ws = "wspd", wd = "wdir", type="seasonName")
dev.off()

# windrose for each season of all years
png('MoodyAFB_seasonal_windroses_allyears.png',width=7.2, height=7.2, units="in",res=1200)
windRose(rawdat, ws = "wspd", wd = "wdir", type="seasonName2")
dev.off()

# data points for each season of each year
num_of_points <- table(rawdat$seasonName)
write.csv(num_of_points, file='num_of_points.')
