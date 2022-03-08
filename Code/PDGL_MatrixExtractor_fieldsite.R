################################################################################
# Extract relevant data for all identified lakes for a selected area
# 
# PDGL_MatrixExtractor_fieldsite.R
#
# ReadMe: 
#
# Input:
#   shp file of lakes
#   SRTM for the region
#   RGI for the region
# 
#
# Created:          2022/08/15
# Latest Revision:  2022/05/05
#
# Jakob F Steiner| ICIMOD | jakob.steiner@icimod.org | x-hydrolab.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# install necessary packages if not available yet via install.packages()
library(pacman)
p_load(rgdal,rgeos,maptools,raster,rasterVis)
library(RColorBrewer)
library(sp)
library(geosphere)
library(secr)

##########################
# SPECIFY FILENAMES AND DESIRED PROJECTION
##########################
# Define Projection
projec_utm <- '+proj=utm +datum=WGS84'
projec_deg <- '+proj=longlat +datum=WGS84 +no_defs'

projec <- projec_utm

##########################
# Load RGI
##########################
path_RGI <- 'D:\\Work\\GeospatialData\\RGI60'                                     # Folder for RGI glacier outlines
RGI15_filename <- '15_rgi60_SouthAsiaEast.shp'

ogrInfo(path_RGI&'\\'&RGI15_filename)
RGI60_15<-readOGR(dsn=path_RGI&'\\'&RGI15_filename)
#RGI60_15 <- spTransform(RGI60_15, projec)

##########################
# Load DEM
##########################
path_DEM <- 'D:\\Work\\GeospatialData\\HMA\\SRTM\\VoidFilled'                      # Folder for RGI glacier outlines
DEM_filename <- 'SRTM_Corrected_HMA.tif'                                          # RGI filename

regionalDEM <- raster(path_DEM&'\\'&DEM_filename)

##########################
# Load Seisimic Data
##########################
path_Seismo <- 'D:\\Work\\ICIMODProjects\\TBWG\\LakeStandardization\\SeismicData' 
Seismo_filename <- 'gshap_globe1.tif'                                          # seismic acceleration filename [m s-2]

regionalSeismic <- raster(path_Seismo&'\\'&Seismo_filename)
regionalSeismic <- crop(regionalSeismic,extent(regionalDEM))

##########################
# Load Precipitation Intensity Data
##########################
path_Precip <- 'D:\\Work\\ICIMODProjects\\TBWG\\LakeStandardization\\PPDIST\\' 
Precip_filename <- 'daily_intensity_point_scale.tif'                           # precipitation intensity with 10yr return period [mm d-1]

regionalPrecip <- raster(path_Precip&'\\'&Precip_filename)
regionalPrecip <- crop(regionalPrecip,extent(regionalDEM))

##########################
# Load Permafrost Data
##########################
path_pf <- 'D:\\Work\\GeospatialData\\HMA\\Permafrost' 
pf_filename <- 'PZI_HMA.tif'                                          # PZI filename

regionalPF <- raster(path_pf&'\\'&pf_filename)
regionalPF <- crop(regionalPF,extent(regionalDEM))
##########################
# Load Lake locations
##########################

path_lakes <- 'D:\\Work\\Code\\GLOFStandardization\\Data\\TrainingExamples\\DudhKosi'

lake_filename_shp <- 'Lakes.shp'                  # lakes filename
ogrInfo(path_lakes&'\\'&lake_filename_shp)
lakes_loc<-readOGR(dsn=path_lakes&'\\'&lake_filename_shp)
lakes_loc_deg <- spTransform(lakes_loc,projec_deg)
lakes_loc <- spTransform(lakes_loc,projec)


for(i in 1:length(lakes_loc$GLIMS_ID)){

#########################
# Evaluate individual lakes
#########################

NumLakes <- length(lakes_loc)

# Crop Products

#if(lakes_loc_deg$Longitude[i]==0){
  lakes_loc_deg$Longitude[i]<-coordinates(lakes_loc_deg[i,])[1]
  lakes_loc_deg$Latitude[i]<-coordinates(lakes_loc_deg[i,])[2]

extentLLoc <- c(lakes_loc_deg$Longitude[i] - 0.02, lakes_loc_deg$Longitude[i] + 0.02, lakes_loc_deg$Latitude[i] - 0.02, lakes_loc_deg$Latitude[i] + 0.02)
LDEM <- crop(regionalDEM,extentLLoc)
LDEM_deg <- LDEM

LSeismic <- crop(regionalSeismic,extentLLoc + c(-0.1,0.1,-0.1,0.1))
LSeismic_deg <- LSeismic

LPrecip <- crop(regionalPrecip,extentLLoc + c(-0.1,0.1,-0.1,0.1))
LPrecip_deg <- LPrecip

LPZI <- crop(regionalPF,extentLLoc)
LPZI_deg <- LPZI

# Restrict the datasets to the domain
sub_RGI15 <- subset(RGI60_15, CenLon >= extent(LDEM_deg)[1] & CenLon <= extent(LDEM_deg)[2] & CenLat >= extent(LDEM_deg)[3] & CenLat <= extent(LDEM_deg)[4])

if(!is.null(tryCatch(mask(LDEM_deg,sub_RGI15), error = function(e) NULL))){
  RGI_selected <- sub_RGI15
}else{
  RGI_selected <- NULL
}

LDEM <- projectRaster(LDEM,crs=projec)
terrain_slope <- terrain(LDEM,'slope',unit="degrees")
terrain_slope_deg <- terrain(LDEM_deg,'slope',unit="degrees")

if(!is.null(RGI_selected)){
# TRansform everything to UTM
RGI_selected_deg <- RGI_selected
RGI_selected <- spTransform(RGI_selected,projec)

glac_locs <- mask(LDEM,RGI_selected)


#if(is.na(cellStats(mask(LDEM_deg,lakes_loc_deg[i,]),'mean'))){
  dist.mat <- geosphere::dist2Line(p = rasterToPoints(mask(disaggregate(LDEM_deg,fact=2),lakes_loc_deg[i,]))[,1:2], line = RGI_selected_deg)
  pts.wit.dist <- cbind(rasterToPoints(mask(disaggregate(LDEM_deg,fact=2),lakes_loc_deg[i,]))[,1:2], dist.mat)
#}else{
#dist.mat <- geosphere::dist2Line(p = rasterToPoints(mask(LDEM_deg,lakes_loc_deg[i,]))[,1:2], line = RGI_selected_deg)
#pts.wit.dist <- cbind(rasterToPoints(mask(LDEM_deg,lakes_loc_deg[i,]))[,1:2], dist.mat)}
# ID and slope of the closest glacier
closeGLID <- pts.wit.dist[1,6]
#Extent close to the nearest glacier tip
tiploc <- pts.wit.dist[which.min(pts.wit.dist[,3]),4:5]
extent_tiploc <- extent(tiploc[1]-0.005,tiploc[1]+0.005,tiploc[2]-0.005,tiploc[2]+0.005)
lakes_loc$slopelowertongue[i] <- mean(extract(mask(terrain_slope_deg,RGI_selected_deg[closeGLID,]),extent_tiploc),na.rm=T)
# Distance to closest glacier
lakes_loc$mindisttoglac[i] <- min(pts.wit.dist[,3])
}else{
  lakes_loc$slopelowertongue[i] <- NA
  lakes_loc$mindisttoglac[i] <- NA
}
# Lake elevation
lake_elev <- cellStats(mask(LDEM_deg,lakes_loc_deg[i,]),'min')

# Seismic Risk
lake_seis <- cellStats(mask(resample(LSeismic_deg,LDEM_deg),lakes_loc_deg[i,]),'mean')

# Precipitation Intensity
lake_prec <- cellStats(mask(resample(LPrecip_deg,LDEM_deg),lakes_loc_deg[i,]),'mean')

# Get buffer around lake for slope extraction
b <- gBuffer(lakes_loc[i,], width=300, quadsegs=100)
e <- erase(b, lakes_loc[i,])                             # using raster:erase
e <- rgeos::gDifference(b, lakes_loc[i,], byid=TRUE)     # using rgeos:gDifference
e_deg <- spTransform(e,CRS("+proj=longlat +datum=WGS84"))

slopeHeadWalls <- terrain_slope_deg
slopeDam <- terrain_slope_deg
slopeHeadWalls[which(LDEM_deg[]<(lake_elev+10))] <- NA
slopeDam[which(LDEM_deg[]>(lake_elev))] <- NA
adjacentSlope <- mask(slopeHeadWalls,e_deg)
slopeDam <- mask(slopeDam,e_deg)
lakes_loc$adjecentheadwalls[i] <- cellStats(adjacentSlope,'mean')
lakes_loc$slopeDam[i] <- cellStats(slopeDam,'mean')

# Avalanche/Landslide Pathways

steeppx <- which(slopeHeadWalls[]>30)   # steep pixels above the lake

if(length(steeppx)>0){
#plot(LDEM_deg,legend=F)
r.fd <- terrain(LDEM_deg, opt='flowdir')
#plot(lakes_loc_deg[i,],add=T,col='blue')
LSsources <- vector()
for(k in 1:length(steeppx)){
  r.p <- flowPath(r.fd, coordinates(LDEM_deg)[steeppx[k],])
  if(is.null(r.p)){}else{
    p.xy <- xyFromCell(r.fd,r.p)
    #lines(p.xy,col='red', lwd=4)
    LSsources[k] <- length(which(pointsInPolygon(p.xy, lakes_loc_deg[i,], logical = T)==T))
  }
}

if(length(which(LSsources>0))>5& quantile(adjacentSlope, probs = c(0.75)) > 30){
  lakes_loc$landslide_score[i] <- 1
}else{
  lakes_loc$landslide_score[i] <- 0
}
}else{
  lakes_loc$landslide_score[i] <- 0
}



# Make output Table

lakes_loc$seismic[i] <- lake_seis
seisScore <- lake_seis/9.81
if(is.na(seisScore)){
  lakes_loc$seismic_score[i] <- NA
}else if(seisScore>0.18){
  lakes_loc$seismic_score[i] <- 1
}else{
  lakes_loc$seismic_score[i] <- 0
}

if(is.na(lake_prec)){
  lakes_loc$precip_score[i] <- NA
}else if(lake_prec<=50){
  lakes_loc$precip_score[i] <- 0
}else{
  lakes_loc$precip_score[i] <- 1
}

lakes_loc$precip[i] <- lake_prec
#}
}

writeOGR(lakes_loc, dsn = path_lakes&'\\'&lake_filename_shp&'updated' , layer = "lakes_loc")