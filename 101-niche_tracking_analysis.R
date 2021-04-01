
# some code for timing how long this takes
knitr::opts_chunk$set(message = FALSE)
start_time <- Sys.time()


library(geosphere)
library(fields)
library(MASS)
library(raster)
library(move)
library(igraph)
library(sp)
library(rgeos)
library(dggridR)
library(rworldmap)
library(mapplots)
library(ebirdst)
library(car)
library(ade4)
library(emdist)
library(reticulate)
library(ecospat)
library(ggplot2)
library(inlabru)
library(tidyverse)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)

dir.create("outputs/101", showWarnings = FALSE, recursive = TRUE)

##  Construct an hexagon grid covering the region of interest  ##

hexgrid <- dgconstruct(projection="ISEA", topology="HEXAGON", res=9, metric=T)
hexgrid_center <- dgSEQNUM_to_GEO(hexgrid, 1:196832)
hexgrid_centroids <- cbind(hexgrid_center$lon_deg, hexgrid_center$lat_deg)
hex_sel <- which(hexgrid_centroids[,1] < -30 & hexgrid_centroids[,1] > -170 & hexgrid_centroids[,2] > -60 & hexgrid_centroids[,2] < 80)
hexgrid2_stem <- dgcellstogrid(hexgrid, hex_sel, frame=F,wrapcells=TRUE)
hexgrid2_stem_centroids <- matrix(unlist(lapply(hexgrid2_stem@polygons, function(x) x@labpt)), byrow=T, ncol=2)
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, proj4string(hexgrid2_stem))
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- gBuffer(newmap, byid=TRUE, width=0)
newmap <- gUnaryUnion(newmap, id=newmap@data$world)
hexgrid3_stem <- gIntersection(hexgrid2_stem, newmap, byid=T)
hexgrid3_stem_centroids <- matrix(unlist(lapply(hexgrid3_stem@polygons, function(x) x@labpt)), byrow=T, ncol=2)
sr <- proj4string(hexgrid3_stem)


##  Load climate data  ##

months <- sprintf("%02d", 1:12)
names(months) <- months
# make a rasterStack of temperature months
Temp_stack <- lapply(
  months,
  function(m) raster(paste0("Climate-data/wc2.1_2.5m_tavg_", m, ".tif"))
) %>%
  stack()
Prec_stack <- lapply(
    months,
    function(m) raster(paste0("Climate-data/wc2.1_2.5m_prec_", m, ".tif"))
  ) %>%
  stack()

Temp.raster.winter <- calc(Temp_stack[[c(11:12,1:4)]], fun = mean)
Temp.raster.summer <- calc(Temp_stack[[6:8]], fun = mean)
Prec.raster.winter <- calc(Prec_stack[[c(11:12,1:4)]], fun = mean)
Prec.raster.summer <- calc(Prec_stack[[6:8]], fun = mean)


##  Extract climate data onto the hexagon grid  ##

Temp.hex.winter <- raster::extract(Temp.raster.winter, hexgrid3_stem, fun=mean, na.rm=T)
Temp.hex.summer <- raster::extract(Temp.raster.summer, hexgrid3_stem, fun=mean, na.rm=T)
Prec.hex.winter <- raster::extract(Prec.raster.winter, hexgrid3_stem, fun=mean, na.rm=T)
Prec.hex.summer <- raster::extract(Prec.raster.summer, hexgrid3_stem, fun=mean, na.rm=T)


##  Standardize climate data  ##

Temp.winter.zscore <- (Temp.hex.winter - mean(c(Temp.hex.winter, Temp.hex.summer), na.rm=T)) / sd(c(Temp.hex.winter, Temp.hex.summer), na.rm=T)
Temp.summer.zscore <- (Temp.hex.summer - mean(c(Temp.hex.winter, Temp.hex.summer), na.rm=T)) / sd(c(Temp.hex.winter, Temp.hex.summer), na.rm=T)
Prec.winter.zscore <- (Prec.hex.winter - mean(c(Prec.hex.winter, Prec.hex.summer), na.rm=T)) / sd(c(Prec.hex.winter, Prec.hex.summer), na.rm=T)
Prec.summer.zscore <- (Prec.hex.summer - mean(c(Prec.hex.winter, Prec.hex.summer), na.rm=T)) / sd(c(Prec.hex.winter, Prec.hex.summer), na.rm=T)


##  Get WIFL range map and convert it to presence/absence on the hexagon grid  ##

WIFL_rangemap <- readOGR("data/WIFL_RangeMap", "WIFLrev", verbose=F)
WIFL_rangemap <- spTransform(WIFL_rangemap, proj4string(hexgrid3_stem))
WIFL_rangemap <- gBuffer(WIFL_rangemap, byid=TRUE, width=0)
WIFL_rangemap <- gUnaryUnion(WIFL_rangemap)
WIFL <- gIntersection(WIFL_rangemap, hexgrid3_stem, byid=T)


##  Get breeding population distribution from genoscape  ##

BreedingPops <- readRDS("stored_results/genoscape_brick.rds")
crs(BreedingPops) <- sr

# Eastern population
EST <- raster::extract(BreedingPops[[2]], hexgrid3_stem, fun=mean, na.rm=T)
EST[which(is.na(EST) == T)] <- 0
EST[which(EST < 0.5)] <- 0
EST[which(EST >= 0.5)] <- 1

# Interior North West population
INW <- raster::extract(BreedingPops[[1]], hexgrid3_stem, fun=mean, na.rm=T)
INW[which(is.na(INW) == T)] <- 0
INW[which(INW < 0.5)] <- 0
INW[which(INW >= 0.5)] <- 1

# Pacific North West population
PNW <- raster::extract(BreedingPops[[3]], hexgrid3_stem, fun=mean, na.rm=T)
PNW[which(is.na(PNW) == T)] <- 0
PNW[which(PNW < 0.5)] <- 0
PNW[which(PNW >= 0.5)] <- 1

# South West population
SSW <- raster::extract(BreedingPops[[4]], hexgrid3_stem, fun=mean, na.rm=T)
SSW[which(is.na(SSW) == T)] <- 0
SSW[which(SSW < 0.5)] <- 0
SSW[which(SSW >= 0.5)] <- 1


##  Location and population assignment of sampled individuals  ##

individual.assignements <- read.csv("data/appendices/app1-Individual_Assignment_All.csv", header=T, sep=",")

# Breeding individuals
breeding.assignments <- individual.assignements[which(individual.assignements$Stage == "Breeding"),]
ptsBR <- as.matrix(cbind(breeding.assignments$Long, breeding.assignments$Lat)) # location of breeding individuals
sptsBR <- SpatialPoints(ptsBR)
proj4string(sptsBR) <- proj4string(hexgrid3_stem)
breedingPoints_inHexagons <- gContains(hexgrid3_stem, sptsBR, byid=T)
breedingHexagons <- apply(breedingPoints_inHexagons, 1, function(x) which(x==TRUE))
toKeep <- which(lapply(breedingHexagons, length) > 0)
breedingHexagons <- unlist(breedingHexagons[toKeep]) # hexagon occupied by breeding individuals
breedingHexagons2 <- unique(breedingHexagons)
popBR <- as.character(breeding.assignments$Assignment)
popBR <- popBR[toKeep] # population assignment of breeding individuals

# Wintering individuals
wintering.assignments <- individual.assignements[which(individual.assignements$Stage == "Wintering"),]
ptsNB <- as.matrix(cbind(wintering.assignments$Long, wintering.assignments$Lat)) # location of wintering individuals
sptsNB <- SpatialPoints(ptsNB)
proj4string(sptsNB) <- proj4string(hexgrid3_stem)
winteringPoints_inHexagons <- gContains(hexgrid3_stem, sptsNB, byid=T)
winteringHexagons <- apply(winteringPoints_inHexagons, 1, function(x) which(x==TRUE))
toKeep <- which(lapply(winteringHexagons, length) > 0)
winteringHexagons <- unlist(winteringHexagons[toKeep]) # hexagon occupied by wintering individuals
winteringHexagons2 <- unique(winteringHexagons)
popNB <- as.character(wintering.assignments$Assignment)
popNB <- popNB[toKeep] # population assignment of wintering individuals
ptsNB <- ptsNB[toKeep,]

##   Estimating seasonal climatic niches  ##

# Function to estimate seasonal climatic niche
nicheDensityRaster <- function(seasonalNiche){
  niche.kernel <- kde2d(seasonalNiche[,1], seasonalNiche[,2], n=50, h=1, lims=c(-1.5,3, -1.5,3))
  niche.kernel$z = niche.kernel$z/max(niche.kernel$z)
  niche.raster <- raster(niche.kernel)
  threshold=0; i=0
  while(threshold <= 0.95 * sum(niche.kernel$z)){
    i=i+1
    threshold = threshold + sort(as.vector(niche.raster), decreasing=T)[i]
  }
  niche.raster[which(as.vector(niche.raster) < sort(as.vector(niche.raster), decreasing=T)[i])] = 0
  niche.raster = niche.raster / sum(as.vector(niche.raster))
  return(niche.raster)
}

# Estimating breeding niches
breeding.niche.INW <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "INW")]], Prec.summer.zscore[breedingHexagons[which(popBR == "INW")]]))
breeding.niche.EST <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "EST")]], Prec.summer.zscore[breedingHexagons[which(popBR == "EST")]]))
breeding.niche.PNW <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "PNW")]], Prec.summer.zscore[breedingHexagons[which(popBR == "PNW")]]))
breeding.niche.SSW <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "SSW")]], Prec.summer.zscore[breedingHexagons[which(popBR == "SSW")]]))
breeding.niche.WIFL <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons], Prec.summer.zscore[breedingHexagons])) # Entire species

# Estimating wintering niches
wintering.niche.INW <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "INW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "INW")]]))
wintering.niche.EST <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "EST")]], Prec.winter.zscore[winteringHexagons[which(popNB == "EST")]]))
wintering.niche.PNW <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "PNW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "PNW")]]))
wintering.niche.SSW <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "SSW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "SSW")]]))
wintering.niche.WIFL <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons], Prec.winter.zscore[winteringHexagons])) # Entire species

# Estimating total niches (year-round)
total.niche.INW <- nicheDensityRaster(cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "INW")]], Temp.winter.zscore[winteringHexagons[which(popNB == "INW")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "INW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "INW")]])))
total.niche.EST <- nicheDensityRaster(cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "EST")]], Temp.winter.zscore[winteringHexagons[which(popNB == "EST")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "EST")]], Prec.winter.zscore[winteringHexagons[which(popNB == "EST")]])))
total.niche.PNW <- nicheDensityRaster(cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "PNW")]], Temp.winter.zscore[winteringHexagons[which(popNB == "PNW")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "PNW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "PNW")]])))
total.niche.SSW <- nicheDensityRaster(cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "SSW")]], Temp.winter.zscore[winteringHexagons[which(popNB == "SSW")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "SSW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "SSW")]])))
total.niche.WIFL <- nicheDensityRaster(cbind(c(Temp.summer.zscore[breedingHexagons], Temp.winter.zscore[winteringHexagons]), c(Prec.summer.zscore[breedingHexagons], Prec.winter.zscore[winteringHexagons]))) # Entire species

# Niche density
density.breeding.INW <- rasterToPoints(breeding.niche.INW)[,3]
density.breeding.EST <- rasterToPoints(breeding.niche.EST)[,3]
density.breeding.PNW <- rasterToPoints(breeding.niche.PNW)[,3]
density.breeding.SSW <- rasterToPoints(breeding.niche.SSW)[,3]
density.breeding.WIFL <- rasterToPoints(breeding.niche.WIFL)[,3]
density.wintering.INW <- rasterToPoints(wintering.niche.INW)[,3]
density.wintering.EST <- rasterToPoints(wintering.niche.EST)[,3]
density.wintering.PNW <- rasterToPoints(wintering.niche.PNW)[,3]
density.wintering.SSW <- rasterToPoints(wintering.niche.SSW)[,3]
density.wintering.WIFL <- rasterToPoints(wintering.niche.WIFL)[,3]
density.total.INW <- rasterToPoints(total.niche.INW)[,3]
density.total.EST <- rasterToPoints(total.niche.EST)[,3]
density.total.PNW <- rasterToPoints(total.niche.PNW)[,3]
density.total.SSW <- rasterToPoints(total.niche.SSW)[,3]
density.total.WIFL <- rasterToPoints(total.niche.WIFL)[,3]


##  Compute seasonal niche size for every populations  ##

breeding.niche.size.INW <- length(which(density.breeding.INW > 0))
breeding.niche.size.EST <- length(which(density.breeding.EST > 0))
breeding.niche.size.PNW <- length(which(density.breeding.PNW > 0))
breeding.niche.size.SSW <- length(which(density.breeding.SSW > 0))
breeding.niche.size.WIFL <- length(which(density.breeding.WIFL > 0))

wintering.niche.size.INW <- length(which(density.wintering.INW > 0))
wintering.niche.size.EST <- length(which(density.wintering.EST > 0))
wintering.niche.size.PNW <- length(which(density.wintering.PNW > 0))
wintering.niche.size.SSW <- length(which(density.wintering.SSW > 0))
wintering.niche.size.WIFL <- length(which(density.wintering.WIFL > 0))

total.niche.size.INW <- length(which(density.total.INW > 0))
total.niche.size.EST <- length(which(density.total.EST > 0))
total.niche.size.PNW <- length(which(density.total.PNW > 0))
total.niche.size.SSW <- length(which(density.total.SSW > 0))
total.niche.size.WIFL <- length(which(density.total.WIFL > 0))


##  Compute seasonal niche overlap (D metric) for every populations  ##

niche.overlap.INW <- 1 - (0.5 * (sum(abs(density.breeding.INW - density.wintering.INW))))
niche.overlap.EST <- 1 - (0.5 * (sum(abs(density.breeding.EST - density.wintering.EST))))
niche.overlap.PNW <- 1 - (0.5 * (sum(abs(density.breeding.PNW - density.wintering.PNW))))
niche.overlap.SSW <- 1 - (0.5 * (sum(abs(density.breeding.SSW - density.wintering.SSW))))
niche.overlap.WIFL <- 1 - (0.5 * (sum(abs(density.breeding.WIFL - density.wintering.WIFL))))



###  Plotting figure  ###

##  Prepare data frame for plotting niche density kernels  ##

breedingNicheRasters <- list(
  WIFL = breeding.niche.WIFL,
  EST = breeding.niche.EST,
  INW = breeding.niche.INW,
  PNW = breeding.niche.PNW,
  SSW = breeding.niche.SSW
)
winteringNicheRasters <- list(
  WIFL = wintering.niche.WIFL,
  EST = wintering.niche.EST,
  INW = wintering.niche.INW,
  PNW = wintering.niche.PNW,
  SSW = wintering.niche.SSW
)
totalNicheRasters <- list(
  WIFL = total.niche.WIFL,
  EST = total.niche.EST,
  INW = total.niche.INW,
  PNW = total.niche.PNW,
  SSW = total.niche.SSW
)

bigList.niche <- list(`Breeding niche` = breedingNicheRasters, `Wintering niche` = winteringNicheRasters, `Total niche` = totalNicheRasters)

tibbleVar.niche <- lapply(bigList.niche, function(x){
  lapply(x, function(y){
    tib <- as.data.frame(rasterToPoints(y))
    names(tib) <- c("x", "y", "val")
    as_tibble(tib)
  }) %>%
    bind_rows(.id = "region")
}) %>%
  bind_rows(.id = "seasonality") %>%
  mutate(
    season_f = factor(seasonality, levels=c("Breeding niche", "Wintering niche", "Total niche")),
    region_f = factor(region, levels=c("WIFL", "EST", "INW", "PNW", "SSW"))
  ) %>%
  mutate(
    forColor = case_when(
      val > 0.012 ~ "col5",
      val > 0.009 ~ "col4",
      val > 0.006 ~ "col3",
      val > 0.003 ~ "col2",
      val > 0     ~ "col1",
      TRUE        ~ "zero"
    ),
    forColor = paste(region, "–", forColor)
  )




##  Prepare data frame for plotting points in niche space  ##

breedingPoints_inHexagons_EST <- over(sptsBR[which(popBR=="EST")], hexgrid3_stem[which(EST>0),], byid=T)
outsidePoint.EST <- which(is.na(breedingPoints_inHexagons_EST)==T)
breedingPoints_inHexagons_INW <- over(sptsBR[which(popBR=="INW")], hexgrid3_stem[which(INW>0),], byid=T)
outsidePoint.INW <- which(is.na(breedingPoints_inHexagons_INW)==T)
breedingPoints_inHexagons_PNW <- over(sptsBR[which(popBR=="PNW")], hexgrid3_stem[which(PNW>0),], byid=T)
outsidePoint.PNW <- which(is.na(breedingPoints_inHexagons_PNW)==T)
breedingPoints_inHexagons_SSW <- over(sptsBR[which(popBR=="SSW")], hexgrid3_stem[which(SSW>0),], byid=T)
outsidePoint.SSW <- which(is.na(breedingPoints_inHexagons_SSW)==T)

# Breeding individuals inside the population polygon to which they are assigned
breedingNichePoints.in <- list(
  WIFL = cbind(Temp.summer.zscore[breedingHexagons], Prec.summer.zscore[breedingHexagons], "in"),
  EST = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "EST")]][-outsidePoint.EST], Prec.summer.zscore[breedingHexagons[which(popBR == "EST")]][-outsidePoint.EST], "in"),
  INW = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "INW")]][-outsidePoint.INW], Prec.summer.zscore[breedingHexagons[which(popBR == "INW")]][-outsidePoint.INW], "in"),
  PNW = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "PNW")]][-outsidePoint.PNW], Prec.summer.zscore[breedingHexagons[which(popBR == "PNW")]][-outsidePoint.PNW], "in"),
  SSW = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "SSW")]][-outsidePoint.SSW], Prec.summer.zscore[breedingHexagons[which(popBR == "SSW")]][-outsidePoint.SSW], "in")
)
# Breeding individuals outside the population polygon to which they are assigned
breedingNichePoints.out <- list(
  WIFL = NULL,
  EST = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "EST")]][outsidePoint.EST], Prec.summer.zscore[breedingHexagons[which(popBR == "EST")]][outsidePoint.EST], "out"),
  INW = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "INW")]][outsidePoint.INW], Prec.summer.zscore[breedingHexagons[which(popBR == "INW")]][outsidePoint.INW], "out"),
  PNW = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "PNW")]][outsidePoint.PNW], Prec.summer.zscore[breedingHexagons[which(popBR == "PNW")]][outsidePoint.PNW], "out"),
  SSW = cbind(Temp.summer.zscore[breedingHexagons[which(popBR == "SSW")]][outsidePoint.SSW], Prec.summer.zscore[breedingHexagons[which(popBR == "SSW")]][outsidePoint.SSW], "out")
)
winteringNichePoints <- list(
  WIFL = cbind(Temp.winter.zscore[winteringHexagons], Prec.winter.zscore[winteringHexagons], "in"),
  EST = cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "EST")]], Prec.winter.zscore[winteringHexagons[which(popNB == "EST")]], "in"),
  INW = cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "INW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "INW")]], "in"),
  PNW = cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "PNW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "PNW")]], "in"),
  SSW = cbind(Temp.winter.zscore[winteringHexagons[which(popNB == "SSW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "SSW")]], "in")
)
TotalNichePoints <- list(
  WIFL = cbind(c(Temp.summer.zscore[breedingHexagons], Temp.winter.zscore[winteringHexagons]), c(Prec.summer.zscore[breedingHexagons], Prec.winter.zscore[winteringHexagons]), "in"),
  EST = cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "EST")]], Temp.winter.zscore[winteringHexagons[which(popNB == "EST")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "EST")]], Prec.winter.zscore[winteringHexagons[which(popNB == "EST")]]), "in"),
  INW = cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "INW")]], Temp.winter.zscore[winteringHexagons[which(popNB == "INW")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "INW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "INW")]]), "in"),
  PNW = cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "PNW")]], Temp.winter.zscore[winteringHexagons[which(popNB == "PNW")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "PNW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "PNW")]]), "in"),
  SSW = cbind(c(Temp.summer.zscore[breedingHexagons[which(popBR == "SSW")]], Temp.winter.zscore[winteringHexagons[which(popNB == "SSW")]]), c(Prec.summer.zscore[breedingHexagons[which(popBR == "SSW")]], Prec.winter.zscore[winteringHexagons[which(popNB == "SSW")]]), "in")
)

breedingNichePoints <- list()
breedingNichePoints$WIFL <- rbind(breedingNichePoints.in$WIFL, breedingNichePoints.out$WIFL)
breedingNichePoints$EST <- rbind(breedingNichePoints.in$EST, breedingNichePoints.out$EST)
breedingNichePoints$INW <- rbind(breedingNichePoints.in$INW, breedingNichePoints.out$INW)
breedingNichePoints$PNW <- rbind(breedingNichePoints.in$PNW, breedingNichePoints.out$PNW)
breedingNichePoints$SSW <- rbind(breedingNichePoints.in$SSW, breedingNichePoints.out$SSW)

bigList.points <- list(`Breeding niche` = breedingNichePoints, `Wintering niche` = winteringNichePoints, `Total niche` = TotalNichePoints)

tibbleVar.points <- lapply(bigList.points, function(x){
  lapply(x, function(y){
    tib <- as.data.frame(y)
    names(tib) <- c("x", "y", "inout")
    tib$x <- as.numeric(as.character(tib$x))
    tib$y <- as.numeric(as.character(tib$y))
    as_tibble(tib)
  }
  ) %>%
    bind_rows(.id = "region")
}) %>%
  bind_rows(.id = "seasonality") %>%
  mutate(
    season_f = factor(seasonality, levels=c("Breeding niche", "Wintering niche", "Total niche")),
    region_f = factor(region, levels=c("WIFL", "EST", "INW", "PNW", "SSW"))
  )


##  Plot niches  ##

pdf("outputs/101/Figure_nicheTracking3.pdf", width=6, height=10)
theme_set(theme_light())
ggplot(tibbleVar.niche, aes(x = x, y = y)) +
  geom_raster(aes(fill = forColor)) +
  scale_x_continuous(limits = c(-0.8, 1.8)) +
  scale_y_continuous(position = "right", limits = c(-1.8, 4.5)) +
  scale_fill_manual(values = c(brewer.pal(5, "Blues"), "white",
                               brewer.pal(5, "Purples"), "white",
                               brewer.pal(5, "Greens"), "white",
                               brewer.pal(5, "Oranges"), "white",
                               brewer.pal(5, "Greys"), "white"
  )
  ) +
  geom_point(data=tibbleVar.points, size=0.8, aes(shape=inout)) +
  geom_contour(data = tibbleVar.niche, aes(x = x, y = y, z = val), col="grey47", binwidth=0.003, size=0.2, lty=1) +
  geom_contour(data = tibbleVar.niche, aes(x = x, y = y, z = val), col="grey47", breaks=0.0000001, size=0.2, lty=2) +
  facet_grid(region_f ~ season_f, switch = "y") +
  labs(x="Temperature", y="Precipitation") +
  theme(legend.position = "none", strip.placement = "outside")
dev.off()


## Plot maps

world <- ne_countries(scale = "medium", returnclass = "sf")

pdf("outputs/101/Figure_nicheTracking_maps2.pdf", width=3, height=10)

theme_set(theme_classic())

## Whole species ##
p1 <- ggplot() +
  geom_sf(data = world, color="grey65", fill="grey65") +
  coord_sf(xlim = c(-140, -55), ylim = c(-5, 62), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")
# Add summer range
p1 <- p1 + geom_polygon(data=WIFL, aes(x=long, y=lat, group=group)) + geom_polygon(colour="grey47", fill="grey47")
# Add breeding and wintering locations
p1 <- p1 + geom_point(aes(x=ptsBR[,1], y=ptsBR[,2]), color="black", size=0.7, pch=1) #dark orange2
p1 <- p1 + geom_point(aes(x=ptsNB[,1], y=ptsNB[,2]), color="black", size=0.7, pch=1) #blue

## EST ##
p1.EST <- ggplot() +
  geom_sf(data = world, color="grey65", fill="grey65") +
  coord_sf(xlim = c(-140, -55), ylim = c(-5, 62), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")
# Add summer range
p1.EST <- p1.EST + geom_polygon(data=hexgrid3_stem[which(EST>0),], aes(x=long, y=lat, group=group)) + geom_polygon(colour="grey47", fill="grey47")
# Add breeding and wintering locations
p1.EST <- p1.EST + geom_point(aes(x=ptsBR[which(popBR=="EST"),][,1][-outsidePoint.EST], y=ptsBR[which(popBR=="EST"),][,2][-outsidePoint.EST]), color="#377eb8", size=0.7, shape=1)
p1.EST <- p1.EST + geom_point(aes(x=ptsBR[which(popBR=="EST"),][,1][outsidePoint.EST], y=ptsBR[which(popBR=="EST"),][,2][outsidePoint.EST]), color="#377eb8", size=0.7, shape=2)
p1.EST <- p1.EST + geom_point(aes(x=ptsNB[which(popNB=="EST"),][,1], y=ptsNB[which(popNB=="EST"),][,2]), color="#377eb8", size=0.7, shape=1)

# INW ##
p1.INW <- ggplot() +
  geom_sf(data = world, color="grey65", fill="grey65") +
  coord_sf(xlim = c(-140, -55), ylim = c(-5, 62), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")
# Add summer range
p1.INW <- p1.INW + geom_polygon(data=hexgrid3_stem[which(INW>0),], aes(x=long, y=lat, group=group)) + geom_polygon(colour="grey47", fill="grey47")
# Add breeding and wintering locations
p1.INW <- p1.INW + geom_point(aes(x=ptsBR[which(popBR=="INW"),][,1][-outsidePoint.INW], y=ptsBR[which(popBR=="INW"),][,2][-outsidePoint.INW]), color="#984ea3", size=0.7, pch=1)
p1.INW <- p1.INW + geom_point(aes(x=ptsBR[which(popBR=="INW"),][,1][outsidePoint.INW], y=ptsBR[which(popBR=="INW"),][,2][outsidePoint.INW]), color="#984ea3", size=0.7, pch=2)
p1.INW <- p1.INW + geom_point(aes(x=ptsNB[which(popNB=="INW"),][,1], y=ptsNB[which(popNB=="INW"),][,2]), color="#984ea3", size=0.7, pch=1)

## PNW ##
p1.PNW <- ggplot() +
  geom_sf(data = world, color="grey65", fill="grey65") +
  coord_sf(xlim = c(-140, -55), ylim = c(-5, 62), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")
# Add summer range
p1.PNW <- p1.PNW + geom_polygon(data=hexgrid3_stem[which(PNW>0),], aes(x=long, y=lat, group=group)) + geom_polygon(colour="grey47", fill="grey47")
# Add breeding and wintering locations
p1.PNW <- p1.PNW + geom_point(aes(x=ptsBR[which(popBR=="PNW"),][,1][-outsidePoint.PNW], y=ptsBR[which(popBR=="PNW"),][,2][-outsidePoint.PNW]), color="#4daf4a", size=0.7, pch=1)
p1.PNW <- p1.PNW + geom_point(aes(x=ptsBR[which(popBR=="PNW"),][,1][outsidePoint.PNW], y=ptsBR[which(popBR=="PNW"),][,2][outsidePoint.PNW]), color="#4daf4a", size=0.7, pch=2)
p1.PNW <- p1.PNW + geom_point(aes(x=ptsNB[which(popNB=="PNW"),][,1], y=ptsNB[which(popNB=="PNW"),][,2]), color="#4daf4a", size=0.7, pch=1)

## SSW ##
p1.SSW <- ggplot() +
  geom_sf(data = world, color="grey65", fill="grey65") +
  coord_sf(xlim = c(-140, -55), ylim = c(-5, 62), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")
# Add summer range
p1.SSW <- p1.SSW + geom_polygon(data=hexgrid3_stem[which(SSW>0),], aes(x=long, y=lat, group=group)) + geom_polygon(colour="grey47", fill="grey47")
# Add breeding and wintering locations
p1.SSW <- p1.SSW + geom_point(aes(x=ptsBR[which(popBR=="SSW"),][,1][-outsidePoint.SSW], y=ptsBR[which(popBR=="SSW"),][,2][-outsidePoint.SSW]), color="#ff7f00", size=0.7, pch=1)
p1.SSW <- p1.SSW + geom_point(aes(x=ptsBR[which(popBR=="SSW"),][,1][outsidePoint.SSW], y=ptsBR[which(popBR=="SSW"),][,2][outsidePoint.SSW]), color="#ff7f00", size=0.7, pch=2)
p1.SSW <- p1.SSW + geom_point(aes(x=ptsNB[which(popNB=="SSW"),][,1], y=ptsNB[which(popNB=="SSW"),][,2]), color="#ff7f00", size=0.7, pch=1)

gridExtra::grid.arrange(p1, p1.EST, p1.INW, p1.PNW, p1.SSW, nrow = 5)

dev.off()





##  Niche similarity tests  ##

dist.mat.NB <- rdist.earth(hexgrid3_stem_centroids[winteringHexagons,], miles=F)
dist.mat.NB[which(dist.mat.NB < 1, arr.ind=T)] <- 0
dist.mat.BR <- rdist.earth(hexgrid3_stem_centroids[breedingHexagons,], miles=F)
dist.mat.BR[which(dist.mat.BR < 1, arr.ind=T)] <- 0

# EST
niche.overlap.randNB.EST <- list()
for(i in 1:1000){
  winteringHexagons.randNB.EST <- sample(winteringHexagons, 1)
  sampleID <- which(names(winteringHexagons) == names(winteringHexagons.randNB.EST))
  #winteringHexagons.randNB.EST <- c(winteringHexagons.randNB.EST, winteringHexagons[-sampleID][order(dist.mat.NB[sampleID, -sampleID])[1:(length(which(popNB == "EST"))-1)]])
  winteringHexagons.randNB.EST <- c(winteringHexagons.randNB.EST, sample(winteringHexagons[-sampleID], length(which(popNB == "EST"))-1, replace=F, prob = 1 / (rank(dist.mat.NB[sampleID, -sampleID])^2) ))  #1 - exp(-1/(1+dist.mat.NB[sampleID, -sampleID]))))
  #plot(hexgrid3_stem_centroids[winteringHexagons,])
  #points(hexgrid3_stem_centroids[winteringHexagons.randNB.EST,], col="red", pch=20)
  wintering.niche.randNB.EST <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons.randNB.EST], Prec.winter.zscore[winteringHexagons.randNB.EST]))
  density.wintering.randNB.EST <- rasterToPoints(wintering.niche.randNB.EST)[,3]
  niche.overlap.randNB.EST[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.EST - density.wintering.randNB.EST))))
}
niche.overlap.randBR.EST <- list()
for(i in 1:1000){
  breedingHexagons.randBR.EST <- sample(breedingHexagons, 1)
  sampleID <- which(names(breedingHexagons) == names(breedingHexagons.randBR.EST))
  breedingHexagons.randBR.EST <- c(breedingHexagons.randBR.EST, sample(breedingHexagons[-sampleID], length(which(popBR == "EST"))-1, replace=F, prob = 1 / (rank(dist.mat.BR[sampleID, -sampleID])^2) ))  # 1 / (1 + dist.mat.BR[sampleID, -sampleID]) ))  #1 - exp(-1/(1+dist.mat.BR[sampleID, -sampleID]))))
  breeding.niche.randBR.EST <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons.randBR.EST], Prec.summer.zscore[breedingHexagons.randBR.EST]))
  density.breeding.randBR.EST <- rasterToPoints(breeding.niche.randBR.EST)[,3]
  niche.overlap.randBR.EST[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.randBR.EST - density.wintering.EST))))
}
niche.overlap.rand.EST <- c(unlist(niche.overlap.randNB.EST), unlist(niche.overlap.randBR.EST))
es.randBR.EST <- (niche.overlap.EST - mean(unlist(niche.overlap.randBR.EST))) / sd(unlist(niche.overlap.randBR.EST))
es.randNB.EST <- (niche.overlap.EST - mean(unlist(niche.overlap.randNB.EST))) / sd(unlist(niche.overlap.randNB.EST))
es.rand.EST <- (niche.overlap.EST - mean(unlist(niche.overlap.rand.EST))) / sd(unlist(niche.overlap.rand.EST))
pval.randBR.EST <- length(which(niche.overlap.randBR.EST > niche.overlap.EST)) / length(niche.overlap.randBR.EST)
pval.randNB.EST <- length(which(niche.overlap.randNB.EST > niche.overlap.EST)) / length(niche.overlap.randNB.EST)
pval.rand.EST <- length(which(niche.overlap.rand.EST > niche.overlap.EST)) / length(niche.overlap.rand.EST)


# INW
niche.overlap.randNB.INW <- list()
for(i in 1:1000){
  winteringHexagons.randNB.INW <- sample(winteringHexagons, 1)
  sampleID <- which(names(winteringHexagons) == names(winteringHexagons.randNB.INW))
  winteringHexagons.randNB.INW <- c(winteringHexagons.randNB.INW, sample(winteringHexagons[-sampleID], length(which(popNB == "INW"))-1, replace=F, prob = 1 / (rank(dist.mat.NB[sampleID, -sampleID])^2) ))  #1 - exp(-1/(1+dist.mat.NB[sampleID, -sampleID]))))
  wintering.niche.randNB.INW <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons.randNB.INW], Prec.winter.zscore[winteringHexagons.randNB.INW]))
  density.wintering.randNB.INW <- rasterToPoints(wintering.niche.randNB.INW)[,3]
  niche.overlap.randNB.INW[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.INW - density.wintering.randNB.INW))))
}
niche.overlap.randBR.INW <- list()
for(i in 1:1000){
  breedingHexagons.randBR.INW <- sample(breedingHexagons, 1)
  sampleID <- which(names(breedingHexagons) == names(breedingHexagons.randBR.INW))
  breedingHexagons.randBR.INW <- c(breedingHexagons.randBR.INW, sample(breedingHexagons[-sampleID], length(which(popBR == "INW"))-1, replace=F, prob = 1 / (rank(dist.mat.BR[sampleID, -sampleID])^2) ))   #1 - exp(-1/(1+dist.mat.BR[sampleID, -sampleID]))))
  breeding.niche.randBR.INW <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons.randBR.INW], Prec.summer.zscore[breedingHexagons.randBR.INW]))
  density.breeding.randBR.INW <- rasterToPoints(breeding.niche.randBR.INW)[,3]
  niche.overlap.randBR.INW[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.randBR.INW - density.wintering.INW))))
}
niche.overlap.rand.INW <- c(unlist(niche.overlap.randNB.INW), unlist(niche.overlap.randBR.INW))
es.randBR.INW <- (niche.overlap.INW - mean(unlist(niche.overlap.randBR.INW))) / sd(unlist(niche.overlap.randBR.INW))
es.randNB.INW <- (niche.overlap.INW - mean(unlist(niche.overlap.randNB.INW))) / sd(unlist(niche.overlap.randNB.INW))
es.rand.INW <- (niche.overlap.INW - mean(unlist(niche.overlap.rand.INW))) / sd(unlist(niche.overlap.rand.INW))
pval.randBR.INW <- length(which(niche.overlap.randBR.INW > niche.overlap.INW)) / length(niche.overlap.randBR.INW)
pval.randNB.INW <- length(which(niche.overlap.randNB.INW > niche.overlap.INW)) / length(niche.overlap.randNB.INW)
pval.rand.INW <- length(which(niche.overlap.rand.INW > niche.overlap.INW)) / length(niche.overlap.rand.INW)


# PNW
niche.overlap.randNB.PNW <- list()
for(i in 1:1000){
  winteringHexagons.randNB.PNW <- sample(winteringHexagons, 1)
  sampleID <- which(names(winteringHexagons) == names(winteringHexagons.randNB.PNW))
  winteringHexagons.randNB.PNW <- c(winteringHexagons.randNB.PNW, sample(winteringHexagons[-sampleID], length(which(popNB == "PNW"))-1, replace=F, prob = 1 / (rank(dist.mat.NB[sampleID, -sampleID])^2) )) #1 - exp(-1/(1+dist.mat.NB[sampleID, -sampleID]))))
  wintering.niche.randNB.PNW <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons.randNB.PNW], Prec.winter.zscore[winteringHexagons.randNB.PNW]))
  density.wintering.randNB.PNW <- rasterToPoints(wintering.niche.randNB.PNW)[,3]
  niche.overlap.randNB.PNW[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.PNW - density.wintering.randNB.PNW))))
}
niche.overlap.randBR.PNW <- list()
for(i in 1:1000){
  breedingHexagons.randBR.PNW <- sample(breedingHexagons, 1)
  sampleID <- which(names(breedingHexagons) == names(breedingHexagons.randBR.PNW))
  breedingHexagons.randBR.PNW <- c(breedingHexagons.randBR.PNW, sample(breedingHexagons[-sampleID], length(which(popBR == "PNW"))-1, replace=F, prob = 1 / (rank(dist.mat.BR[sampleID, -sampleID])^2) )) #1 - exp(-1/(1+dist.mat.BR[sampleID, -sampleID]))))
  breeding.niche.randBR.PNW <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons.randBR.PNW], Prec.summer.zscore[breedingHexagons.randBR.PNW]))
  density.breeding.randBR.PNW <- rasterToPoints(breeding.niche.randBR.PNW)[,3]
  niche.overlap.randBR.PNW[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.randBR.PNW - density.wintering.PNW))))
}
niche.overlap.rand.PNW <- c(unlist(niche.overlap.randNB.PNW), unlist(niche.overlap.randBR.PNW))
es.randBR.PNW <- (niche.overlap.PNW - mean(unlist(niche.overlap.randBR.PNW))) / sd(unlist(niche.overlap.randBR.PNW))
es.randNB.PNW <- (niche.overlap.PNW - mean(unlist(niche.overlap.randNB.PNW))) / sd(unlist(niche.overlap.randNB.PNW))
es.rand.PNW <- (niche.overlap.PNW - mean(unlist(niche.overlap.rand.PNW))) / sd(unlist(niche.overlap.rand.PNW))
pval.randBR.PNW <- length(which(niche.overlap.randBR.PNW > niche.overlap.PNW)) / length(niche.overlap.randBR.PNW)
pval.randNB.PNW <- length(which(niche.overlap.randNB.PNW > niche.overlap.PNW)) / length(niche.overlap.randNB.PNW)
pval.rand.PNW <- length(which(niche.overlap.rand.PNW > niche.overlap.PNW)) / length(niche.overlap.rand.PNW)


# SSW
niche.overlap.randNB.SSW <- list()
for(i in 1:1000){
  winteringHexagons.randNB.SSW <- sample(winteringHexagons, 1)
  sampleID <- which(names(winteringHexagons) == names(winteringHexagons.randNB.SSW))
  winteringHexagons.randNB.SSW <- c(winteringHexagons.randNB.SSW, sample(winteringHexagons[-sampleID], length(which(popNB == "SSW"))-1, replace=F, prob = 1 / (rank(dist.mat.NB[sampleID, -sampleID])^2) )) #1 - exp(-1/(1+dist.mat.NB[sampleID, -sampleID]))))
  wintering.niche.randNB.SSW <- nicheDensityRaster(cbind(Temp.winter.zscore[winteringHexagons.randNB.SSW], Prec.winter.zscore[winteringHexagons.randNB.SSW]))
  density.wintering.randNB.SSW <- rasterToPoints(wintering.niche.randNB.SSW)[,3]
  niche.overlap.randNB.SSW[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.SSW - density.wintering.randNB.SSW))))
}
niche.overlap.randBR.SSW <- list()
for(i in 1:1000){
  breedingHexagons.randBR.SSW <- sample(breedingHexagons, 1)
  sampleID <- which(names(breedingHexagons) == names(breedingHexagons.randBR.SSW))
  breedingHexagons.randBR.SSW <- c(breedingHexagons.randBR.SSW, sample(breedingHexagons[-sampleID], length(which(popBR == "SSW"))-1, replace=F, prob = 1 / (rank(dist.mat.BR[sampleID, -sampleID])^2) )) #1 - exp(-1/(1+dist.mat.BR[sampleID, -sampleID]))))
  breeding.niche.randBR.SSW <- nicheDensityRaster(cbind(Temp.summer.zscore[breedingHexagons.randBR.SSW], Prec.summer.zscore[breedingHexagons.randBR.SSW]))
  density.breeding.randBR.SSW <- rasterToPoints(breeding.niche.randBR.SSW)[,3]
  niche.overlap.randBR.SSW[[i]] <- 1 - (0.5 * (sum(abs(density.breeding.randBR.SSW - density.wintering.SSW))))
}
niche.overlap.rand.SSW <- c(unlist(niche.overlap.randNB.SSW), unlist(niche.overlap.randBR.SSW))
es.randBR.SSW <- (niche.overlap.SSW - mean(unlist(niche.overlap.randBR.SSW))) / sd(unlist(niche.overlap.randBR.SSW))
es.randNB.SSW <- (niche.overlap.SSW - mean(unlist(niche.overlap.randNB.SSW))) / sd(unlist(niche.overlap.randNB.SSW))
es.rand.SSW <- (niche.overlap.SSW - mean(unlist(niche.overlap.rand.SSW))) / sd(unlist(niche.overlap.rand.SSW))
pval.randBR.SSW <- length(which(niche.overlap.randBR.SSW > niche.overlap.SSW)) / length(niche.overlap.randBR.SSW)
pval.randNB.SSW <- length(which(niche.overlap.randNB.SSW > niche.overlap.SSW)) / length(niche.overlap.randNB.SSW)
pval.rand.SSW <- length(which(niche.overlap.rand.SSW > niche.overlap.SSW)) / length(niche.overlap.rand.SSW)


niche.similarity.results <- data.frame(population = c("EST", "INW", "PNW", "SSW"),
                                       effect_size = c(es.randBR.EST, es.randBR.INW, es.randBR.PNW, es.randBR.SSW),
                                       pval = c(pval.randBR.EST, pval.randBR.INW, pval.randBR.PNW, pval.randBR.SSW))
write.csv(niche.similarity.results, "outputs/101/niche_similarity_results.csv", row.names = F)


# Stuff for recording how long it took to run this:
td <- Sys.time() - start_time
tdf <- hms::as_hms(ceiling(hms::as_hms(td)))
dir.create("stored_run_times", showWarnings = FALSE, recursive = TRUE)
write_rds(tdf, "stored_run_times/101.rds")
tdf

