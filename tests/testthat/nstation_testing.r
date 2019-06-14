
"	DATA TO DEAL WITH

	real<lower=0> dt;		// length (in minutes) of a time step
	vector<lower=0> [nSites] DOinitial;
	vector<lower = 0> [nSites] latWeight;
	vector<lower=0> [nSites] latInputDO;

	ALSO TEST 
		transformed_parameters
"

library(RSQLite)
library(sp)
library(maptools)
library(raster)
library(lubridate)
library(data.table)
library(ggplot2)
# library(WatershedTools)
# library(NSmetabolism)
library(reshape2)
library(rstan)
setwd("~/work/packages/NSmetabolism")
devtools::load_all("../WatershedTools")

# dbPath <- "~/work/projects/metabolism/metabolismDB/metabolism.sqlite"
# gisPath <- "~/work/projects/metabolism/catchment_delineations/vjosa/res/"
dbPath <- "~/work/projects/metabolismDB/metabolism.sqlite"
gisPath <- "~/work/projects/catchment_delineations/vjosa/res/"
lightPath <- "/Volumes/Data/vjosa_sun/"

vjosa <- readRDS(file.path(gisPath, "vjosaWatershedSpring2018.rds"))

metabDB <- dbConnect(RSQLite::SQLite(), dbPath)
dbExecute(metabDB, "PRAGMA foreign_keys = ON")

expedition <- "VjosaSpring2018"
qry <- paste0("SELECT DISTINCT siteID, siteName, x, y FROM allSiteData WHERE expedition LIKE '", 
	expedition, "'")
sites <- data.table(dbGetQuery(metabDB, qry))
coordinates(sites) <- c('x', 'y')
proj4string(sites) <- CRS("+init=epsg:4326")

qry <- "SELECT * FROM siteAttributeView WHERE variable LIKE 'elevation'"
siteElev <- data.table(dbGetQuery(metabDB, qry))



# plot the river
riverChannel <- readRDS(file.path(gisPath, "stream_vector.rds"))
sites <- spTransform(sites, CRS=proj4string(riverChannel))
plot(riverChannel, col='#377eb8')
cols <- rep('#984ea3', nrow(sites))
cols[grep('trekking_hellas|ws_albania', sites$siteName)] <- "#ff7f00"
plot(sites, add=TRUE, col=cols, pch=16, cex=0.8)
maptools::pointLabel(coordinates(sites),labels=sites$siteName, cex=0.6)


# choose a smaller subcatchment to work with
dsSite <- "37"
dsSiteInd <- match(dsSite, sites$siteName)

### starting/ending times
startDate <- dmy("24-04-18")
endDate <- dmy("29-04-18")



sarantopouros <- subcatchment(vjosa, WatershedTools::extract(vjosa, sites[dsSiteInd,]))
sites <- sites[!is.na(WatershedTools::extract(sarantopouros, sites)),]
sites$pixelID <- WatershedTools::extract(sarantopouros, sites)

# get all data measurements for desired sites
qry <- paste0("SELECT timestamp, variable, value, siteID, siteName FROM allSiteData 
	WHERE expedition LIKE '", expedition, "'")
allDat <- data.table(dbGetQuery(metabDB, qry))
allDat$time <- as_datetime(allDat$timestamp, tz="Europe/Tirane")
allDat <- allDat[date(time) >= startDate & date(time) >= endDate]
allDat$minutes <- 1 + (allDat$timestamp - min(allDat$timestamp)) / 60

prDat <- allDat[variable == 'pressure']

allDat <- merge(allDat, sites, by = 'siteID', all.x = FALSE, all.y = TRUE)
allDat$siteID <- factor(allDat$siteID)


#################
#
#  DISSOLVED OXYGEN
#
#################

doDat <- allDat[allDat$variable %in% c('water temperature', 'dissolved oxygen')]
ggplot(doDat, aes(x = time, y = value, col=siteID)) + geom_line(size = 0.5) + 
	facet_grid(variable ~ .)
doDat <- dcast(doDat, siteName.x + siteID + timestamp + x + y + time + minutes ~ variable)
doDat$pixelID <- sites[match(doDat$siteID, sites$siteID),]$pixelID

#################
#
#  WATER TEMPERATURE
#
#################

maxTime <- max(doDat$minutes)
# water temperature needs to be every minute
watTempDat <- acast(doDat, siteName.x ~ minutes, value.var = 'water temperature')
watTempDatList <- apply(watTempDat, 1, function(x) {
	approx(as.integer(colnames(watTempDat)), x, 1:maxTime, rule = 2)
})
watTempDatInterp <- do.call(rbind, lapply(watTempDatList, function(x) x$y))
colnames(watTempDatInterp) <- watTempDatList[[1]][['x']]
rownames(watTempDatInterp) <- names(watTempDatList)
waterTempSiteIDs <- sites$pixelID[match(rownames(watTempDatInterp), sites$siteName)]
wtDmat <- wsDistance(sarantopouros, waterTempSiteIDs)
wtNBs <- nearestNeighbors(sarantopouros, sarantopouros[,'id'], wtDmat, selfAdjacency = TRUE)


## note - this might be a useful addition to nearestNeighbors
wtNBMat <- cbind(wtNBs$downstream[sort(wtNBs$downstream[,1]),2], NA, NA)
rownames(wtNBMat) <- sort(wtNBs$downstream[,1])
for(si in unique(wtNBs$upstream[,1]))
	wtNBMat[si,2:3] <- wtNBs$upstream[wtNBs$upstream[,1] == si,2]
wtNBMat[wtNBMat[,2] == wtNBMat[,3], 3] <- NA
wtDistanceMatrix <- matrix(NA, nrow = nrow(wtNBMat), ncol = ncol(wtNBMat))
for(i in 1:nrow(wtNBMat)) {
	j <- match(as.character(wtNBMat[i,]), rownames(wtDmat))
	wtDistanceMatrix[i,] <- wtDmat[j, i]
}
wtDistanceMatrix <- abs(wtDistanceMatrix)
waterTempIndices <- apply(wtNBMat, 2, function(x) match(x, waterTempSiteIDs))

#################
#
#  PRESSURE
#
#################
prDatMat <- acast(prDat, siteName ~ minutes, value.var = 'value')
prDatList <- apply(prDatMat, 1, function(x) {
	approx(as.integer(colnames(prDatMat)), x, 1:maxTime, rule = 2)
})
prDatInterp <- do.call(rbind, lapply(prDatList, function(x) x$y))
colnames(prDatInterp) <- prDatList[[1]][['x']]
rownames(prDatInterp) <- names(prDatList)
prElev = siteElev[match(rownames(prDatInterp), siteElev$siteName),value]
qry <- "SELECT x,y,siteName FROM sites"
prCoords <- data.table(dbGetQuery(metabDB, qry))
prCoords <- prCoords[siteName %in% rownames(prDatMat), .(x, y)]

#################
#
#  LIGHT
#
#################

# parse the day numbers into dates to look for light files
# ltFiles <- as.vector(sapply(paste0("vj_irradiance_d", unique(yday(doDat$time)), "_.+\\.tif"), 
# 	function(x) list.files(lightPath, pattern = x, full.names = TRUE)))
# light <- stack(ltFiles)
# # grab times out of file names, construct into proper date objects
# dt <- date("2017-12-31") + as.integer(sub(".+_d([0-9]+)_.+\\.tif", "\\1", ltFiles))
# hr <- sub(".+_t([0-9]+)_.+\\.tif", "\\1", ltFiles)
# mn <- ceiling(60*(as.integer(
# 	sub("^0$", "00", sub(".+_t[0-9]+_([0-9]+)\\.tif", "\\1", ltFiles)))/100))
# lightTimes <- ymd_hm(paste(dt, paste(hr, mn, sep=':')), tz = "Europe/Tirane")
# lightMins <- 1 + (as.integer(lightTimes) - min(allDat$timestamp)) / 60
# lightVals <- raster::extract(light, sarantopouros$data[,c('x', 'y')])
# rownames(lightVals) <- sarantopouros[,'id']
# colnames(lightVals) <- lightMins

# # need to add 0s during night
# # miuntes digit
# mDigit <- min(lightMins) - round(min(lightMins), -1)
# minsExpected <- expand.grid(id = sarantopouros[,'id'], 
# 	minute = seq(mDigit, max(doDat$minutes), 10))
# minsExpected$noLight <- 0

# lvTall <- melt(lightVals, varnames = c("id", "minute"))
# lvTallMerge <- merge(lvTall, minsExpected, all = TRUE, by = c('id', 'minute'))
# lvTallMerge$value[is.na(lvTallMerge$value)] <- lvTallMerge$noLight[is.na(lvTallMerge$value)]
# lightVals <- acast(lvTallMerge, id ~ minute, value.var = "value")
# saveRDS(lightMins, "inst/testdata/lightMins.rds")
# saveRDS(lightVals, "inst/testdata/lightVals.rds")

lightVals <- readRDS(system.file('testdata/lightVals.rds', package="NSmetabolism"))
lightMins <- readRDS(system.file('testdata/lightMins.rds', package="NSmetabolism"))


#############
#
# PIXEL NEIGHBORS & LATERAL INPUT
#
##############
adjTall <- (which(sarantopouros$adjacency == 1, arr.ind = TRUE))
usNb <- t(sapply(sarantopouros[,'id'], function(i) {
	vals <- adjTall[adjTall[,1] == i, 2]
	if(length(vals) == 0) {
		c(NA, NA)
	} else if(length(vals) == 1) {
		c(vals, NA)
	} else
		vals
}))
usWeight <- usNb
usWeight[,1] <- sarantopouros[usNb[,1], 'discharge']
usWeight[,2] <- sarantopouros[usNb[,2], 'discharge']
usWeight[is.na(usWeight)] <- 0

# we compute the weight for lateral input as the difference between focal and upstream discharge
latWeight <- sarantopouros[,'discharge'] - rowSums(usWeight)


#############
#
# STARTING VALUES FOR DO
#
##############
DOinitial <- NA * numeric(nrow(sarantopouros$data))
doi <- t(sapply(unique(doDat$pixelID), function(id) {
	i <- doDat$pixelID == id
	c(id, doDat[['dissolved oxygen']][i & doDat$minutes == min(doDat$minutes[i])])
}))

# we already have distance matrix and nearest neighbors, so use those
# first, handle sites that actually have observations
DOinitial[doi[,1]] <- doi[,2]

# now anything with only one neighbor; just set to neighbor value
i <- apply(wtDistanceMatrix, 1, function(x) is.na(x[2] & is.na(x[3]))) & is.na(DOinitial)
DOinitial[i] <- wtDistanceMatrix[i,1]

# for the rest, inverse distance weighting with discharge ratio
i <- is.na(DOinitial)
wts <- apply(wtNBMat, 2, function(x) sarantopouros[,'discharge'][x])
wts <- apply(wts, 2, function(x) x/sarantopouros[,'discharge'])
# upstream is always on the numerator
wts[,1] <- 1/wts[,1]
wts <- wts / (wtDistanceMatrix^2)
wts[is.na(wts)] <- 0
wts <- wts / rowSums(wts)
vals <- wts * apply(wtNBMat, 2, function(x) doi[match(x, doi[,1]),2])
vals[is.na(vals)] <- 0
vals <- rowSums(vals)
DOinitial[i] <- vals[i]



## take care of NAs in data; this is SUPER hacky and I should make this better
for(i in 2:3) {
	j <- is.na(wtNBMat[,i])
	wtNBMat[j,i] <- wtNBMat[j,1]
	j <- is.na(wtDistanceMatrix[,i])
	wtDistanceMatrix[j,i] <- wtDistanceMatrix[j,1]
}
usNb[is.na(usNb)] <- 1
stanDat <- list(
	nDO = nrow(doDat),
	nSites = nrow(watTempDatInterp),
	nPressure = nrow(prDatInterp),
	nPixels = nrow(sarantopouros$data),
	nReaches = length(unique(sarantopouros[,'reachID'])),
	nLightTimes = ncol(lightVals),
	coords = sarantopouros[,c('x', 'y')],
	elevation = sarantopouros[,'elevation'],
	Q = sarantopouros[,'discharge'],
	reachID = sarantopouros[,'reachID'],
	depth = sarantopouros[,'depth'],
	area = sarantopouros[,'width']*sarantopouros[,'depth'],
	dx = sarantopouros[,'length'],
	DOinitial = DOinitial,
	usWeight = usWeight,
	usNb = usNb,
	latWeight = latWeight,
	latInputDO = rep(0, length(latWeight)),
	light = lightVals,
	lightTimes = as.integer(colnames(lightVals)),
	logP1_pr_mu = 9,
	logP2_pr_mu = 9,
	logP1_pr_sd = 1,
	logP2_pr_sd = 1,
	slope = tapply(sarantopouros[,'slope'], sarantopouros[,'reachID'], mean),
	velocity = tapply(sarantopouros[,'velocity'], sarantopouros[,'reachID'], mean),
	maxTime = ncol(watTempDatInterp),
	DO = doDat[['dissolved oxygen']],
	DOpixels = doDat[['pixelID']],
	DOtimes = doDat[['minutes']],
	waterTempMeasured = watTempDatInterp,
	waterTempSiteIDs = waterTempSiteIDs,
	waterTempNbs = wtNBMat,
	waterTempDist = wtDistanceMatrix,
	waterTempIndices = waterTempIndices,
	pressure = prDatInterp,
	prCoords = prCoords,
	prElev = as.array(prElev),
	dummyOsat = matrix(rnorm(length(prDatInterp)), nrow = nrow(prDatInterp), 
		ncol=ncol(prDatInterp)),
	dt = 1)


devtools::load_all(".")
modCode <- stanc_builder(file = 
	system.file("stan/nstation_reach_test.stan", package="NSmetabolism"))
mod <- stan_model(stanc_ret = modCode)
fit <- sampling(mod, data = stanDat, chains = 1, iter = 100)
