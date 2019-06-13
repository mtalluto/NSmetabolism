
"	DATA TO DEAL WITH

	int<lower = 1> nReaches;
	real<lower=0> dt;		// length (in minutes) of a time step

	vector<lower=0> [nSites] DOinitial;
	vector<lower=0> [nSites] area;
	vector<lower=0> [nSites] dx;
	vector<lower=0> [nSites] depth;
	int<lower=1, upper = nReaches> reachID [nSites];

	vector<lower = 0> [nReaches] slope;
	vector<lower = 0> [nReaches] velocity;

	
	matrix<lower=0> [nSites, 2] usWeight;
	int<lower = 1, upper = nSites> usNb [nSites, 2]; // neighbor indices for upstream pixels

	vector<lower = 0> [nSites] latWeight;
	vector<lower=0> [nSites] latInputDO;


	int <lower = 1> nLightTimes;
	matrix [nSites, nLightTimes] light;
	vector [nLightTimes] lightTimes;

	real logP1_pr_mu; // suggested from 1 station - mean of 9
	real logP1_pr_sd; // suggested from 1 station - sd of 1
	real logP2_pr_mu; 
	real logP2_pr_sd;

"

"	functions to test
	real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx);
	real computeGPP(real PAR, real lP1, real lP2);
	real computeER(real temp, real ER24_20);

	real approx(vector x, vector y, xnew) {


	ALSO TEST 
		transformed_data (can do in line with functions where appropriate)
		transformed_parameters
		finally model
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

dbPath <- "~/work/projects/metabolism/metabolismDB/metabolism.sqlite"
gisPath <- "~/work/projects/metabolism/catchment_delineations/vjosa/res/"

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

stanDat <- list(
	nDO = nrow(doDat),
	nSites = nrow(watTempDatInterp),
	nPressure = nrow(prDatInterp),
	nPixels = nrow(sarantopouros$data),
	nReaches = 1,
	coords = sarantopouros[,c('x', 'y')],
	elevation = sarantopouros[,'elevation'],
	Q = sarantopouros[,'discharge'],
	maxTime = ncol(watTempDatInterp),
	DO = doDat[['dissolved oxygen']],
	DOpixels = doDat[['pixelID']],
	DOtimes = doDat[['minutes']],
	waterTempMeasured = watTempDatInterp,
	waterTempSiteIDs = waterTempSiteIDs,
	waterTempNbs = wtNBMat,
	waterTempDist = wtDistanceMatrix,
	pressure = prDatInterp,
	prCoords = prCoords,
	prElev = as.array(prElev),
	dummyOsat = matrix(rnorm(length(prDatInterp)), nrow = nrow(prDatInterp), 
		ncol=ncol(prDatInterp)))


devtools::load_all(".")
modCode <- stanc_builder(file = 
	system.file("stan/nstation_reach_test.stan", package="NSmetabolism"))
mod <- stan_model(stanc_ret = modCode)
fit <- sampling(mod, data = stanDat, chains = 1, iter = 100)
