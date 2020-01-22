context("Pixel functions")
library("NSmetabolism")
suppressWarnings(library("WatershedTools"))

test_that("Pixel-based metabolism", {
	ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))

	## eventually
	## ws$light and ws$temp needs to be a thing

	## for now, use light and other data from the db
	metabDat <- readRDS(system.file("testdata/site_metab_data/site_metab_data.rds", 
		package="NSmetabolism"))

	light <- metabDat$light
	temp <- metabDat$water_temp
	pressure <- metabDat$pressure
	pixes <- as.data.frame(ws$data)[1:2,]
	pixes$DO_i <- metabDat$DO_i
	lp1 <- metabDat$lP1
	lp2 <- metabDat$lP2
	er24_20 <- metabDat$er24_20
	k600 <- metabDat$k600

	expect_error(pixgpp <- pixelGPP(pixes, light, temp, pressure, lp1, lp2, er24_20, k600, dt=1),
		regex=NA)
	reggpp <- sapply(light[1,], computeGPP, lP1=lp1[1], lP2=lp2[1])
	expect_identical(reggpp, pixgpp[1,])

	expect_error(pixer <- pixelER(pixes, light, temp, pressure, lp1, lp2, er24_20, k600, dt=1),
		regex=NA)
	reger <- sapply(temp[1,], inSituER, ER24_20=er24_20[1])
	expect_identical(reger, pixer[1,])

	compDO <- readRDS(system.file("testdata/site_metab_data/DO_predicted_model.rds", 
		package="NSmetabolism"))
	compER <- readRDS(system.file("testdata/site_metab_data/er_predicted_model.rds", 
		package="NSmetabolism"))
	compGPP <- readRDS(system.file("testdata/site_metab_data/gpp_predicted_model.rds", 
		package="NSmetabolism"))

	expect_error(pixmetab <- pixelMetabolism(pixes, light, temp, pressure, lp1, lp2, er24_20, 
		k600, dt=1), regex=NA)

	meanPctDif <- function(x, y) mean(abs(x - y) / ((x+y)/2))
	tol <- 0.005

	# note; the tolerance is quite tight here, but some error should be expected
	# the comparison fits are from a bayesian model, so some error/uncertainty can be expected
	# in theory, they are comparisons from the deterministic part of the model, so should
	# be identical, but I expect randomness and also a bit of numerical error works its way in
	expect_lt(meanPctDif(pixmetab$DO[1,], compDO[1,]), tol)
	expect_lt(meanPctDif(pixmetab$DO[2,], compDO[2,]), tol)
	expect_lt(meanPctDif(pixmetab$GPP[1,], compGPP[1,]), tol)
	expect_lt(meanPctDif(pixmetab$GPP[2,], compGPP[2,]), tol)
	expect_lt(meanPctDif(pixmetab$ER[1,], compER[1,]), tol)
	expect_lt(meanPctDif(pixmetab$ER[2,], compER[2,]), tol)

})
