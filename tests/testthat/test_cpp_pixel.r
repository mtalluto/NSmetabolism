context("Pixel functions")
library("NSmetabolism")
library("WatershedTools")

test_that("Pixel-based metabolism", {
	ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))

	## eventually
	## ws$light and ws$temp needs to be a thing

	## for now, use light and other data from the db
	light <- matrix(c(rep(0,10), 1:10, rep(0,10), 11:20), nrow=2, byrow=TRUE)
	temp <- matrix(c(6:15, 16:7, rep(6, 20)), nrow=2, byrow=TRUE)
	pixes <- as.data.frame(ws$data)[1:2,]
	lp1 <- c(5, 4)
	lp2 <- c(2, NA)
	er24_20 <- c(-5, -8)


	expect_error(pixgpp <- pixelMetabolism(pixes, light, temp, pressure, do_init, 
		lp1, lp2, er24_20, k600, "gp"), 
		class="std::range_error")

	expect_error(pixgpp <- pixelMetabolism(pixes, light, temp, pressure, do_init, 
		lp1, lp2, er24_20, k600, "gpp"), regex=NA)
	reggpp <- sapply(light[1,], computeGPP, lP1=lp1[1], lP2=lp2[1])
	expect_identical(reggpp, pixgpp[1,])

	expect_error(pixer <- pixelMetabolism(pixes, light, temp, pressure, do_init, 
		lp1, lp2, er24_20, k600, 'er'), regex=NA)
	reger <- sapply(temp[1,], inSituER, ER24_20=er24_20[1])
	expect_identical(reger, pixer[1,])




})
