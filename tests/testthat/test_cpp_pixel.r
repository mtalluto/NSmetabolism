context("Pixel functions")
library("NSmetabolism")
library("WatershedTools")

test_that("Pixel-based metabolism", {
	ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))

	## eventually
	## ws$light needs to be a thing

	## for now, generate a fake light curve for two pixels
	light <- matrix(c(rep(0,10), 1:10, rep(0,10), 11:20), nrow=2, byrow=TRUE)
	pixes <- as.data.frame(ws$data)[1:2,]
	lp1 <- c(5, 4)
	lp2 <- c(2, NA)

	expect_error(pixgpp <- pixelGPP(pixes, light, lp1, lp2), regex=NA)
	reggpp <- sapply(light[1,], computeGPP, lP1=lp1, lP2=lp2)
	expect_identical(reggpp, pixgpp[1,])
})
