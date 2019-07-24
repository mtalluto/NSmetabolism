context("DO functions")
library("NSmetabolism")


test_that("Computing input mass flux", {
	upDO <- c(10, 8)
	upQ <- c(15, 20)
	latQ <- 5
	latDO <- 0

	# input checking
	expect_error(computeInputDOFlux(-1 * upQ, upDO, latQ, latDO))
	expect_error(computeInputDOFlux(upQ, -1 * upDO, latQ, latDO))
	expect_error(computeInputDOFlux(upQ, upDO, -1 * latQ, latDO))
	expect_error(computeInputDOFlux(upQ, upDO, latQ, -5))
	expect_error(computeInputDOFlux(upQ, upDO[1], latQ, latDO))

	# test with 0, 1, and two upstream neighbors
	expect_error(computeInputDOFlux(upQ, upDO, latQ, latDO), regex=NA)
	expect_error(computeInputDOFlux(upQ[1], upDO[1], latQ, latDO), regex=NA)
	expect_error(computeInputDOFlux(numeric(0), numeric(0), latQ, latDO), regex=NA)

	#higher discharge and DO concentration means higher flux
	expect_lt(computeInputDOFlux(upQ[1], upDO[1], latQ, latDO), 
		computeInputDOFlux(upQ[1]*2, upDO[1], latQ, latDO))
	expect_lt(computeInputDOFlux(upQ[1], upDO[1], latQ, latDO), 
		computeInputDOFlux(upQ[1], upDO[1]*2, latQ, latDO))

})

test_that("Advection", {
	dat <- c(Q = 10, area = 20, dx = 5)
	expect_error(computeAdvection(10, 10, dat), regex=NA)
	expect_gt(computeAdvection(10, 10, dat), computeAdvection(0, 10, dat))

	# input checking
	expect_error(computeAdvection(-10, 10, dat))
	expect_error(computeAdvection(10, -10, dat))
	dat['area'] <- -1
	expect_error(computeAdvection(10, 10, dat))
	dat[c('area', 'Q')] <- c(20, -2)
	expect_error(computeAdvection(10, 10, dat))
	dat[c('Q', 'dx')] <- c(10, -5)
	expect_error(computeAdvection(10, 10, dat))

})

test_that("Rearation flux produces sensible values", {
	T <- 20
	p <- 1025
	do <- 9
	k <- 2/(24*60)

	# input checking
	expect_error(kT(T, -5))
	expect_error(kT(T, c(k, 2*k)))
	expect_error(osat(c(T, T*2), p))

	# reference value for kT
	expect_equal(kT(T, k), 0.001477129, tol=1e-5)

	# reference value for Osat
	CstarO <- 10.084
	T <- 15
	P <- 1013.2500
	Patm <- 1
	ref <- CstarO * Patm
	expect_equal(osat(T, P), ref, tol=1e-3)

	## basic input checking
	expect_error(computeRF(c(T, T*2), p, do, k))
	expect_error(computeRF(T, p, do, c(k, 2*k)))
	expect_error(computeRF(T, p, -5, k))
	expect_error(computeRF(T, p, do, -2))
	expect_error(computeRF(T, -1, do, k))

	## DO below saturation
	expect_gt(computeRF(T, p, do, k), 0)

	## DO above saturation
	expect_lt(computeRF(T, p, 20, k), 0)
})


test_that("GPP and ER", {
	expect_error(computeGPP(-10, 5, 2))
	# vector sizes match
	expect_error(computeGPP(c(100, 120), 5, 2))

	# no light, no GPP
	expect_equal(computeGPP(0, 5, 2), 0)
	expect_equal(computeGPP(0, 5, NA), 0)

	# linear curves are equal no matter how we compute
	expect_equal(computeGPP(1000, 5, -Inf), computeGPP(1000, 5, NA))

	# GPP is positive when there is light
	expect_gt(computeGPP(1000, 5, 2), 0)
	expect_gt(computeGPP(1000, 5, NA), 0)

	# Linear curve is strictly greater than saturating
	expect_gt(computeGPP(1000, 5, NA), computeGPP(1000, 5, 2))

	# saturating curve saturates
	expect_gt(computeGPP(1000, 5, 2) - computeGPP(1, 5, 2), 
		computeGPP(11000, 5, 2) - computeGPP(10001, 5, 2))

	# er24 must be negative
	expect_error(inSituER(15, 5))

	# at 20 degrees we should get the same value back
	expect_equal(inSituER(20, -5), -5, tol = 1e-2)

	# otherwise make sure ER is negative, and is 0 if er24 is 0
	expect_equal(inSituER(15, 0), 0)
	expect_lt(inSituER(15, -5), 0)

	# test using a vector of temps with one or more ERs
	expect_error(inSituER(c(5,10,15), c(-5, -6)))
	expect_error(inSituER(c(5,10,15), c(-5, -6, -10)), regex=NA)
	expect_error(inSituER(c(5,10,15), c(-5)), regex=NA)
})



# test_that("ddodt", {
# 	defaultParams <- c(lP1 = 10, lP2 = 5, k600 = 2/(24*60), ER24 = -10)
# 	defaultData <- c(Q = 20, area = 45, dx = 24, z = 1)
# 	defaultInputDOMass <- 9 * defaultData['Q']
# 	defaultDOPrev <- 8
# 	defaultlight <- 800
# 	defaultWaterTemp = 15
# 	defaultPressure = 1015

# 	# check for errors on invalid input
# 	pars <- defaultParams
# 	pars['ER24'] <- 5
# 	expect_error(dDOdt(pars, defaultData, defaultInputDOMass, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	data <- defaultData
# 	data['Q'] <- -10
# 	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	data <- defaultData
# 	data['area'] <- -10
# 	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	data <- defaultData
# 	data['dx'] <- -10
# 	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	data['z'] <- -10
# 	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	expect_error(dDOdt(defaultParams, defaultData, -10, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, -10, defaultlight,
# 		defaultWaterTemp, defaultPressure))
# 	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, defaultDOPrev, -10,
# 		defaultWaterTemp, defaultPressure))
# 	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, defaultDOPrev, defaultlight,
# 		defaultWaterTemp, -10))

# 	# no error on valid input
# 	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, defaultDOPrev, 
# 		defaultlight, defaultWaterTemp, defaultPressure), regex = NA)
# 	# no light and no input means loss of oxygen
# 	pars <- defaultParams
# 	pars['k600'] <- 0
# 	expect_lt(dDOdt(pars, defaultData, 0, defaultDOPrev, 0, defaultWaterTemp, defaultPressure), 0)
# 	# lots of light OR lots of input and no respiration and output means increasing oxygen
# 	pars['ER24'] <- 0
# 	expect_gt(dDOdt(pars, defaultData, 0, 0, 1000, defaultWaterTemp, 
# 		defaultPressure), 0)
# 	expect_gt(dDOdt(pars, defaultData, 50, 0, 0, defaultWaterTemp, defaultPressure), 0)

# })
