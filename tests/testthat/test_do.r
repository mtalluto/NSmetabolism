context("DO functions")
library("NSmetabolism")

test_that("Rearation flux produces sensible values", {
	T <- 20
	p <- 1025
	do <- 9
	k <- 2/(24*60)

	## basic input checking
	expect_error(computeRF(T, p, -5, k))
	expect_error(computeRF(T, p, do, -2))
	expect_error(computeRF(T, -1, do, k))

	## DO below saturation
	expect_gt(computeRF(T, p, do, k), 0)

	## DO above saturation
	expect_lt(computeRF(T, p, 20, k), 0)
})


test_that("GPP and ER", {
	# no light, no GPP
	expect_equal(computeGPP(0, 5, 2), 0)
	expect_equal(computeGPP_linear(0, 5), 0)

	# linear curves are equal no matter how we compute
	expect_equal(computeGPP(1000, 5, -Inf), computeGPP_linear(1000, 5))

	# GPP is positive when there is light
	expect_gt(computeGPP(1000, 5, 2), 0)
	expect_gt(computeGPP_linear(1000, 5), 0)

	# Linear curve is strictly greater than saturating
	expect_gt(computeGPP_linear(1000, 5), computeGPP(1000, 5, 2))

	# saturating curve saturates
	expect_gt(computeGPP(1000, 5, 2) - computeGPP(1, 5, 2), 
		computeGPP(11000, 5, 2) - computeGPP(10001, 5, 2))

	# er24 must be negative
	expect_error(computeER(15, 5))

	# at 20 degrees we should get the same value back
	expect_equal(computeER(20, -5)*24*60, -5, tol = 1e-2)

	# otherwise make sure ER is negative, and is 0 if er24 is 0
	expect_equal(computeER(15, 0), 0)
	expect_lt(computeER(15, -5), 0)
})

test_that("Advection", {
	# input checking
	expect_error(computeAdvection(-10, 10, 20, 10, 5))
	expect_error(computeAdvection(10, -10, 20, 10, 5))
	expect_error(computeAdvection(10, 10, -20, 10, 5))
	expect_error(computeAdvection(10, 10, 20, -10, 5))
	expect_error(computeAdvection(10, 10, 20, 10, -5))
	expect_error(computeAdvection(10, 10, 20, 10, 5), regex=NA)

	expect_gt(computeAdvection(10, 10, 1, 1, 1), computeAdvection(0, 10, 1, 1, 1))
})

test_that("ddodt", {
	defaultParams <- c(lP1 = 10, lP2 = 5, k600 = 2/(24*60), ER24 = -10)
	defaultData <- c(Q = 20, area = 45, dx = 24, z = 1)
	defaultInputDOMass <- 9 * defaultData['Q']
	defaultDOPrev <- 8
	defaultlight <- 800
	defaultWaterTemp = 15
	defaultPressure = 1015

	# check for errors on invalid input
	pars <- defaultParams
	pars['ER24'] <- 5
	expect_error(dDOdt(pars, defaultData, defaultInputDOMass, defaultDOPrev, defaultlight,
		defaultWaterTemp, defaultPressure))
	data <- defaultData
	data['Q'] <- -10
	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
		defaultWaterTemp, defaultPressure))
	data <- defaultData
	data['area'] <- -10
	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
		defaultWaterTemp, defaultPressure))
	data <- defaultData
	data['dx'] <- -10
	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
		defaultWaterTemp, defaultPressure))
	data['z'] <- -10
	expect_error(dDOdt(defaultParams, data, defaultInputDOMass, defaultDOPrev, defaultlight,
		defaultWaterTemp, defaultPressure))
	expect_error(dDOdt(defaultParams, defaultData, -10, defaultDOPrev, defaultlight,
		defaultWaterTemp, defaultPressure))
	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, -10, defaultlight,
		defaultWaterTemp, defaultPressure))
	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, defaultDOPrev, -10,
		defaultWaterTemp, defaultPressure))
	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, defaultDOPrev, defaultlight,
		defaultWaterTemp, -10))

	# no error on valid input
	expect_error(dDOdt(defaultParams, defaultData, defaultInputDOMass, defaultDOPrev, 
		defaultlight, defaultWaterTemp, defaultPressure), regex = NA)
	# no light and no input means loss of oxygen
	pars <- defaultParams
	pars['k600'] <- 0
	expect_lt(dDOdt(pars, defaultData, 0, defaultDOPrev, 0, defaultWaterTemp, defaultPressure), 0)
	# lots of light OR lots of input and no respiration and output means increasing oxygen
	pars['ER24'] <- 0
	expect_gt(dDOdt(pars, defaultData, 0, 0, 1000, defaultWaterTemp, 
		defaultPressure), 0)
	expect_gt(dDOdt(pars, defaultData, 50, 0, 0, defaultWaterTemp, defaultPressure), 0)

})