

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
