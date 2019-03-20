context("One Station")
library("NSmetabolism")
library(lubridate)

disso <- read_minidot(system.file("testdata/metab/minidot", package="NSmetabolism"))
light <- read_hobo(system.file("testdata/metab/hobo.csv", package="NSmetabolism"))
## the elevation and depth at this site
elev <- 447.3237
z <- 0.722
# filter to only two days, april 23 and 24
startTime <- ymd_hms("2018-04-23 00:00:00", tz = tz(disso$timestamp[1]))
endTime <-  ymd_hms("2018-04-25 00:00:00", tz = tz(disso$timestamp[1]))
disso <- disso[disso$timestamp >= startTime & disso$timestamp <  endTime,]
light <- light[light$timestamp >= startTime & light$timestamp <  endTime,]

## approximate conversion for light
## from http://bccp.berkeley.edu/o/Academy/workshop08/08%20PDFs/Inv_Square_Law.pdf
light$light <- light$light * 0.0079

# add a time column in elapsed minutes
disso$minutes <- (as.integer(disso$timestamp) - as.integer(startTime)) / 60
light$minutes <- (as.integer(light$timestamp) - as.integer(startTime)) / 60

initParams <- c(P1=69911.063207, P2=220.062287, k600 = 2/(24*60), ER24_20=-5.386056, sd=1)

test_that("Pressure is computed right from elevation", {
	expect_equal(pressureFromElevation(0), 1)
	expect_equal(pressureFromElevation(500), 0.94, tol=0.01)
	expect_equal(pressureFromElevation(1000), 0.89, tol=0.01)
	expect_equal(pressureFromElevation(2000), 0.78, tol=0.01)
	expect_equal(pressureFromElevation(6000), 0.47, tol=0.01)
})

i <- which.max(light$light)
j <- which(disso$minutes > light$minutes[i] - 5 & disso$minutes < light$minutes[i] + 5)

test_that("dDOdt function working for one station", {
	DOdata <- list(
		PAR = approxfun(x = light$minutes, y = light$light, rule = 2),
		temp = approxfun(x = disso$minutes, y = disso$temperature, rule = 2),
		P = pressureFromElevation(elev),
		z = z)

	# dDOdt is negative at night when light == 0; the first point works
	expect_lt(ddo <- oneStation_dDOdt(disso$minutes[1], disso$DO[1], 
		initParams, DOdata), 0)
	# expect an absolute value less than the max in the data
	expect_lt(abs(ddo), max(disso$DO))
	# expect a positive value when light is at it's maximum
	expect_gt(ddo <- oneStation_dDOdt(disso$minutes[j], disso$DO[j], 
		initParams, DOdata), 0)
	# expect an absolute value less than the max in the data
	expect_lt(ddo, max(disso$DO))

})


test_that("DOpredict working for one station", {
	DOdata <- list(
		PAR = light[,c('light', 'minutes')],
		temp = disso[,c('temperature', 'minutes')],
		P = pressureFromElevation(elev),
		z = z)

	expect_error(doPr_euler <- oneStation_DOPredict(disso$DO[1], c(0, disso$minutes[1],   
		disso$minutes[j]), initParams, DOdata, dt = 0.1, method='euler'), regex=NA)
	expect_error(doPr_lsoda <- oneStation_DOPredict(disso$DO[1], c(0, disso$minutes[1],
		 disso$minutes[j]), initParams, DOdata, method='lsoda'), regex=NA)

	# quantitative predictions don't make a ton of sense, but we can try qualitative ones
	# similar to each other
	expect_lt(sum(abs(doPr_euler - doPr_lsoda)), 1e-2)
	# should be close to the initial value (within 5x for sure!)
	expect_true(all((doPr_euler / doPr_euler[2]) > 1/5 & (doPr_euler / doPr_euler[2]) < 5))

	# strictly nonnegative
	expect_true(all(doPr_euler >= 0))
	expect_true(all(doPr_lsoda >= 0))
})


