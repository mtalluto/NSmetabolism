context("Metabolism Functions")
library("NSmetabolism")

test_that("Pressure is computed right from elevation", {
	expect_equal(pressureFromElevation(0), 1013.25, tol=0.01)
	expect_equal(pressureFromElevation(500), 954.6, tol=0.01)
	expect_equal(pressureFromElevation(1000), 898.75, tol=0.01)
	expect_equal(pressureFromElevation(2000), 794.95, tol=0.01)
	expect_equal(pressureFromElevation(6000), 471.81, tol=0.01)
})

test_that("Conversion to atm works", {
	expect_equal(hPaToAtm(1013.25), 1, tol = 0.01)
})

test_that("Pressure conversion", {
	pr_meas = c(1000, 1100)
	elev_meas = c(100, 200)
	elev_out = 1000
	
	# errors
	expect_error(pressureCorrection(pr_meas, elev_meas[1], elev_out), regex = "length")
	
	# sigle point should work fine
	expect_error(one <- pressureCorrection(pr_meas[1], elev_meas[1], elev_out), regex=NA)
	expect_equal(length(one), 1)
	expect_true(one > 875 & one < 900)
	
	# multiple pts
	expect_error(two <- pressureCorrection(pr_meas, elev_meas, elev_out), regex = NA)
	expect_equal(length(two), length(pr_meas))
	expect_true(all(two < pr_meas))
	
	# same elevation yields same pressure
	expect_error(same <- pressureCorrection(pr_meas[1], elev_meas[1], elev_meas[1]), regex=NA)
	expect_equal(same, pr_meas[1])
})