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
