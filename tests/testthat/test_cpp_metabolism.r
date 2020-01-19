context("Metabolism functions")
library("NSmetabolism")

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
