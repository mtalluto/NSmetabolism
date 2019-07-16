context("DO functions")
library("NSmetabolism")

test_that("Rearation flux produces sensible values", {
	T <- 20
	p <- 1025
	do <- 9
	k <- 2/(24*60)

	## basic input checking
	expect_error(compute_rf(T, p, -5, k))
	expect_error(compute_rf(T, p, do, -2))
	expect_error(compute_rf(T, -1, do, k))

	## DO below saturation
	expect_gt(compute_rf(T, p, do, k), 0)

	## DO above saturation
	expect_lt(compute_rf(T, p, 20, k), 0)
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
