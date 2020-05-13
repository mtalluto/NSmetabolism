context("Priors")
library("NSmetabolism")

test_that("Generation of priors for k", {
	sl = 0.05
	v = 0.5
	q = 1
	d = 0.5

	# default case, all equations
	expect_error(res_all <- k_prior(sl, v, d, q, nsim=2), regex=NA)
	res_all = do.call(rbind, res_all)
	expect_equal(nrow(res_all), 7)
	expect_error(res_1 <- k_prior(sl, v, d, q, eq=1, nsim=2), regex=NA)
	expect_equal(res_all[1,1], res_1[,1])
	
	# degenerate cases
	expect_error(res_no_d <- k_prior(sl, v, discharge = q, nsim=2), regex=NA)
	expect_error(res_no_q <- k_prior(sl, v, depth = d, nsim=2), regex=NA)
	expect_error(res_no_v <- k_prior(sl, v = NA, discharge = q, depth = d, nsim=2), regex=NA)
	expect_error(res_sl_0 <- k_prior(slope = 0, velocity = v, discharge = q, depth = d, nsim=2), regex=NA)
	res_no_d = do.call(rbind, res_no_d)
	res_no_q = do.call(rbind, res_no_q)
	res_no_v = do.call(rbind, res_no_v)
	res_sl_0 = do.call(rbind, res_sl_0)
	
	expect_true(all(is.na(res_no_d[c(1,2,7),])))
	expect_true(all(is.finite(res_no_d[-c(1,2,7),])))
	expect_true(all(is.na(res_no_q[c(6,7),])))
	expect_true(all(is.finite(res_no_q[-c(6,7),])))
	expect_true(all(is.na(res_no_v)))
	expect_true(all(res_sl_0[-5,] == 0))
	
	# errors
	expect_error(k_prior(sl, v, d, q, eqn = 8, nsim=2), regex = "equation")
	
})