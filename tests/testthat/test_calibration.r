context("Calibration")
library("NSmetabolism")
library(rstan)
# library(lubridate)

test_that("One station calibration working", {
	skip_on_cran()
	dat = readRDS(system.file("testdata/do_calib.rds", package="NSmetabolism"))
	expect_error(mod <- onestation(dat$DO, dat$temp, dat$light, dat$pressure, dat$z, dat$delta_t, nsamples = 1000), regex=NA)
	pars = as.matrix(mod)
	expect_true(all(pars[,'er'] <= 0))
	expect_true(all(pars[, 'gpp'] >= 0))
})

