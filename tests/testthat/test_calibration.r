context("Calibration")
library("NSmetabolism")
library(lubridate)

disso <- read_minidot(system.file("testdata/metab/minidot", package="NSmetabolism"))
lightdat <- read_hobo(system.file("testdata/metab/hobo.csv", package="NSmetabolism"))
elev <- 447.3237
z <- 0.722
# filter to only two days, april 23 and 24
startTime <- ymd_hms("2018-04-23 00:00:00", tz = tz(disso$timestamp[1]))
endTime <-  ymd_hms("2018-04-25 00:00:00", tz = tz(disso$timestamp[1]))
disso <- disso[disso$timestamp >= startTime & disso$timestamp <  endTime,]
lightdat <- lightdat[lightdat$timestamp >= startTime & lightdat$timestamp <  endTime,]
teeZero <- disso$timestamp[which.min(disso$timestamp)]

## approximate conversion for light
## from http://bccp.berkeley.edu/o/Academy/workshop08/08%20PDFs/Inv_Square_Law.pdf
lightdat$light <- lightdat$light * 0.0079

# add a time column in elapsed minutes
disso$minutes <- (as.integer(disso$timestamp) - as.integer(teeZero)) / 60
lightdat$minutes <- (as.integer(lightdat$timestamp) - as.integer(teeZero)) / 60

#' 
#' `data` is a named list of (constant) data items, including the following:
#' * `DO` 2-column data frame; first column is dissolved oxygen, second is time of observation
#' * `PAR`  2-column data frame; first column is light, second is time of observation
#' * `temp` 2-column data frame; first column is temperature, second is time of observation
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters


initParams <- c(logP1=5, logP2=3, logMinusER24_20=1.5, logsd=1)

test_that("DO Calibration working for one station", {
	skip_on_cran()
	DOdata <- list(
		DO = as.data.frame(disso[,c("DO", "minutes")]),
		PAR = as.data.frame(lightdat[,c('light', 'minutes')]),
		temp = as.data.frame(disso[,c('temperature', 'minutes')]),
		P = pressureFromElevation(elev),
		z = z)
	ns <- 1000
	expect_error(samps <- DOCalibration(initParams, DOdata, nsamples = ns), regex=NA)
	expect_equal(nrow(samps), ns) 
	expect_equal(ncol(samps), length(initParams))
})


