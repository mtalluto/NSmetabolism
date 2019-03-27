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
## https://physics.stackexchange.com/questions/135618/rm-lux-and-w-m2-relationship
lightdat$light <- lightdat$light * 0.0079

# add a time column in elapsed minutes
disso$minutes <- (as.integer(disso$timestamp) - as.integer(teeZero)) / 60
lightdat$minutes <- (as.integer(lightdat$timestamp) - as.integer(teeZero)) / 60

initParams <- c(logP1=5, logP2=3, logk600 = -6.6, logMinusER24_20=1.5, logsd=1)

test_that("DO Calibration working for one station", {
	skip_on_cran()
	DOdata <- list(
		DO = as.data.frame(disso[,c("DO", "minutes")]),
		PAR = as.data.frame(lightdat[,c('light', 'minutes')]),
		temp = as.data.frame(disso[,c('temperature', 'minutes')]),
		P = pressureFromElevation(elev),
		z = z)
	ns <- 1000
	expect_error(calib <- DOCalibration(initParams, DOdata, nsamples = ns), regex=NA)
	expect_equal(nrow(calib$params), ns) 
	expect_equal(ncol(calib$params), length(initParams))
})


# quartz()
# par(mfrow=c(2,1))
# plot(DOdata$DO$minutes, DOdata$DO$DO, pch=16, xlab='time', ylab='DO', cex=0.5, bty='n')
# lines(calib$time, rowMeans(calib$DO), col='#993333')
# polygon(c(calib$time, rev(calib$time)), c(apply(calib$DO, 1, quantile, 0.05), rev(apply(calib$DO, 1, quantile, 0.95))), col='#99333366', border=NA)
# par(new=TRUE)
# plot(DOdata$temp$minutes, DOdata$temp$temperature, col='#3333aa', axes=FALSE, bty='n', type='l', xlab='', ylab='')
# axis(side=4)
# mtext("Temperature", side=4)


# plot(DOdata$PAR$minutes, DOdata$PAR$light, pch=16, xlab='time', ylab='PAR', cex=0.5, bty='n')
# par(new=TRUE)
# plot(calib$time, rowMeans(calib$GPP), col='#993333', axes=FALSE, bty='n', type='l', xlab='', ylab='')
# polygon(c(calib$time, rev(calib$time)), c(apply(calib$GPP, 1, quantile, 0.05), rev(apply(calib$GPP, 1, quantile, 0.95))), col='#99333366', border=NA)
# par(new=TRUE)
# plot(calib$time, rowMeans(calib$GPP), col='#993333', axes=FALSE, bty='n', type='l', xlab='', ylab='')
