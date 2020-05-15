context("IDW")
library("NSmetabolism")
library(lubridate)
library(sf)
test_that("Pressure interpolation", {
	pressure = data.frame(pressure = c(1000, 1100, 900, 950), station = c(1,1,2,2), 
			time = ymd_hm(c("2020_01-01 10:00", "2020_01-01 10:10", "2020_01-01 10:00", "2020_01-01 10:10")))
	stations = data.frame(station = c(1, 2), x = c(1, 2), y = c(1, 2), elevation = c(100, 200))
	stations = st_as_sf(stations, coords=c('x', 'y'))
	
	out_sites = data.frame(site = c('o1', 'o2'), x = c(0, 3), y = c(0, 3), elevation = c(150, 175))
	out_sites = st_as_sf(out_sites, coords = c('x', 'y'))
	out_times =  list(ymd_hm(c("2020_01-01 10:05", "2020_01-01 10:06", "2020_01-01 10:07")))
	out_times[[2]] = out_times[[1]]

	## check error conditions
	expect_error(interpolate_air_pressure(pressure[, -1], stations, out_sites, out_times), regex='pressure')
	#expect_error(interpolate_air_pressure(pressure[, -2], stations, out_sites, out_times), regex='station')
	expect_error(interpolate_air_pressure(pressure[, -3], stations, out_sites, out_times), regex='time')
	#expect_error(interpolate_air_pressure(pressure, stations[, -1], out_sites, out_times), regex='station')
	expect_error(interpolate_air_pressure(pressure, stations[, -2], out_sites, out_times), regex='elevation')
	expect_error(interpolate_air_pressure(pressure, st_drop_geometry(stations), out_sites, out_times), regex='geometry')
	expect_error(interpolate_air_pressure(pressure, stations, st_drop_geometry(out_sites), out_times), regex='geometry')
	expect_error(interpolate_air_pressure(pressure, stations, out_sites, lapply(out_times, as.integer)), regex='POSIXct')
	
	## only basic validation is done here
	expect_error(test <- interpolate_air_pressure(pressure, stations, out_sites, out_times, parallel=FALSE) , regex=NA)
	expect_equal(nrow(test), sum(sapply(out_times, length)))	
	expect_equal(length(unique(test$site)), nrow(out_sites))
})
