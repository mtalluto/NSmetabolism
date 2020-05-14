## data prep functions


#' Interpolate air pressure for metabolism models
#' 
#' In many cases we don't have air pressure precisely at the locations where other stream 
#' metabolism measurements are made. This function takes as input pressure readings from weather
#' stations, along with desired output locations and times. Output is matrix, with output sites
#' in the rows, output times in the columns.
#' 
#' `pressure` must be a data.frame with at least 3 columns, named:
#' 		* 'station', a station name or id number, corresponding to the same column in `stations`
#' 		* 'pressure', pressure observations (hPa)
#' 		* 'time', a POSIXct column giving the times of measurements
#' 
#' `stations`: Spatial object, in addition to locations has at least 2 columns, named:
#' 		* 'station', a station name or id number, corresponding to the same column in `pressure`
#' 		* 'elevation', Elevation (m) of the station
#'
#' It is recommended that the input air pressure data represent a continuous record with no
#' substantial gaps in the time series. If this is not the case, interpolation may be unreliable.
#' 
#' Elevation is required for both control and output points so that corrections to air pressure
#' can be applied.
#'
#' @param pressure Data.frame required columns are 'station' 'pressure' (hPa), 'time' (see details)
#' @param stations Spatial (sp or sf), station locations where pressure was recorded (see details)
#' @param out_sites Spatial (sp or sf) object with output locations, 'elevation' is a required
#'  	column
#' 
#' @return A matrix, one row per output site and one column per output time
#' @export
interpolate_air_pressure = function(pressure, stations, out_sites, out_times) {
	if(!is(stations, "sf"))
		stations = sf::st_as_sf(stations)
	if(!is(out_sites, "sf"))
		out_sites = sf::st_as_sf(out_sites)
	.check_pressure_args(pressure, stations, out_sites, out_times)

	do_pressure = function(stn) 
		.pr_site(pressure, stations, out_sites[out_sites$station == stn,], out_times)

	res = if(requireNamespace("parallel")) {
		parallel::mclapply(out_sites$station, do_pressure)
	} else {
		lapply(out_sites$station, do_pressure)
	}
	res = do.call(rbind, res)
	if("station" %in% colnames(out_sites))
		rownames(res) = out_sites$station
	colnames(res) = as.character(out_times)
	res
}


#' get a vector of pressure for a single site
#' @keywords internal
.pr_site = function(pr, st, si, times) {
	days = unique(lubridate::date(times))
	pr = pr[lubridate::date(pr$time) %in% days,]
	pr$elevation = st$elevation[match(pr$station, st$station)]

	# correct observed pressure to the elevation of the site
	pr$pressure_corr = pressureCorrection(pr$pressure, pr$elevation, si$elevation)
	pr_out = by(pr, pr$station, function(x) 
		approx(as.integer(x$time), x$pressure_corr, as.integer(times))$y)
	pr_out = do.call(rbind, pr_out)
	dmat = sf::st_distance(st, si)
	rownames(dmat) = st$station
	dmat = dmat[rownames(pr_out),]
	idw_matrix(pr_out, dmat)
}

#' Validate input data for the air pressure function
#' @keywords internal
.check_pressure_args = function(pressure, stations, out_sites, out_times) {
	if(!all(c("pressure", 'station', 'time') %in% names(pressure)))
		stop("pressure must be a data.frame with columns 'station', 'time', and 'pressure'") 
	if(!all(c("station", 'elevation') %in% names(stations)))
		stop("stations must be a spatial object with columns 'elevation' and 'station'") 
	if(!is.POSIXct(pressure$time) || !is.POSIXct(out_times))
		stop("pressure$time and out_times must be of class POSIXct")

}

