## data prep functions


#' Interpolate air pressure for metabolism models
#' 
#' In many cases we don't have air pressure precisely at the locations where other stream 
#' metabolism measurements are made. This function takes as input pressure readings from weather
#' stations, along with desired output locations and times. Output is matrix, with output sites
#' in the rows, output times in the columns.
#' 
#' `pressure` must be a data.frame with at least 3 columns, named:
#'      * 'station' (by default, or the value assigned to stn_id), station name
#' 		* 'pressure', pressure observations (hPa)
#' 		* 'time', a POSIXct column giving the times of measurements
#' 
#' `stations`: Spatial object, in addition to locations has at least 2 columns, named:
#'      * 'station' (by default, or the value assigned to stn_id), station name
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
#' @param out_times A list, with one element per out_site; each element is a vector of type 
#'		POSIXct giving times at which to compute pressure
#' @param parallel Boolean, should we use multiple threads to speed computation? Mac/Linux only
#' @param stn_id Character, the column name where pressure station names are stored in both
#' 		`pressure` and `stations`
#' @param site_id Character or numeric, the column in out_sites where site names are stored
#' @param ... Additional parameters to pass to [approx()].
#' 
#' @return A data frame with three columns: site name, time, and interpolated pressure
#' @export
interpolate_air_pressure = function(pressure, stations, out_sites, out_times, 
			parallel = TRUE, stn_id = "station", site_id = "site", ...) {
	if(!is(stations, "sf"))
		stations = sf::st_as_sf(stations)
	if(!is(out_sites, "sf"))
		out_sites = sf::st_as_sf(out_sites)
	.check_pressure_args(pressure, stations, out_sites, out_times)


	do_pressure = function(i) 
		.pr_site(pressure, stations, out_sites[i,], out_times[[i]], stn_id, site_id, ...)

	res = if(parallel && requireNamespace("parallel")) {
		parallel::mclapply(1:nrow(out_sites), do_pressure)
	} else {
		lapply(1:nrow(out_sites), do_pressure)
	}
	res = do.call(rbind, res)
	res
}


#' get a vector of pressure for a single site
#' @keywords internal
.pr_site = function(pr, st, si, times, stn_id, site_id, ...) {
	days = unique(lubridate::date(times))
	pr = pr[lubridate::date(pr$time) %in% days,]
	pr$elevation = st$elevation[match(pr[[stn_id]], st[[stn_id]])]

	# correct observed pressure to the elevation of the site
	pr$pressure_corr = pressureCorrection(pr$pressure, pr$elevation, si$elevation)
	pr_out = by(pr, pr[[stn_id]], function(x) 
		approx(as.integer(x$time), x$pressure_corr, as.integer(times), ...)$y)
	pr_out = do.call(rbind, pr_out)
	dmat = sf::st_distance(st, si)
	rownames(dmat) = st[[stn_id]]
	dmat = dmat[rownames(pr_out),]
	res = data.frame(time= times, pressure = idw_matrix(pr_out, dmat))
	res[[site_id]] = si[[site_id]]
	res
}

#' Validate input data for the air pressure function
#' @keywords internal
.check_pressure_args = function(pressure, stations, out_sites, out_times) {
	if(!all(c("pressure", 'time') %in% names(pressure)))
		stop("pressure must be a data.frame with columns 'time' and 'pressure'") 
	if(!all(c('elevation') %in% names(stations)))
		stop("stations must be a spatial object with 'elevation' column")
	if(length(out_times) != nrow(out_sites))
		stop("out_times must be a list with one element per site in out_sites")
	if(!is.POSIXct(pressure$time) | !all(sapply(out_times, is.POSIXct)))
		stop("pressure$time and all out_times elements must be of class POSIXct")
}

