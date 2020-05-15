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



#' Interpolate light for metabolism models
#' 
#' Ideally, light data is present for all observations. However, sensor failures, logistical
#' failures, curious locals, etc often leave us with usable DO data, but not light. This
#' function uses spatial interpolation to try to salvage these data. Note that the results should
#' be used with caution, as the light record may be wrong.
#' 
#' This function is intended for use with sites spanning a single time period in a single
#' watershed (i.e., one field expedition). To use in multiple locations or over disconnected
#' time periods, please split your data first.
#' 
#' The interpolation proceeds in several steps. First, irradiance is computed (in W/m^2) using
#' [irradiance()] for every site at time intervals determined by `t_res`. Then,
#' a randomForest model is used to with the calibration data
#' 
#' The irradiance model is quite slow to run, so it is inefficient and unnecessary to run it for 
#' every desired output time. The parameter `t_res` controls how often the model is run; the 
#' irradiance will be computed at this interval (in minutes) starting from the first value in 
#' out_times and ending with the last. Linear interpolation will then be used to fill in the gaps.
#' 
#' @param data Data.frame with calibration light data, should contain columns named 'light' and
#'     'time'; the model will predict light for all observations where light is na
#' @param sites Spatial (sp or sf), locations corresponding to sites, site names
#'     names must match those in data
#' @param elev Raster digital elevation model for the study area
#' @param tz Integer, timezone of the study area, in hours offset from UTC
#' @param site_id Character, the column name where site names are stored;  same name
#' 		must be in data and sites
#' @param t_res Temporal resolution, in minutes, for the irradiance model. See 'details'
#' @param skip_irradiance Boolean; if TRUE, will skip the irradiance computation, in which case
#'      a column named 'irradiance' is required in data
#' @param ... Additional parameters to pass to [irradiance()]
#' 
#' @return A matrix, one row per output site and one column per output time
#' @export
interpolate_light = function(data, sites, elev, tz, site_id = "site", t_res = 30, 
	skip_irradiance = FALSE, ...) {

	if(!skip_irradiance) {
		# irradiance portion; compute estimates of W/m^2 for each site at each time
		irrad = .multisite_irradiance(data, sites, elev, tz, site_id, t_res, ...)
	} else {
		if(! "irradiance" %in% colnames(data))
			stop("If skip_irradiance is TRUE an irradiance column is required in data")

		# now convert the irradiance to whatever units original light data were in
		full_sun = .light_convert(data)

		# finally adjust for local conditions with interpolation
		.light_interp
	}
}

#' Multisite irradiance computation
#' @keywords internal
.multisite_irradiance = function(data, sites, elev, tz, site_id, t_res, ...) {
	ir_times = seq(min(data$time), max(data$time), by = paste(t_res, "min"))
	irr = irradiance(elev, sites, ir_times, timezone = tz, newSession = TRUE, ...)
	rownames(irr) = sites[[site_id]]



	irr_tall = as.data.frame.table(irr, stringsAsFactors = FALSE)
	colnames(irr_tall) = c(site_id, "time", "irradiance")
	irr_tall$time = as.integer(irr_tall$time)
	irr_pr = by(irr_tall, irr_tall[[site_id]], function(x) 
			data.frame(approx(as.integer(x$time), x$irradiance, 
			unique(as.integer(data$time)), rule=2)))
	irr_pr = data.table::rbindlist(irr_pr, idcol = site_id)
		
	# add irradiance into data
	irr_data = merge(data, irr_pr, by.x = c("time", site_id), by.y = c("x", site_id), all.x = TRUE)
	colnames(irr_data)[colnames(irr_data) == "y"] = "irradiance"
	calib = irr_data[light > 0 & data$irradiance > 0]
	calib$irradiance = log(calib$irradiance)
	calib$light = log(calib$light)

	mod = glm(light ~ irradiance, data = calib)
	preds = t(apply(irr, 1, function(x) {
		res = predict(mod, newdata=data.frame(irradiance=log(x)))
		res = exp(res)
		res[x == 0] = 0
		res
	}))

	preds = predict(mod, newdata = data.frame(irradiance = log(irr)))
}

#' convert light from modelled irradiance to a realistic predicted light
#' @keywords internal
.light_convert = function(data, sites, site_id) {

	calSites = sites[sites[[site_id]] %in% unique(calib[[site_id]]), ]
	pred = data[is.na(data$light), ]
	prSites = sites[sites[[site_id]] %in% unique(pred[[site_id]]), ]


	prdat = data[, "irradiance", drop = FALSE]
	prdat$irradiance = log(prdat$irradiance)
	light_predict = predict(mod, newdata = prdat)
	light_predict = exp(light_predict)
	light_predict[data$irradiance == 0] = 0

	## compute a local per-minute multiplier
	light_mult = data$light / light_predict
	light_mult[light_predict == 0] = 0

	# interpolate a multiplier for sites with missing data, then compute final light data
	
	# compute distance-based weights for each site
	dmat = st_distance(prSites, calSites)
	dwts = t(apply(dmat, 1, function(x) {x = 1/(x^2); x/sum(x)}))
	
	for(i in seq_len(nrow(calSites))) {
		dat = calib[calib[[site_id]] == sites[[site_id]][i], ]
		mod = glm(light ~ irradiance, data = calib)
	}


}

#' Adjust light for local conditions
#' @keywords internal
.light_interp(data, sites, light) {

}
