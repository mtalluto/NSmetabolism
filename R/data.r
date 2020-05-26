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


#' Produce a summary from a list of one-station fits
#' 
#' If many fits are involved, memory usage can be extreme when providing `fits` instead of
#' `files.` It is recommended instead to save fits one at a time
#' 
#' If `diss_o` is provided, then the model root mean square error and the pearson correlation
#' (predicted vs calibration) will also be included
#' 
#' @param fits A list of onestation fits, required if `files` is not specified
#' @param files A vector or list of filenames with .rds files of onestation fits, required if
#' 		`fits` is not specified
#' @param diss_o Optional, a list of vectors with one element per element in fits and files. 
#' 		Contains the dissolved oxygen time series used to fit each model.
#' @export
osm_summarize = function(fits, files, diss_o) {
	if((missing(fits) & missing(files)) | (!missing(fits) & !missing(files)))
		stop("Exactly one of fits or files must be provided")
	
	if(missing(diss_o))
		diss_o = NA
	
	if(missing(fits)) {
		res = mapply(function(x, oxy) {
			fit = readRDS(x)
			.osm_summarize_single(fit, basename(x), oxy)
		}, files, diss_o, SIMPLIFY = FALSE)
	} else {
		res = mapply(.osm_summarize_single, fit = fits, name = names(fits), diss_o = diss_o,
					 SIMPLIFY = FALSE)
	}
	data.table::rbindlist(res)
}

#' Produces a variety of plots for one station models
#' @param fit A onestation model fit
#' @param oxy Dissolved oxygen data, required for `type=fit`
#' @param times Optional, vector of times at which oxy was measured
#' @param type The type of plot to produce
#' @return A ggplot object
osm_plot = function(fit, oxy, times, type = c("fit", "trace", "density")) {
	type = match.arg(type)
	if(is.na(nrow(fit))) {
		warning(modname, " contains no samples")
		return(NULL)	
	}
	if(type == "fit") {
		samples = as.matrix(fit, pars="DO_pr")
		if(missing(times))
			times = seq_along(oxy)
		do_pts = data.frame(time = times, oxy = oxy)
		pl_data = data.table::data.table(median = apply(samples, 2, median), 
										 lower = apply(samples, 2, quantile, 0.05),
										 upper = apply(samples, 2, quantile, 0.95), 
										 times = times)
		pl = ggplot(pl_data, aes(x = times, y = median)) + 
			geom_ribbon(aes(min = lower, max = upper), fill = '#a6cee3', alpha = 0.6) + 
			geom_line(col = '#1f78b4', size = 0.3) + 
			geom_point(data = do_pts, aes(x = times, y = oxy), size = 0.4, col = "#fb9a99") + 
			ylab(expression(Dissolved~Oxygen~(g/m^3))) + xlab("Time")
	} else {
		if(!requireNamespace("bayesplot", quitely=TRUE))
			stop("Package 'bayesplot' is required for this function, ",
				 "please install it and try again")
		params = c("lP1", "lP2", "ER24_20", "k600", "gpp", "er", "sigma")
		samples = as.array(fit, pars = params)
		if(type == 'trace') {
			pl = bayesplot::mcmc_trace(samples)
		} else if(type == "density") {
			pl = bayesplot::mcmc_dens(samples)
		} else {
			stop("Type ", type, " is not implemented")
		}
	}
	return(pl)
}



#' Summarize a single onestation fit
#' @param fit A onestation fit
#' @param name Model name
#' @param diss_o Optional, DO time series used to fit the model
#' @keywords internal
.osm_summarize_single = function(fit, name, diss_o) {
	if(is.na(nrow(fit))) {
		warning(name, " contains no samples")
		return(NULL)	
	}	
	
	params = c("lP1", "lP2", "ER24_20", "k600", "gpp", "er", "sigma")
	samples = as.matrix(fit, pars = params)
	samples[, c("lP1", "lP2")] = exp(samples[, c("lP1", "lP2")])
	colnames(samples)[colnames(samples) %in% c("lP1", "lP2")] = c("P1", "P2")
	
	if(!is.na(diss_o)) {
		do_pr = as.matrix(fit, pars="DO_pr")
		do_rmse = apply(do_pr, 1, .rmse, fit = diss_o)
		do_cor = apply(do_pr, 1, cor, y = diss_o, use = "complete.obs")
		samples = cbind(samples, rmse = do_rmse, pearson = do_cor)
	}
	
	data.table::data.table(
		file = name,
		parameter = colnames(samples),
		mean = apply(samples, 2, mean),
		median = apply(samples, 2, median),
		se = apply(samples, 2, sd),
		q0.05 = apply(samples, 2, quantile, 0.05, na.rm = TRUE),
		q0.95 = apply(samples, 2, quantile, 0.95, na.rm = TRUE))
}

#' Root mean square error
#' @keywords internal
.rmse = function(predict, fit) sqrt(mean((predict - fit)^2))


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

