#' Calibrate a single station metabolism model
#' @param oxy Dissolved oxygen time series (g/L)
#' @param temp Water temperature time series (degrees); one per DO observation
#' @param light Light, (lux or W/m^2); one per DO observation
#' @param pressure Air pressure (hPa); single value or one per DO observation
#' @param depth Water depth, in meters
#' @param dt Time interval of observations (minutes)
#' @param nsamples The number of posterior samples to return
#' @param prior Priors for the function chosenl, see 'details'
#' @param ... Additional arguments to pass to the rstan
#' 
#' @details Implements a one-station stream metaoblism model; much code from Fuß et al 2017. 
#' 
#' This model requires that the observations consist of a single 24-hour period. Start time
#' may be arbitrary, as long as an entire 24-hours is covered. No checking is done for this, 
#' it is up to the user to ensure this is correct.
#' 
#' For a prior, you may specify a named list, with each item a vector of parameters to use for
#' the prior. Currently, the following parameters are supported:
#' * lp1: the log of the inverse of the slope of the photosynthesis-irradiance curve; default 
#' 		`c(0, 9)`.
#' * lp2: the log of the saturation term of the photosynthesis-irradiance curve; default 
#' 		`c(0, 9)`.
#' * k: gas transfer velocity, default c(0, 10)
#' * gpp: daily gross primary productivity, default c(0, 10)
#' * er: daily in-situ ecosystem respiration, default c(0, 10) (note first item must be <= 0)
#' 
#' Priors all use a half normal distribution. 
#' @references Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
#' 		 dissolved organic matter and nutrient regimes. *Freshw Biol* **62**:582–599. 
#' @return A fitted rstan model
#' @export
onestation <- function(oxy, temp, light, pressure, depth, dt, 
	nsamples = 1000, prior = list(), ...) {

	# basic input validation
	if(length(oxy) != length(temp) | length(oxy) != length(light))
		stop("oxy, light, and temp must have the same length")
	if(length(pressure) == 1)
		pressure = rep(pressure, length(oxy))
	if(length(oxy) != length(pressure))
		stop("length(pressure) must be 1 or equal to length(oxy)")

	if(! "lp1" %in% prior)
		prior$lp1 <- c(0,9)
	if(! "lp2" %in% prior)
		prior$lp2 <- c(0,9)
	if(! "k" %in% prior)
		prior$k = c(0, 10)
	if(! "gpp" %in% prior)
		prior$gpp = c(0, 10)
	if(! "er" %in% prior)
		prior$er = c(0, 10)

	data = list(
		nDO = length(oxy), 
		DO = oxy,
		temp = temp,
		light = light,
		pressure = pressure,
		depth = depth,
		delta_t = dt,
		lp1_pr = prior$lp1,
		lp2_pr = prior$lp2,
		k_pr = prior$k,
		gpp_pr = prior$gpp,
	 	er_pr = prior$er)

	file = system.file("stan/onestation.stan", package="NSmetabolism")
	modcode = rstan::stanc_builder(file = file)
	mod = rstan::stan_model(stanc_ret = modcode)
	calib = rstan::sampling(mod, data = data, iter = nsamples, ...)
	return(calib)
}

