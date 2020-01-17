#' Calibrate a dissolved oxygen model 
#' @param initial a named vector with one element per parameter giving the initial values 
#' @param data Dataset to be passed along to the model chosen
#' @param model The model to use; see 'details'
#' @param method Calibration method; either `laplace` approximation or `mcmc` (not implemented)
#' @param nsamples The number of posterior samples to return
#' @param prior Priors for the function chosenl, see 'details'
#' @param ... Additional arguments to pass to the sampler
#' 
#' @details The `model` parameter selects which model will be used; one of:
#' * onestation: The one-station model from Fuß et al 2017; params/data from
#' 		[oneStation_DOlogprob()]
#' * twostation: Not implemented
#' * nstation: Not implemented
#' 
#' For a prior, you may specify a named list, with each item a vector of parameters to use for
#' the prior. Currently, the following parameters are supported:
#' * lp1: the log of the inverse of the slope of the photosynthesis-irradiance curve; default 
#' 		`c(0, 10)`.
#' * lp2: the log of the saturation term of the photosynthesis-irradiance curve; default 
#' 		`c(0, 10)`.
#' @references Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
#' 		 dissolved organic matter and nutrient regimes. *Freshw Biol* **62**:582–599. 
#' @return A matrix of posterior samples of all parameters
#' @export
DOCalibration <- function(initial, data, model = c('onestation', 'twostation', 'nstation'), 
			method = c('laplace', 'stan'), nsamples = 1000, prior = list(), dt = 1, ...) {
	model <- match.arg(model)
	method <- match.arg(method)

	if(method == 'laplace') {
		stop("Laplace approximation is deprecated; use method = 'stan' instead")
		calib <- calibLaplace(initial, data, model, nsamples, prior)
	} else if(method == 'stan') {
		if(! "lp1" %in% prior)
			prior$lp1 <- c(0,9)
		if(! "lp2" %in% prior)
			prior$lp2 <- c(0,9)
		calib <- calibStan(data, nsamples, model, dt, prior, ...)
	} else {
		stop("Method ", method, " is not implemented yet")
	}

	return(calib)
}



#' keywords @internal
calibStan <- function(data, nsamples, model, dt, prior, ...) {
	if(!require("rstan"))
		stop("Method 'stan' requires the rstan package; please install it and try again")

	if(model == 'onestation') {
		file <- system.file("stan/onestation.stan", package="NSmetabolism")
		# if(is.list(data)) {
			stanDat <- data
			stanDat$dt = dt
			stanDat$lp1mu <- prior$lp1[1]
			stanDat$lp1sd <- prior$lp1[2]
			stanDat$lp2mu <- prior$lp2[1]
			stanDat$lp2sd <- prior$lp2[2]
	} else if(model == 'twostation') {
		file <- system.file("stan/twostation.stan", package="NSmetabolism")
		stanDat <- data
		stanDat$dt = dt
	} else {
		stop("Model ", model, " is not implemented yet")
	}

	modCode <- stanc_builder(file = file)
	mod <- stan_model(stanc_ret = modCode)
	fit <- sampling(mod, data = stanDat, iter = nsamples, ...)
	return(fit)
}
