#' Calibrate a dissolved oxygen model 
#' @param initial a named vector with one element per parameter giving the initial values 
#' @param data Dataset to be passed along to the model chosen
#' @param model The model to use; see 'details'
#' @param method Calibration method; either `laplace` approximation or `mcmc` (not implemented)
#' @param nsamples The number of posterior samples to return
#' @param prior Priors for the function chosen
#' @param ... Additional arguments to pass to the sampler
#' 
#' @details The `model` parameter selects which model will be used; one of:
#' * onestation: The one-station model from Fuß et al 2017; params/data from
#' 		[oneStation_DOlogprob()]
#' * twostation: Not implemented
#' * nstation: Not implemented
#' @references Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
#' 		 dissolved organic matter and nutrient regimes. *Freshw Biol* **62**:582–599. 
#' @return A matrix of posterior samples of all parameters
#' @export
DOCalibration <- function(initial, data, model = c('onestation', 'twostation', 'nstation'), 
			method = c('laplace', 'stan'), nsamples = 1000, prior = list(), dt = 1, ...) {
	model <- match.arg(model)
	method <- match.arg(method)

	if(method == 'laplace') {
		calib <- calibLaplace(initial, data, model, nsamples, prior)
	} else if(method == 'stan') {
		calib <- calibStan(data, nsamples, model, dt, ...)
	} else {
		stop("Method ", method, " is not implemented yet")
	}

	return(calib)
}


#' @keywords internal
calibLaplace <- function(initial, data, model, nsamples, prior) {
	if(!requireNamespace("mvtnorm"))
		stop("Method 'laplace' requires the mvtnorm package; please install it and try again")

	if(model == 'onestation') {
		logProb <- oneStation_DOlogprob
		parTransform <- oneStation_parTransform 
		simFct <- oneStation_DOSim
	} else {
		stop("Model ", model, " is not implemented yet")
	}

	fit <- optim(initial, logProb, control = list(fnscale = -1), hessian = TRUE, data = data, 
		prior = prior)
	vcvMat <- solve(-fit$hessian)
	
	# check for positive semidefiniteness; if not, find the nearest positiveSD matrix
	eigenvalues <- eigen(vcvMat, only.values = TRUE)$values
	eigenvalues[abs(eigenvalues) < 1e-8] <- 0 
	if(any(eigenvalues <= 0)) {
		if(requireNamespace("Matrix", quietly=TRUE)) {
			vcvMat <- as.matrix(Matrix::nearPD(vcvMat)$mat)
		} else {
			warning("The parameter variance-covariance matrix is not positive semidefinite. The samples may not reflect the posterior distribution and should be considered approximate. Please install the 'Matrix' package for a better estimate.")
		}
	}

	samples <- mvtnorm::rmvnorm(nsamples, fit$par, vcvMat)
	samples <- parTransform(samples)
	sims <- simFct(samples, data)

	list(params = samples, DO = sims$DO, GPP = sims$GPP, time = sims$times, fit = fit)
}

#' keywords @internal
calibStan <- function(data, nsamples, model, dt, ...) {
	if(!require("rstan"))
		stop("Method 'stan' requires the rstan package; please install it and try again")

	if(model == 'onestation') {
		file <- system.file("stan/onestation.stan", package="NSmetabolism")
		# if(is.list(data)) {
			stanDat <- data
			stanDat$dt = dt
		# } else {
		# 	stanDat <- list(
		# 		nDO = nrow(data$DO),
		# 		timesDO = data$DO$minutes,
		# 		dt = dt,
		# 		z = data$z,
		# 		pressure = data$P,
		# 		DO = data$DO$DO,
		# 		slope = data$slope,
		# 		velocity = data$velocity
		# 	)		
		# 	parFun <- approxfun(x = data$PAR[,2], y = data$PAR[,1], rule = 2)
		# 	tempFun = approxfun(x = data$temp[,2], y = data$temp[,1], rule = 2)

		# 	stanDat$timesInt <- seq(stanDat$timesDO[1], stanDat$timesDO[length(stanDat$timesDO)], dt)
		# 	stanDat$nTime <- length(stanDat$timesInt)
		# 	stanDat$PAR <- parFun(stanDat$timesInt)
		# 	stanDat$temp <- tempFun(stanDat$timesInt)
		# 	if(any(stanDat$timesDO != as.integer(stanDat$timesDO)) | dt != as.integer(dt))
		# 		stop("Non-integer times are not supported; please ensure that dt and all times are integers")

		# 	if(!all(stanDat$timesDO %in% stanDat$timesInt))
		# 		stop("Please choose dt such that seq(times[1], times[length(times)], dt) includes all entries in times")
	} else if(model == 'twostation') {
		file <- system.file("stan/twostation.stan", package="NSmetabolism")
		stanDat <- data
		stanDat$dt = dt
	} else {
		stop("Model ", model, " is not implemented yet")
	}


	fit <- rstan::stan(file, data = stanDat, iter=nsamples, ...)
	return(fit)
}
