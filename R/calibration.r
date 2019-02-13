#' Calibrate a dissolved oxygen model 
#' @param initial a named vector with one element per parameter giving the initial values 
#' @param data Dataset to be passed along to the model chosen
#' @param model The model to use; see 'details'
#' @param method Calibration method; either `laplace` approximation or `mcmc` (not implemented)
#' @param nsamples The number of posterior samples to return
#' 
#' @details The `model` parameter selects which model will be used; one of:
#' * onestation: The one-station model from Fuß et al 2017; params/data from
#' 		[oneStation_DOlogprob()]
#' * twostation: Not implemented
#' * nstation: Not implemented
#' @references Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
#' 		 dissolved organic matter and nutrient regimes. *Freshw Biol* **62**:582–599. 
#' @return A matrix of posterior samples of all parameters
DOCalibration <- function(initial, data, model = c('onestation', 'twostation', 'nstation'), 
			method = c('laplace', 'mcmc'), nsamples = 1000) {
	model <- match.arg(model)
	method <- match.arg(method)

	if(model == 'onestation') {
		logProb <- oneStation_DOlogprob
		parTransform <- oneStation_parTransform 
	} else {
		stop("Model ", model, " is not implemented yet")
	}

	if(method == 'laplace') {
		if(!requireNamespace("mvtnorm"))
			stop("Method 'laplace' requires the mvtnorm package; please install it and try again")
		fit <- optim(initial, logProb, control = list(fnscale = -1), hessian = TRUE, data = data)
		vcvMat <- solve(-fit$hessian)
		samples <- mvtnorm::rmvnorm(nsamples, fit$par, vcvMat)
		samples <- parTransform(samples)
	} else {
		stop("Method ", method, " is not implemented yet")
	}

	return(samples)
}
