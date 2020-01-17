

#' Transfornm parameters for the one station model
#' @param params The vector of parameters, see [oneStation_DOlogprob()], or a matrix of
#' 		samples (in which case the transformation is applied to app samples)
#' @param reverse Should the reverse transformation be applied?
#' @details This function transforms the parameters of the one station model from the
#' 	scale at which the parameters are calibrated to the scale of modelling. If
#'  `reverse == TRUE` the opposite transformation is performed.
#' @return The transformed parameters
#' @keywords internal
oneStation_parTransform <- function(params, reverse = FALSE) {
	if(reverse) {
		if(is.matrix(params)) {
			logP1 <- log(params[,'P1'])
			logP2 <- log(params[,'P2'])
			logk600 <- log(params[,'k600'])
			logMinusER24_20 <- log(-1 * params[,'ER24_20'])
			logsd <- log(params[,'sd'])
			result <- cbind(logP1, logP2, logk600, logMinusER24_20, logsd)
		} else {
			logP1 <- log(params['P1'])
			logP2 <- log(params['P2'])
			logk600 <- log(params['k600'])
			logMinusER24_20 <- log(-1 * params['ER24_20'])
			logsd <- log(params['sd'])
			result <- c(logP1 = logP1, logP2 = logP2, logk600 = logk600, 
				logMinusER24_20 = logMinusER24_20, logsd = logsd)
		}
	} else {
		if(is.matrix(params)) {
			P1 <- exp(params[,'logP1'])
			P2 <- exp(params[,'logP2'])
			k600 <- exp(params[,'logk600'])
			ER24_20 <- -1 * exp(params[,'logMinusER24_20'])
			sd <- exp(params[,'logsd'])
			result <- cbind(P1, P2, k600, ER24_20, sd)
		} else {
			P1 <- exp(params['logP1'])
			P2 <- exp(params['logP2'])
			k600 <- exp(params['logk600'])
			ER24_20 <- -1 * exp(params['logMinusER24_20'])
			sd <- exp(params['logsd'])
			result <- c(P1, P2, k600, ER24_20, sd)
			names(result) <- c("P1", "P2", "k600", "ER24_20", "sd")
		}
	}
	return(result)
}

#' Produce posterior simulations for DO curve and GPP
#' @param params Parameter posterior matrix returned (e.g.) from [DOCalibration()].
#' @param data Data list as required by [oneStation_DOPredict()].
#' @return list including predicted DO and GPP at each time step as well as the prediction times.
oneStation_DOSim <- function(params, data) {
	initial <- data$DO[1,1]
	times <- data$DO[,2]

	sims <- apply(params, 1, function(par) {
		c(oneStation_DOPredict(initial, times, par, data, dt = 1, gpp=TRUE))
	})
	list(DO = sims[1:length(times),], GPP = sims[(length(times)+1):nrow(sims),], times = times)
}
