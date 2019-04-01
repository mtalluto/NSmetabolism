#' Time derivative of DO concentration
#' 
#' @param t Time step for the ode integration
#' @param y Current state (dissolved oxygen concentration)
#' @param parms a named vector of parameters; see details
#' @param data a list of data for the ode; see details
#' @param return_list logical, should a list (suitable for deSolve) be returned?
#' @details `params` must be a vector including the following named items:
#' * `P1` $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
#' 		 light intensity
#' * `P2` $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
#' 		  GPP is linear with light intensity instead of saturating
#' * `k600` coefficient of gas exchange for a gas with a Schmidt number of 600 (see [kT()]).
#' * `ER24_20` daily ecosystem respiration rate, standardized at 20 degrees C
#' 
#' `data` is a named list of (constant) data items, including the following:
#' * `PAR`  a function that takes time as a parameter and returns light intensity $(W/m^2)$
#' * `temp` a function that takes time as a parameter and returns water temperature (degrees C)
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters
#' 
#' @references Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
#' 		 dissolved organic matter and nutrient regimes. *Freshw Biol* **62**:582–599. \newline
#'		 Van de Bogert MC et al. 2007. Assessing pelagic and benthic metabolism using free
#' 		 water measurements. *Limnol. Oceanogr. Methods* **5**:145–155. 
#' 
#' @return Time derivative of dissolved oxygen
#' @export
oneStation_dDOdt <- function(t, y, parms, data, return_list=FALSE)
{
	# gross primary production; Fuß et al eq 2
	GPP <- gpp(data$PAR(t), parms['P1'], parms['P2'])

	# rearation flux, Van de Bogert et al eqn 4; note that k600 is fixed
	RF <- kT(data$temp(t), parms['k600']) * (osat(data$temp(t), data$P) - y)

	# ecosystem respiration; Fuß et al eq 4
	ER <- (parms['ER24_20'] / (60*24)) * (1.045^(data$temp(t) - 20))

	val <- (GPP + ER + RF) / data$z
	if(return_list) {
		return(list(val))
	} else {
		return(val)
	}
}

#' Predict a time series of dissolved oxygen concentration
#' @param initial initial dissolved oxygen concentration
#' @param times The times at which to solve the system
#' @param params a named vector of parameters; see details
#' @param data a list of data for the ode; see details
#' @param dt Time step for integration; only used if `method == 'euler'`
#' @param method The integration method to use; default is euler
#' @param gpp Logical, should GPP at each time step be returned as well?
#' @details Light and temperature time series will be approximated using linear interpolation
#' 		at the desired time steps
#' 
#' `params` must be a vector including the following named items:
#' * `P1` $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
#' 		 light intensity
#' * `P2` $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
#' 		  GPP is linear with light intensity instead of saturating
#' * `k600` coefficient of gas exchange for a gas with a Schmidt number of 600; (see [kT()])
#' * `ER24_20` daily ecosystem respiration rate, standardized at 20 degrees C
#' 
#' `data` is a named list of (constant) data items, including the following:
#' * `PAR`  2-column data frame; first column is light, second is time of observation
#' * `temp` 2-column data frame; first column is temperature, second is time of observation
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters
#' @return Time series of dissolved oxygen concentrations
#' @export
oneStation_DOPredict <- function(initial, times, params, data, dt = 1, 
			method=c('euler', 'lsoda'), gpp = FALSE) {
	method <- match.arg(method)
	if(method == "lsoda" && !(requireNamespace("deSolve")))
		stop("Package deSolve is required for method lsoda; please install it and try again")
	if(times[1] != 0)
		stop('first time step must be 0')
	if(method == 'euler' && sum(times %% dt) != 0)
		stop('dt must divide evenly into all values in times')

	dDOData <- list(
		PAR = approxfun(x = data$PAR[,2], y = data$PAR[,1], rule = 2),
		temp = approxfun(x = data$temp[,2], y = data$temp[,1], rule = 2),
		P = data$P,
		z = data$z)

	if(method == 'euler') {
		state <- numeric(length(times))
		state[1] <- statet <- initial
		t <- times[1]
		for(i in 2:length(state)) {
			while(t < times[i]) {
				statet <- statet + dt * oneStation_dDOdt(t, statet, params, dDOData)
				t <- t + dt
			}
			state[i] <- statet
		}
	} else {
		state <- deSolve::ode(initial, times = times, func = oneStation_dDOdt, parms = params, 
			data = dDOData, return_list=TRUE)
		state <- state[,2]
	}
	if(gpp) {
		state <- cbind(DO = state, GPP = gpp(dDOData$PAR(times), params['P1'], params['P2']))
	}
	return(state)
}

#' Log posterior probability for the DO curve given empirical data
#' @param params A named vector of model parameters; see details
#' @param data A list of model data, see details
#' @param prior A list of prior hyperparameters; see details
#' @param ... Additional parameters to be passed to [DO_predict()]
#' @details #' `params` must be a vector including the following named items:
#' * `logP1` $(W min g^{-1} O_2)$ log(inverse of the slope of a photosynthesis–irradiance curve)
#' 		 at low light intensity
#' * `logP2` $(m^2 min g^{-1} O_2)$ log(inverse maximum photosynthesis rate)
#' * `logk600` coefficient of gas exchange for a gas with a Schmidt number of 600; (see [kT()])
#' * `logMinusER24_20` log(-1 * daily ecosystem respiration rate), standardized at 20 degrees C
#' * `logsd` log(Error standard deviation)
#' 
#' `data` is a named list of (constant) data items, including the following:
#' * `DO` 2-column data frame; first column is dissolved oxygen, second is time of observation
#' * `PAR`  2-column data frame; first column is light, second is time of observation
#' * `temp` 2-column data frame; first column is temperature, second is time of observation
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters
#' 
#' `prior` is a list of hyperparameters for the priors on parameters listed under `params`.
#' All priors are normal. Thus, each list item should be named (following the names from 
#' `params`), and each item should be a vector of length two giving the mean and standard
#' deviation of the prior for that parameter. Note that each parameter is optional 
#' (as is the entire list); missing items will get a minimally informative normal(0,10)
#'  distribution.
#' @return The unnormalized log probability of the model given the data
#' @export
oneStation_DOlogprob <- function(params, data, prior = list(), ...) {
	# check for each parameter, assign default priors if not specified
	trParNames <- c('logP1', 'logP2', 'logk600', 'logMinusER24_20', 'logsd')
	parNames <- c('P1', 'P2', 'k600', 'ER24_20', 'sd')
	for(nm in trParNames) {
		if(! nm %in% names(prior))
			prior[[nm]] <- c(0,10)
	}

	initial <- data$DO[1,1]
	times <- data$DO[,2]

	# parameters are log transformed to make sampling nicer
	doParams <- oneStation_parTransform(params)
	predicted <- oneStation_DOPredict(initial, times, doParams, data)

	logprob <- sum(dnorm(data$DO[,1], predicted, doParams['sd'], log = TRUE)) + 
		dnorm(params['logP1'], prior[['logP1']][1], prior[['logP1']][2], log = TRUE) + 
		dnorm(params['logP2'], prior[['logP2']][1], prior[['logP2']][2], log = TRUE) + 
		dnorm(params['logk600'], prior[['logk600']][1], prior[['logk600']][2], log = TRUE) + 
		dnorm(params['logMinusER24_20'], prior[['logMinusER24_20']][1], prior[['logMinusER24_20']][2], log = TRUE) + 
		dnorm(params['logsd'], prior[['logsd']][1], prior[['logsd']][2], log = TRUE)

	return(logprob)
}


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
