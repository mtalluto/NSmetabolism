#' Time derivative of DO concentration
#' 
#' @param t Time step for the ode integration
#' @param y Current state (dissolved oxygen concentration) at both sites
#' @param parms a vector of parameters; see details
#' @param data a list of data for the ode; see details
#' @param return_list logical, should a list (suitable for deSolve) be returned?
#' @details `params` must be a vector. The indices of the vector must map to model parameters
#'       as follows; (where parameters with length two are fit per site):
#' * [1:2] -- `P1` $(W min g^{-1} O_2)$ -- inverse of the slope of a photosynthesis–irradiance
#' 			curve at low light intensity
#' * [3:4] -- `P2` $(m^2 min g^{-1} O_2)$ -- inverse maximum photosynthesis rate; can be zero 
#' 			to assume GPP is linear with light intensity instead of saturating
#' * [5:6] -- `k600` -- coefficient of gas exchange for a gas with a Schmidt number of 600
#' * [7:8] -- `ER24_20` -- daily ecosystem respiration rate, standardized at 20 degrees C
#' 
#' `data` is a list of length 2 (one per site); each item is a list of (constant) data items, 
#' 				including the following:
#' * `PAR`  a function that takes time as a parameter and returns light intensity $(W/m^2)$
#' * `temp` a function that takes time as a parameter and returns water temperature (degrees C)
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters
#' 
#' @references Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
#' 		 dissolved organic matter and nutrient regimes. *Freshw Biol* **62**:582–599. \\
#'		 Van de Bogert MC et al. 2007. Assessing pelagic and benthic metabolism using free
#' 		 water measurements. *Limnol. Oceanogr.: Methods* **5**:145–155. \\
#' 		 Uehlinger, U., et al. 2000. Variability of photosynthesis‐irradiance curves and 
#' 		 ecosystem respiration in a small river. *Freshw Biol* **44**:493–507.
#' 
#' @return Time derivative of dissolved oxygen
#' @export
twoStation_dDOdt <- function(t, y, parms, data, return_list=FALSE)
{
	nms <- c('P1', 'P2', 'k600', 'ER24_20')
	parlist <- list(parms[c(1,3,5,7)], parms[c(2,4,6,8)])
	names(parlist[[1]]) <- names(parlist[[2]]) <- nms

	val <- mapply(function(pa, da) oneStation_dDOdt(t, y, pa, da, return_list = FALSE), 
		pa = parlist, da = data, SIMPLIFY = !return_list)
	return(val)
}

#' Predict a time series of dissolved oxygen concentration
#' @param ws A [WatershedTools::Watershed()] object, describing ONLY the two sites
#' @param initial initial dissolved oxygen concentration in each station
#' @param lateral Vector  of length two giving the DO concentration of lateral input 
#' 		(useful) if discharge varies significantly between the sites, defaults to 0
#' @param times The times at which to solve the system
#' @param params a named vector of parameters; see details
#' @param data a list of data for the ode; see details
#' @param dt Time step for integration; only used if `method == 'euler'`
#' @param method The integration method to use; default is euler

#' @details Light and temperature time series will be approximated using linear interpolation
#' 		at the desired time steps
#' 
#' `params` must be a vector. The indices of the vector must map to model parameters
#'       as follows; (where parameters with length two are fit per site):
#' * [1:2] -- `P1` $(W min g^{-1} O_2)$ -- inverse of the slope of a photosynthesis–irradiance
#' 			curve at low light intensity
#' * [3:4] -- `P2` $(m^2 min g^{-1} O_2)$ -- inverse maximum photosynthesis rate; can be zero 
#' 			to assume GPP is linear with light intensity instead of saturating
#' * [5:6] -- `k600` -- coefficient of gas exchange for a gas with a Schmidt number of 600
#' * [7:8] -- `ER24_20` -- daily ecosystem respiration rate, standardized at 20 degrees C
#' 
#' `data` is a list of length 2 (one per site); each item is a list of (constant) data items, 
#' 				including the following:
#' * `PAR`  2-column data frame; first column is light, second is time of observation
#' * `temp` 2-column data frame; first column is temperature, second is time of observation
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters
#' @return Time series of dissolved oxygen concentrations
#' @export
twoStation_DOPredict <- function(ws, initial, times, params, data, dt = 1, 
			method=c('euler', 'lsoda')) {
	method <- match.arg(method)
	if(method == "lsoda" && !(requireNamespace("deSolve")))
		stop("Package deSolve is required for method lsoda; please install it and try again")
	if(times[1] != 0)
		stop('first time step must be 0')
	if(method == 'euler' && sum(times %% dt) != 0)
		stop('dt must divide evenly into all values in times')

	dDOData <- lapply(data, function(x) {
		list(PAR = approxfun(x = x$PAR[,2], y = x$PAR[,1], rule = 2),
			temp = approxfun(x = x$temp[,2], y = x$temp[,1], rule = 2),
			P = x$P,
			z = x$z)
	})

	ddodtParams <- list(parms = params, data = dDOData)
	state <- WatershedTools::transport(ws, initial, lateral, times, method = method, dt = dt, 
			rxn = twoStation_dDOdt, rxnParams = ddodtParams)

	return(state)
}


#' Log posterior probability for the DO curve given empirical data
#' @param params A named vector of model parameters; see details
#' @param data A list of model data, see details
#' @param prior A list of prior hyperparameters; see details
#' @param ... Additional parameters to be passed to [DO_predict()]
#' @details `params` must be a vector. The indices of the vector must map to model parameters
#'       as follows; (where parameters with length two are fit per site):
#' * [1:2] -- `logP1` $(W min g^{-1} O_2)$ -- inverse of the slope of a photosynthesis–irradiance
#' 			curve at low light intensity
#' * [3:4] -- `logP2` $(m^2 min g^{-1} O_2)$ -- inverse maximum photosynthesis rate; can be zero 
#' 			to assume GPP is linear with light intensity instead of saturating
#' * [5:6] -- `logk600` -- coefficient of gas exchange for a gas with a Schmidt number of 600
#' * [7:8] -- `logER24_20` -- daily ecosystem respiration rate, standardized at 20 degrees C
#' * [9] -- `logsd` -- log(Error standard deviation)
#' 
#' `data` is a list of length 2 (one per site); each item is a list of (constant) data items, 
#' 				including the following:
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
twoStation_DOlogprob <- function(params, data, prior = list(), ...) {
	logprob <- sum(mapply(oneStation_DOlogprob, params = params, data = data, prior = prior))
	return(logprob)
}
