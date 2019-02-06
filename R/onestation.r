#' Time derivative of DO concentration
#' 
#' @param t Time step for the ode integration
#' @param y Current state (dissolved oxygen concentration)
#' @param parms a list of parameters; see details
#' @param data a list of data for the ode; see details
#' @details `params` must be a list including the following named items:
#' * `P1` $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
#' 		 light intensity
#' * `P2` $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
#' 		  GPP is linear with light intensity instead of saturating
#' * `k600` coefficient of gas exchange for a gas with a Schmidt number of 600; not supported,
#' 		currently fixed to 2/(24*60) (see metabolism_functions.r/kT())
#' * `ER24_20` daily ecosystem respiration rate, standardized at 20 degrees C
#' 
#' `data` is a named list of (constant) data items, including the following:
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
dDOdt <- function(t, y, parms, data)
{
	# gross primary production; Fuß et al eq 2
	# Uehlinger et al 2000 eq 3b
	GPP <- (data$PAR(t) / (parms$P1 + parms$P2 * data$PAR(t)))

	# rearation flux, Van de Bogert et al eqn 4; note that k600 is fixed
	RF <- kT(data$temp(t)) * (osat(data$temp(t), data$P) - y)

	# ecosystem respiration; Fuß et al eq 4
	ER <- (parms$ER24_20 / (60*24)) * (1.045^data$temp(t) - 20)

	(GPP + ER + RF) / data$z
}

#' Predict a time series of dissolved oxygen concentration
#' @param initial initial dissolved oxygen concentration
#' @param times The times at which to solve the system
#' @param params a list of parameters; see details
#' @param data a list of data for the ode; see details
#' @param dt Time step for integration; only used if `method == 'euler'`
#' @param method The integration method to use; default is euler

#' @details Light and temperature time series will be approximated using linear interpolation
#' 		at the desired time steps
#' 
#' `params` must be a list including the following named items:
#' * `P1` $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
#' 		 light intensity
#' * `P2` $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
#' 		  GPP is linear with light intensity instead of saturating
#' * `k600` coefficient of gas exchange for a gas with a Schmidt number of 600; not supported,
#' 		currently fixed to 2/(24*60) (see metabolism_functions.r/kT())
#' * `ER24_20` daily ecosystem respiration rate, standardized at 20 degrees C
#' 
#' `data` is a named list of (constant) data items, including the following:
#' * `PAR`  2-column data frame; first column is light, second is time of observation
#' * `temp` 2-column data frame; first column is temperature, second is time of observation
#' * `P`    pressure, in atmospheres
#' * `z`    Depth, in meters
#' @return Time series of dissolved oxygen concentrations
#' @export
DO_predict <- function(initial, times, params, data, dt = 1, method=c('euler', 'lsoda')) {
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
				statet <- statet + dt * dDOdt(t, statet, params, dDOData)
				t <- t + dt
			}
			state[i] <- statet
		}
	} else {
		state <- deSolve::ode(initial, times = times, funct = dDOdt, parms = params, 
			data = dDOData)
	}
	return(state)
}

#' Log posterior probability for the DO curve given empirical data
DO_logprob <- function() {

}

