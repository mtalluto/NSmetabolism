#' Compute Gross Primary Productivity
#' @param par light intensity ($W/m^2$)
#' @param P1 $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
#' 		 light intensity
#' @param P2 $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
#' 		  GPP is linear with light intensity instead of saturating
#' @references Uehlinger, U., et al. 2000. Variability of photosynthesis‐irradiance curves and 
#' 		 ecosystem respiration in a small river. *Freshw Biol* **44**:493–507.
#' @return Photosynthetic rate
#' @export
gpp <- function(par, P1, P2) {
	# Uehlinger et al 2000 eq 3b
	par / (P1 + P2 * par)
}




#' Compute pressure (in hPa) from elevation
#' @param elev Vector of elevations
#' @return Pressure for given elevations in hPa
#' @export
pressureFromElevation <- function(elev) {
	# barometric eqn
	g <- 9.80665 # acceleration of gravity
	M <- 0.0289644 # molar mass of air
	Lb <- -0.0065 # temperature lapse rate
	Pb <- 101325 # static pressure
	Tb <- 288.15 # standard temperature
	Rs <- 8.31432 # gas constant
	presPascals <- Pb * (Tb / (Tb + Lb * elev))^(g * M / (Rs * Lb))
	presPascals / 100
}

#' Convert hPa to atmospheres
#' @param P Vector of pressures (in hPa)
#' @return Pressure in atm
#' @export
hPaToAtm <- function(P) {
	P * 9.86923266e-4
}


