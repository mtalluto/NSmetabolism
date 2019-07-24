




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


