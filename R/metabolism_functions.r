#' Compute dissolved oxygen saturation given temperature and pressure
#' @param temp Water temperature, in degrees C
#' @param P pressure, in atmospheres
#' @references Benson BB and Krause D, Jr. 1980. The concentration and isotopic fractionation of
#' 		 gases dissolved in freshwater in equilibrium with the atmosphere. 1. Oxygen. *Limnol.
#' 		 Oceanogr.* **25**(4):662-671. \\
#' 		 Benson, B.B. and Krause, D. 1984. The concentration and isotopic fractionation of 
#' 		 oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. 
#' 		 *Limnol. Oceanogr.*, **29**:620–632.
#' @return vector of DO saturation values
#' @keywords internal
osat <- function(temp, P) {	
	tempK <- temp + 273.15

	# C*o, the unit standard atmospheric concentration by volume for oxygen; in mg/kg
	# eqn 31 from Benson and Krause 1984.
	Cstaro <- exp(1.3874202e2 + (1.572288e5 / tempK) - (6.637149e7 / tempK^2) + 
				(1.243678e10 / tempK^3) - (8.621061e11 / tempK^4))

	# original equation 25 from Benson and Krause 1980
	# Cstaro <- exp(-1.3934411e2 - (1.575701e5 / tempK) - (6.642308e7 / tempK^2) - 
	#			(1.243800e10 / tempK^3) - (8.621949e11 / tempK^4))

	# eqn 13, 
	# the negative of the second pressure coefficient in the virial expansion for the real gas 
	# behavior of oxygen.
	theta <- 0.000975 - 1.426e-5 * temp + 6.436e-8 * temp^2

	# eqn 23
	# the saturated vapor pressure of water in atmospheres at the temperature of equilibrium.
	# in atmospheres
	Pwv <- exp(11.8571 - (3840.7 / tempK) - (216961 / tempK^2))

	# eqn 28
	Cstaro * ((P - Pwv) * (1 - theta * P)) / ((1 - Pwv) * (1 - theta))
}

#' Compute the temperature-corrected reaction coefficient
#' @param temp Temperature in degrees C
#' @param K600 gas exchange coefficient
#' @references Van de Bogert MC et al. 2007. Assessing pelagic and benthic metabolism using free
#' 		 water measurements. *Limnol. Oceanogr.: Methods* **5**:145–155. \\
#' 		 Wanninkhof, R. 1992. Relationship between wind speed and gas exchange over the ocean. 
#' 		 *Journal of Geophysical Research* **97**:7373.
#' @return vector of kT values
#' @keywords internal
kT <- function(temp, k600 = 2/(24*60)) {
	# compute Schmidt number for oxygen
	# parameters from Wanninkhof 1992. appendix
	Sc <- 1800.6 - 120.10 * temp + 3.7818 * temp^2 - 0.047608 * temp^3
	
	# NOTE - original params from Thomas' code below
	# Sc <- 1568 - 86.04*temp + 2.142*temp^2 - 0.0216*temp^3
	
	# Van de Bogert et al eqn 5
	k600 * (Sc / 600)^-0.5
}


#' Compute pressure in atmospheres elevation
#' @param elev Vector of elevations
#' @return Pressure for given elevations
#' @export
pressureFromElevation <- function(elev) {
	# barometric eqn
	g <- 9.80665 # acceleration of gravity
	M <- 0.0289644 # molar mass of air
	Lb <- -0.0065 # temperature lapse rate
	Pb <- 101325 # static pressure
	Tb <- 288.15 # standard temperature
	Rs <- 8.31432 # gas constant
	presPascals <- Pb * (Tb / (Tb - Lb * elevation))^(g * M / (Rs * Lb))
	presPascals * 9.86923266*10^-6
}