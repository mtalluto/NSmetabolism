
#' Generate prior distributions for k using stream geometry
#' 
#' @details Generates means and standard deviations for half-normal priors for k600. This 
#' function uses the 7 equations from Raymond et al 2012. Which equation to use is a user
#' choice, though there are some recommendations from Raymond to guide the choice (see below,
#' and see Raymond et al 2012).
#' 
#' Rayment et al provide standard deviations for the parameters, but do NOT provide information
#' on correlations among parameters; therefore this function treats them as independent (which
#' will generally produce larger prior standard deviations). The standard deviation is chosen
#' by simulating across independent parameter combinations.
#' 
#' Slope and velocity are required for all equations. Discharge and depth are required for
#' equations 6 & 7 (discharge) and 1, 2, and 7 (depth); for equations 3-5 they can safely be
#' left set to NA.
#' 
#' Guidelines for k:
#' Best for small streams:
#' includes depth, more preditive, but less general, not in accordance with theory
#' 		eq 1, 2
#'
#' More general (includes only slope & velocity), but maybe less precise
#' 		eq 3-5
#' 
#' Equation 6 also includes discharge, predictability is lower, an alt to 3-5
#' Equation 7 includes depth and discharge; theoretically problematic, but
#' 		prediction might be better for small streams
#' 
#' @param slope Stream slope (unitless; meters of height/meters of length)
#' @param velocity In meters/second
#' @param discharge In m^3/sec
#' @param depth In meters
#' @param eqn Which equation(s) to use
#' @param nsim Number of simulations used for the standard deviation
#' 
#' @return If a single equations is given, a 2-column matrix giving the mean and standard
#'   deviation for k for each value in slope, velocity, etc. Otherwise a list of such matrices,
#'   one for each equation in eqn.
#' @references Raymond, P.A. et al. (2012). Scaling the gas transfer velocity and hydraulic
#' 		 geometry in streams and small rivers. Limnol. Oceanogr., 2, 41–53.
#' @examples
#' sl = 0.05
#' v = 0.5
#' q = 1
#' d = 0.5
#' eq1 = k_prior(sl, v, q, d, eqn = 1, nsim=3) ## normally run thousands of sims
#' @export
k_prior = function(slope, velocity, depth = NA, discharge = NA, eqn = 1:7, nsim = 5000) {
	if(!all(eqn %in% 1:7))
		stop("Unknown equation; ensure all values in eqn are in 1:7")
	if(length(eqn) > 1) {
		lapply(eqn, function(e) k_prior(slope, velocity, depth, discharge, e, nsim))
	} else {
		params = .k_parameters()
		simpars = .k_generate_pars(params$pars[[eqn]], params$pars_se[[eqn]], 5000)
		sds = .k_sd_sim(slope, velocity, depth, discharge, simpars, eqn)
		means = .k_compute(slope, velocity, depth, discharge, params$pars[[eqn]], eqn)
		cbind(mean = means, stdev = sds)
	}
}





#' generate parameter sets (assuming independence) for k sims
#' @keywords internal
.k_generate_pars = function(par, se, n) {
	mapply(function(mu, stdev, N) rnorm(N, mu, stdev), par, se, n)
}

#' generate standard dev of k for a given set of stream variables
#' @keywords internal
.k_sd_sim = function(S, V, D, Q, pars, eqn) {
	getsd = function(S, V, D, Q, p, eqn) sd(apply(p, 1, function(x)
		.k_compute(S, V, D, Q, x, eqn)))
	if(requireNamespace("parallel")) {
		parallel::mcmapply(FUN = getsd, S = S, V = V, D = D, 
			Q = Q, MoreArgs = list(p = pars, eqn = eqn))
	} else {
		mapply(FUN = getsd, S = S, V = V, D = D, 
			Q = Q, MoreArgs = list(p = pars, eqn = eqn))
	}
}

# compute k
#' @keywords internal
.k_compute = function(S, V, D, Q, pars, eq) {
	eq = paste0('eq', eq)
	# magic number is gravitational acceleration
	fr = V / sqrt(9.806*D)
	switch(eq, 
		eq1 = (V*S)^(pars[1]) * D^(pars[2]) * pars[3],
		eq2 = pars[1] * (1 - pars[2] * fr^2) * (V*S)^(pars[3]) * D^(pars[4]),
		eq3 = pars[1] * S^pars[2] * V^pars[3],
		eq4 = (V*S)^pars[1] * pars[2],
		eq5 = V*S*pars[1] + pars[2],
		eq6 = pars[1] * (V*S)^(pars[2]) * Q^pars[3],
		eq7 = pars[1] * (V*S)^(pars[2]) * Q^pars[3] * D^pars[4]
	)
}

#' Constant for k parameters
#' @keywords internal
.k_parameters = function() {
	## parameters from Raymond
	pars = list(
		eq1 = c(0.89, 0.54, 5037),
		eq2 = c(5937, 2.54, 0.89, 0.58),
		eq3 = c(1162, 0.77, 0.85),
		eq4 = c(0.76, 951.5),
		eq5 = c(2841, 2.02),
		eq6 = c(929, 0.75, 0.011),
		eq7 = c(4725, 0.86, -0.14, 0.66))
	pars_se = list(
		eq1 = c(0.020, 0.030, 604),
		eq2 = c(606, 0.223, 0.017, 0.027),
		eq3 = c(192, 0.028, 0.045),
		eq4 = c(0.027, 144),
		eq5 = c(107, 0.209),
		eq6 = c(141, 0.027, 0.016),
		eq7 = c(445, 0.016, 0.012, 0.029))
	list(pars = pars, pars_se = pars_se)
}
