#' Simulation to develop a prior distribution for k600
#' 
#' @details This function is retained for historical purposes, it was used to derive the
#' 	variance for the informative prior for k used in the model. Parameters from Raymond et al 2012
#' 	equation 3.
#' 
#' The prior for k is defined as:
#' 
#' 		\deqn {\mu = a \times S^{b_s} \times V^{b_v} }
#' 		$\sigma = a_{se} + b_{sse} \times S + b_{vse} \times V + b_{svse} \times S \times V
#' where S is the slope (unitless, rise over run) and V is the velocity (m/s).
#' 
#' @references Raymond, P.A. et al. (2012). Scaling the gas transfer velocity and hydraulic
#' 		 geometry in streams and small rivers. Limnol. Oceanogr., 2, 41–53.
#' @return A list with parameters for the mean and standard deviation of the prior distribution
#' 	 for k.
#' @keywords internal
k_prior_sim <- function()
{
	a_mu <- 1162
	a_se <- 192/1.96
	bs_mu <- 0.77
	bs_se <- 0.028/1.96
	bv_mu <- 0.85
	bv_se <- 0.045/1.96


	# ak_mu <- 0.91 
	# ak_se <- 0.24/1.96
	# bk_mu <- 0.91
	# bk_se <- 0.036/1.96

	getK_pre <- function(S, V, a, bs, bv) {
		a * S^bs * V^bv
	}

	# getK_tru <- function(k, ak, bk) {
	# 	ak + bk * k
	# }
	VV <- seq(0.05, 1, length.out = 100)
	SS <- seq(0.0001, 0.1, length.out = 100)

	vs <- as.matrix(expand.grid(V=VV, S=SS))


	a <- rnorm(1000, a_mu, a_se)
	bs <- rnorm(1000, bs_mu, bs_se)
	bv <- rnorm(1000, bv_mu, bv_se)
	# ak <- rnorm(1000, ak_mu, ak_se)
	# bk <- rnorm(1000, bk_mu, bk_se)


	res <- matrix(NA, nrow = nrow(vs), ncol = length(a))
	for(i in 1:nrow(vs)) {
		res[i,] <- getK_pre(S = vs[i, 'S'], V = vs[i,'V'], a = a, bs = bs, bv = bv)
	}

	library(data.table)
	res2 <- data.table(S = vs[,'S'], V = vs[,'V'], kmu = rowMeans(res), ksd = apply(res, 1, sd))

	plot(log(ksd) ~ log(S), data = res2)
	plot(log(ksd) ~ log(V), data = res2)

	mod <- lm(log(ksd) ~ log(S)*log(V), data = res2)
	means <- c(a_mu, bs_mu, bv_mu)
	ses <- coefficients(mod)
	names(ses) <- c("a_se", "b_sse", "b_vse", "b_svse")

	# # compare k predicted this way to a few of Thomas' streams
	# streams <- data.frame(S=c(0.0003, 0.0013, 0.0076, 0.0019, 0.0006, 0.0045), 
	# 	depth = c(0.28, 0.30, 0.36, 0.29, 0.38, 0.18),
	# 	V = c(0.15, 0.42, 0.36, 0.27, 0.41, 0.16),
	# 	K = c(0.008, 0.0013, 0.029, 0.007, 0.002, 0.016))

	# test <- cbind(getK_pre(streams$S, streams$V, a_mu, bs_mu, bv_mu)/ (24*60), predict(mod, newdata = streams), streams$K/streams$depth)

	# plot(1:6, test[,1], ylim=c(0,max(test[,1]+1.96*test[,2], test[,3])), pch=16)
	# segments(1:6, test[,1]+1.96*test[,2], 1:6, test[,1]-1.96*test[,2])
	# points(test[,3], col='blue', pch=16)


	list(mean = means, stdev = ses, note = "Eqn for SE: log(k_se) = a_se + b_sse*log(slope) + b_vse * log(velocity) + bsvse * log(slope)*log(velocity)")
}



