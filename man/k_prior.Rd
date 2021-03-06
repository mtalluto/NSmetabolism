% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.r
\name{k_prior}
\alias{k_prior}
\title{Generate prior distributions for k using stream geometry}
\usage{
k_prior(slope, velocity, depth = NA, discharge = NA, eqn = 1:7, nsim = 5000)
}
\arguments{
\item{slope}{Stream slope (unitless; meters of height/meters of length)}

\item{velocity}{In meters/second}

\item{depth}{In meters}

\item{discharge}{In m^3/sec}

\item{eqn}{Which equation(s) to use}

\item{nsim}{Number of simulations used for the standard deviation}
}
\value{
If a single equations is given, a 2-column matrix giving the mean and standard
deviation for k for each value in slope, velocity, etc. Otherwise a list of such matrices,
one for each equation in eqn. Units for k are meters/day
}
\description{
Generate prior distributions for k using stream geometry
}
\details{
Generates means and standard deviations for half-normal priors for k600. This
function uses the 7 equations from Raymond et al 2012. Which equation to use is a user
choice, though there are some recommendations from Raymond to guide the choice (see below,
and see Raymond et al 2012).

Rayment et al provide standard deviations for the parameters, but do NOT provide information
on correlations among parameters; therefore this function treats them as independent (which
will generally produce larger prior standard deviations). The standard deviation is chosen
by simulating across independent parameter combinations.

Slope and velocity are required for all equations. Discharge and depth are required for
equations 6 & 7 (discharge) and 1, 2, and 7 (depth); for equations 3-5 they can safely be
left set to NA.

Guidelines for k:
Best for small streams:
includes depth, more preditive, but less general, not in accordance with theory
eq 1, 2

More general (includes only slope & velocity), but maybe less precise
eq 3-5

Equation 6 also includes discharge, predictability is lower, an alt to 3-5
Equation 7 includes depth and discharge; theoretically problematic, but
prediction might be better for small streams
}
\examples{
sl = 0.05
v = 0.5
q = 1
d = 0.5
eq1 = k_prior(sl, v, q, d, eqn = 1, nsim=3) ## normally run thousands of sims
}
\references{
Raymond, P.A. et al. (2012). Scaling the gas transfer velocity and hydraulic
geometry in streams and small rivers. Limnol. Oceanogr., 2, 41–53.
}
