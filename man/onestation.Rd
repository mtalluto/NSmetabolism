% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.r
\name{onestation}
\alias{onestation}
\title{Calibrate a single station metabolism model}
\usage{
onestation(
  oxy,
  temp,
  light,
  pressure,
  depth,
  dt,
  nsamples = 1000,
  prior = list(),
  ...
)
}
\arguments{
\item{oxy}{Dissolved oxygen time series (g/L)}

\item{temp}{Water temperature time series (degrees); one per DO observation}

\item{light}{Light, (lux or W/m^2); one per DO observation}

\item{pressure}{Air pressure (hPa); single value or one per DO observation}

\item{depth}{Water depth, in meters}

\item{dt}{Time interval of observations (minutes)}

\item{nsamples}{The number of posterior samples to return}

\item{prior}{Priors for the function chosenl, see 'details'}

\item{...}{Additional arguments to pass to the rstan}
}
\value{
A fitted rstan model
}
\description{
Calibrate a single station metabolism model
}
\details{
Implements a one-station stream metaoblism model; much code from Fuß et al 2017.

This model requires that the observations consist of a single 24-hour period. Start time
may be arbitrary, as long as an entire 24-hours is covered. No checking is done for this,
it is up to the user to ensure this is correct.

For a prior, you may specify a named list, with each item a vector of parameters to use for
the prior. Currently, the following parameters are supported:
\itemize{
\item lp1: the log of the inverse of the slope of the photosynthesis-irradiance curve; default
\code{c(0, 9)}.
\item lp2: the log of the saturation term of the photosynthesis-irradiance curve; default
\code{c(0, 9)}.
\item k: gas transfer velocity, default c(0, 10)
\item gpp: daily gross primary productivity, default c(0, 10)
\item er: daily in-situ ecosystem respiration, default c(0, 10) (note first item must be <= 0)
}

Priors all use a half normal distribution.
}
\references{
Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
dissolved organic matter and nutrient regimes. \emph{Freshw Biol} \strong{62}:582–599.
}
