% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onestation.r
\name{oneStation_DOPredict}
\alias{oneStation_DOPredict}
\title{Predict a time series of dissolved oxygen concentration}
\usage{
oneStation_DOPredict(
  initial,
  times,
  params,
  data,
  dt = 1,
  method = c("euler", "lsoda"),
  gpp = FALSE
)
}
\arguments{
\item{initial}{initial dissolved oxygen concentration}

\item{times}{The times at which to solve the system}

\item{params}{a named vector of parameters; see details}

\item{data}{a list of data for the ode; see details}

\item{dt}{Time step for integration; only used if \code{method == 'euler'}}

\item{method}{The integration method to use; default is euler}

\item{gpp}{Logical, should GPP at each time step be returned as well?}
}
\value{
Time series of dissolved oxygen concentrations
}
\description{
Predict a time series of dissolved oxygen concentration
}
\details{
Light and temperature time series will be approximated using linear interpolation
at the desired time steps

\code{params} must be a vector including the following named items:
\itemize{
\item \code{P1} $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
light intensity
\item \code{P2} $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
GPP is linear with light intensity instead of saturating
\item \code{k600} coefficient of gas exchange for a gas with a Schmidt number of 600; (see \code{\link[=kT]{kT()}})
\item \code{ER24_20} daily ecosystem respiration rate, standardized at 20 degrees C
}

\code{data} is a named list of (constant) data items, including the following:
\itemize{
\item \code{PAR}  2-column data frame; first column is light, second is time of observation
\item \code{temp} 2-column data frame; first column is temperature, second is time of observation
\item \code{P}    pressure, in atmospheres
\item \code{z}    Depth, in meters
}
}
